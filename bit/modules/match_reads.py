"""
Core logic for `bit match-reads`.

By default, each input fasta file is one reference made up of all its sequences.
With `per_seq=True`, each sequence is its own reference. Reads are matched against
individual sequences, and those hits are then collapsed to references, so a read
matching two sequences of the same reference is still assigned 'uniquely' to that reference.

A read matches a sequence if the read, or its reverse complement, is an exact
substring of that sequence. With `circular=True`, sequences are treated as circles
(e.g., plasmids), so reads spanning the origin are matched too.

A canonical k-mer prefilter narrows the candidate sequences for each read before
the exact check. Every k-mer of an exactly matching read must also be a k-mer of
that sequence, so a sequence is excluded only if it lacks a sampled k-mer, and the
prefilter can't cause a true exact match to be missed.

With `max_edits` > 0, reads are instead assigned to the reference(s) they align to
with the lowest edit distance (substitutions + indels), if within `max_edits`, where
a reference's distance is the lowest among its sequences. Ties are ambiguous.
That mode uses its own prefilter, which likewise can't cause a read within
`max_edits` of a reference to be missed.

The matching, reading and writing run in bit's compiled core (rust/src/match_reads.rs, the
ReadMatcher in rust/src/lib.rs); this module loads the references, runs it, and writes the
summaries.

For paired-end input, a pair is assigned to the references both mates match (and
with edits allowed, ranked by the mates' combined edit distance).
"""

import gzip
import io
import os
import sys
from collections import Counter, defaultdict
from contextlib import nullcontext
from pathlib import Path

from tqdm import tqdm # type: ignore

from bit.modules.general import (check_if_output_dir_exists, is_gzipped, log_command_run,
                                 notify_premature_exit, report_message, spinner)


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# spacing of the sparse positional index of the compiled matcher (rust/src/match_reads.rs)
POSITION_STEP = 16

# fixed settings (not exposed as parameters)
KMER_SIZE = 31
NUM_KMER_SAMPLES = 25   # k-mers sampled per read by the exact-match prefilter
SORT_BUFFER_ROWS = 500_000  # read-hits rows held in memory before a sorted chunk goes to disk

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".fsa", ".ffn", ".fsta")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def plural(n, word):
    """Count with the word made plural if needed, e.g. '1 file', '8 files', '1,024 sequences'."""
    return f"{n:,} {word if n == 1 else word + 's'}"


### refs and k-mer index ###

def ref_name_from_path(path):
    """File name minus any .gz and fasta extension, e.g. 'genome-1.fasta.gz' -> 'genome-1'."""
    name = os.path.basename(str(path))
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in FASTA_EXTENSIONS:
        if name.endswith(ext) and len(name) > len(ext):
            return name[:-len(ext)]
    return name


def load_refs(ref_paths, per_seq=False):
    """
    Returns (refs, seqs):
        refs: [(reference name, source file basename), ...]
        seqs: [(sequence name, reference index, uppercase sequence), ...]

    Reads are matched against sequences, and their hits are then collapsed to references.
    By default, each input file is one reference made up of all of its sequences (e.g., a
    genome's chromosome and plasmids), named by its file name minus fasta/gzip extensions.
    With per_seq, each sequence is its own reference, named by the sequence name.
    """
    refs, seqs = [], []
    seen_refs = {}

    def add_ref(name, path):
        if name in seen_refs:
            kind = "sequence name" if per_seq else "reference name (from the file name)"
            print(f"\n    Duplicate {kind} '{name}' from '{path}' (also from '{seen_refs[name]}').")
            print("    Reference names need to be unique.")
            notify_premature_exit()
        seen_refs[name] = path
        refs.append((name, os.path.basename(str(path))))
        return len(refs) - 1

    for path in ref_paths:
        file_ref_idx = None if per_seq else add_ref(ref_name_from_path(path), path)
        seen_seqs = set()
        reader = FastxReader(path)
        try:
            for name, seq, _ in reader:
                if not seq:
                    print(f"\n    Sequence '{name}' in '{path}' has no sequence.")
                    notify_premature_exit()
                if name in seen_seqs:
                    print(f"\n    Duplicate sequence name '{name}' in '{path}'.")
                    notify_premature_exit()
                seen_seqs.add(name)
                ref_idx = add_ref(name, path) if per_seq else file_ref_idx
                seqs.append((name, ref_idx, seq.upper()))
        finally:
            reader.close()

        if not seen_seqs:
            print(f"\n    No sequences were found in '{path}'.")
            notify_premature_exit()

    return refs, seqs


def find_identical_seqs(seqs, circular):
    """
    Returns index pairs of sequences that are identical (in either orientation, and up
    to rotation if circular).
    """
    by_len = defaultdict(list)
    for i, (_, _, seq) in enumerate(seqs):
        by_len[len(seq)].append(i)

    identical = []
    for idxs in by_len.values():
        for pos, a in enumerate(idxs):
            a_seq = seqs[a][2]
            a_search = a_seq * 2 if circular else a_seq
            for b in idxs[pos + 1:]:
                b_seq = seqs[b][2]
                if b_seq in a_search or revcomp(b_seq) in a_search:
                    identical.append((a, b))
    return identical


class FastxReader:
    """
    Iterates (name, seq, qual or None) from fasta or fastq, gzipped or not, and reports
    how many bytes of the file on disk have been consumed (for progress bars).
    """

    def __init__(self, path):
        self.path = str(path)
        self.total_bytes = os.path.getsize(self.path)
        self._raw = open(self.path, "rb")
        stream = gzip.GzipFile(fileobj=self._raw) if is_gzipped(self.path) else self._raw
        self._handle = io.TextIOWrapper(stream)

    def bytes_read(self):
        return self._raw.tell()

    def close(self):
        self._handle.close()
        self._raw.close()

    def __iter__(self):
        f = self._handle
        line = f.readline()
        while line and not line.strip():
            line = f.readline()
        if not line:
            return

        if line.startswith(">"):
            name, chunks = _parse_name(line), []
            for line in f:
                if line.startswith(">"):
                    yield name, "".join(chunks), None
                    name, chunks = _parse_name(line), []
                else:
                    chunks.append(line.strip())
            yield name, "".join(chunks), None

        elif line.startswith("@"):
            while line:
                if line.strip():
                    if not line.startswith("@"):
                        self._bad_fastq(line)
                    name = _parse_name(line)
                    seq = f.readline().strip()
                    plus = f.readline()
                    qual = f.readline().strip()
                    if not plus.startswith("+") or len(qual) != len(seq):
                        self._bad_fastq(line)
                    yield name, seq, qual
                line = f.readline()

        else:
            print(f"\n    '{self.path}' doesn't look like fasta or fastq format.")
            notify_premature_exit()

    def _bad_fastq(self, line):
        print(f"\n    Unexpected fastq format in '{self.path}' near: {line.strip()[:60]}")
        print("    (multi-line fastq records aren't supported)")
        notify_premature_exit()


def _parse_name(header_line):
    parts = header_line[1:].split()
    return parts[0] if parts else ""


def output_paths(output_dir, output_prefix=""):
    """Output locations within `output_dir`, with `output_prefix` prepended to each name."""
    out = Path(output_dir)
    return {
        "hits": str(out / f"{output_prefix}read-hits.tsv"),
        "summary": str(out / f"{output_prefix}summary.txt"),
        "ref_summary": str(out / f"{output_prefix}ref-summary.tsv"),
        "seq_summary": str(out / f"{output_prefix}seq-summary.tsv"),
        "reads_dir": str(out / f"{output_prefix}reads"),
        "log": str(out / f"{output_prefix}command-execution-info.txt"),
    }


def summary_count_lines(unit, total, unique, ambiguous):
    """
    The overall counts, as printed to the terminal and written to summary.txt, e.g.:
        Total reads:                      1,000
        Uniquely assigned:                800 (80.00%)
        Ambiguous (multiple refs):        50 (5.00%)
    """
    pct = lambda x: f"{100 * x / total:.2f}%" if total else "NA"
    return [
        f"{'Total ' + unit + ':':<34}{total:,}",
        f"{'Uniquely assigned:':<34}{unique:,} ({pct(unique)})",
        f"{'Ambiguous (multiple refs):':<34}{ambiguous:,} ({pct(ambiguous)})",
    ]


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed, inputs=()):
    """
    Stops if `output_dir` exists, unless forced: then this program's own earlier outputs there
    are removed (so, e.g., reads written with --write-reads can't be left beside new ones) and
    anything else in the directory is kept. Then creates it and logs the command run there.
    """
    paths = output_paths(output_dir, output_prefix)
    check_if_output_dir_exists(output_dir, force_overwrite, [os.path.basename(p) for p in paths.values()],
                               os.path.basename(paths["log"]), inputs)
    os.makedirs(output_dir, exist_ok=True)
    log_command_run(full_cmd_executed, output_dir, paths["log"])


def read_file_base(ref_name, mate, paired):
    """ a written-reads file's name, less its extension (.fastq.gz or .fasta.gz) """
    safe = ref_name.replace("/", "_")
    return f"{safe}_R{mate}" if paired else safe


def match_reads(ref_paths, read_1, read_2=None, output_dir="match-reads", output_prefix="", per_seq=False,
                 circular=False, max_edits=0, min_frac_of_seq=0.0,
                 write_reads=False, jobs=1, show_progress=True):
    """
    Runs the full assignment and writes outputs. Returns a dict of summary counts.
    """

    paired = read_2 is not None
    paths = output_paths(output_dir, output_prefix)
    os.makedirs(output_dir, exist_ok=True)
    k = KMER_SIZE

    # loading refs and building the indexes
    refs, seqs = load_refs(ref_paths, per_seq)
    seq_to_ref = [ri for _, ri, _ in seqs]

    print(f"\n    Loaded {plural(len(refs), 'reference')} ({plural(len(seqs), 'sequence')}) from "
          f"{plural(len(ref_paths), 'file')}, treated as {'circular' if circular else 'linear'}")
    if max_edits > 0:
        print(f"    Allowing up to {plural(max_edits, 'edit')} per {'mate' if paired else 'read'}")

    if not per_seq and len(ref_paths) == 1 and len(seqs) > 1:
        report_message(f"Note: the single reference file holds {len(seqs):,} sequences, which are all "
                       "being treated as one reference. Add `--per-seq` if each sequence should be "
                       "its own reference (e.g., a mix of different plasmids in one fasta).",
                       initial_indent="    ", subsequent_indent="    ")

    def seq_label(si):
        name, ri, _ = seqs[si]
        return name if per_seq else f"{refs[ri][0]}:{name}"

    for a, b in find_identical_seqs(seqs, circular):
        if seq_to_ref[a] == seq_to_ref[b]:
            # identical sequences within one reference don't affect assignment
            continue
        report_message(f"Note: sequences '{seq_label(a)}' and '{seq_label(b)}' are identical"
                       f"{' (up to rotation/strand)' if circular else ' (up to strand)'}, "
                       "so reads matching them can't be assigned uniquely.",
                       initial_indent="    ", subsequent_indent="    ")

    short_seqs = [si for si, (_, _, seq) in enumerate(seqs) if len(seq) < k]
    if short_seqs:
        report_message(f"Note: {plural(len(short_seqs), 'sequence')} "
                       f"{'is' if len(short_seqs) == 1 else 'are'} shorter than {k} bases and can "
                       f"only be matched by reads shorter than that (e.g., '{seq_label(short_seqs[0])}').",
                       initial_indent="    ", subsequent_indent="    ")

    not_sequence = next((si for si, (_, _, seq) in enumerate(seqs) if not seq.isascii()), None)
    if not_sequence is not None:
        print(f"\n    The reference sequence '{seq_label(not_sequence)}' holds characters that aren't a nucleotide sequence's.")
        notify_premature_exit()

    print()
    building = spinner("Building indexes...", "Built indexes ") if show_progress else nullcontext()
    with building:
        from bit._native import ReadMatcher
        native = ReadMatcher([seq for _, _, seq in seqs], seq_to_ref, k, circular, max_edits,
                             NUM_KMER_SAMPLES, POSITION_STEP, jobs)

    unit = "pairs" if paired else "reads"
    n_total = n_unique = n_ambig = 0
    unique_counts, ambig_counts = Counter(), Counter()
    seq_unique_counts, seq_ambig_counts = Counter(), Counter()

    header = ["read_id", "read_lengths" if paired else "read_length", "status", "num_hits",
              "matching_refs", "edit_distance", "next_best_edit_distance"]
    if not per_seq:
        header.append("matching_seqs")
    # (the compiled core writes the table, its rows sorted by read length (combined for pairs), longest first, then by
    # edit distance)
    read_paths = [read_1] + ([read_2] if paired else [])
    total_bytes = sum(os.path.getsize(p) for p in read_paths)
    print(f"\n    Assigning {unit} with {plural(jobs, 'job')}...\n")
    pbar = tqdm(total=total_bytes, unit="B", unit_scale=True, unit_divisor=1024,
                ncols=80, disable=not show_progress, file=sys.stdout,
                bar_format="    {l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]")
    pbar.set_postfix_str(f"0 {unit}")
    last_pos = 0

    def on_progress(n_in_batch, pos):
        """ after each batch the compiled core matches: the progress bar """
        nonlocal n_total, last_pos
        n_total += n_in_batch
        pbar.update(min(pos, total_bytes) - last_pos)
        last_pos = min(pos, total_bytes)
        pbar.set_postfix_str(f"{n_total:,} {unit}", refresh=False)

    try:
        # reading, matching and writing the reads all in the compiled core, across its threads
        reads_out = None
        if write_reads:
            os.makedirs(paths["reads_dir"], exist_ok=True)
            reads_out = [[os.path.join(paths["reads_dir"], read_file_base(name, mate, paired))
                          for mate in ((1, 2) if paired else (None,))] for name, _ in refs]
        try:
            counts = native.assign_files(
                read_1, read_2, min_frac_of_seq, paths["hits"], "\t".join(header),
                [name for name, _ in refs], None if per_seq else [seq_label(si) for si in range(len(seqs))],
                SORT_BUFFER_ROWS, on_progress, reads_out)
        except ValueError as e:  # a problem with the read files
            print(f"\n\n    {e}")
            notify_premature_exit()
        n_total, n_unique, n_ambig = counts[:3]
        for counter, values in zip((unique_counts, ambig_counts, seq_unique_counts, seq_ambig_counts), counts[3:]):
            counter.update({i: v for i, v in enumerate(values) if v})
        # finishing the bar cleanly
        pbar.update(total_bytes - last_pos)
    finally:
        pbar.close()

    ref_num_seqs, ref_lengths = Counter(), Counter()
    for _, ri, seq in seqs:
        ref_num_seqs[ri] += 1
        ref_lengths[ri] += len(seq)

    with open(paths["ref_summary"], "w") as out:
        out.write(f"ref\tsource_file\tnum_seqs\ttotal_length\tunique_{unit}\tambiguous_{unit}\n")
        for i, (name, src) in enumerate(refs):
            out.write(f"{name}\t{src}\t{ref_num_seqs[i]}\t{ref_lengths[i]}\t"
                      f"{unique_counts[i]}\t{ambig_counts[i]}\n")

    # the per-sequence summary only adds anything when some reference has multiple sequences
    # (never the case with per_seq, where every reference is one sequence)
    write_seq_summary = any(n > 1 for n in ref_num_seqs.values())
    if write_seq_summary:
        with open(paths["seq_summary"], "w") as out:
            out.write(f"ref\tseq\tseq_length\tunique_{unit}\tambiguous_{unit}\n")
            for si, (name, ri, seq) in enumerate(seqs):
                out.write(f"{refs[ri][0]}\t{name}\t{len(seq)}\t"
                          f"{seq_unique_counts[si]}\t{seq_ambig_counts[si]}\n")

    with open(paths["summary"], "w") as out:
        for line in summary_count_lines(unit, n_total, n_unique, n_ambig):
            out.write(line + "\n")

    return {
        "unit": unit,
        "total": n_total,
        "unique": n_unique,
        "ambiguous": n_ambig,
        "paths": paths,
        "wrote_seq_summary": write_seq_summary,
        "write_reads": write_reads,
    }
