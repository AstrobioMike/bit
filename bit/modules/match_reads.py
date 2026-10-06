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
That mode uses its own prefilter (see EditMatcher), which likewise
can't cause a read within `max_edits` of a reference to be missed.

For paired-end input, a pair is assigned to the references both mates match (and
with edits allowed, ranked by the mates' combined edit distance).
"""

import gzip
import heapq
import io
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import nullcontext
from itertools import chain, islice
from operator import itemgetter
from pathlib import Path

import edlib # type: ignore
from tqdm import tqdm # type: ignore

from bit.modules.general import (check_if_output_dir_exists, is_gzipped, log_command_run,
                                 notify_premature_exit, report_message, spinner)


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# spacing of the sparse positional index (see build_position_index())
POSITION_STEP = 16

# fixed settings (not exposed as parameters)
KMER_SIZE = 31
NUM_KMER_SAMPLES = 25   # k-mers sampled per read by the exact-match prefilter
BATCH_SIZE = 500        # reads (or pairs) per batch sent to each worker
SORT_BUFFER_ROWS = 500_000  # read-hits rows held in memory before a sorted chunk goes to disk

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".fsa", ".ffn", ".fsta")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def plural(n, word):
    """Count with the word made plural if needed, e.g. '1 file', '8 files', '1,024 sequences'."""
    return f"{n:,} {word if n == 1 else word + 's'}"


### refs and k-mer index ###

def canon(kmer):
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


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


def build_kmer_index(refs, k, circular):
    """
    One combined index: canonical k-mer -> frozenset of ref indices containing it.
    Identical sets are shared, which keeps memory low when refs are highly similar.
    """
    index = defaultdict(set)
    for i, (_, _, seq) in enumerate(refs):
        L = len(seq)
        if circular:
            # one full cycle of k-mers, including those spanning the origin
            circ = (seq * (k // L + 2))[: L + k - 1]
            n_kmers = L
        else:
            circ = seq
            n_kmers = L - k + 1
        for j in range(max(n_kmers, 0)):
            index[canon(circ[j:j + k])].add(i)

    shared = {}
    return {km: shared.setdefault(frozenset(s), frozenset(s)) for km, s in index.items()}


def get_candidates(seq, index, k, n_samples):
    """
    Starts with every ref as a candidate and excludes any ref missing one of the
    k-mers sampled evenly across the read (including both ends). Returns the refs left.
    """
    last = len(seq) - k
    if n_samples > 1:
        positions = {round(j * last / (n_samples - 1)) for j in range(n_samples)}
    else:
        positions = {0}

    cands = None
    for p in positions:
        hit = index.get(canon(seq[p:p + k]))
        if not hit:
            return frozenset()
        cands = hit if cands is None else cands & hit
        if not cands:
            return cands
    return cands


def build_position_index(refs, k, circular, step=POSITION_STEP):
    """
    Sparse positional index: forward k-mer -> [(ref index, position), ...] for k-mers
    starting every `step` bases along each ref, with no gap between indexed positions
    larger than `step` (including across the origin if circular, and the ref's end if
    linear). So any `step` consecutive k-mers of an exactly matching read include at
    least one indexed k-mer, which gives the read's exact offset in that ref.
    """
    pos_index = defaultdict(list)
    for i, (_, _, seq) in enumerate(refs):
        L = len(seq)
        if L < k:
            continue
        if circular:
            circ = seq + seq[:k - 1]
            positions = range(0, L, step)
        else:
            circ = seq
            last = L - k
            positions = list(range(0, last + 1, step))
            if positions[-1] != last:
                positions.append(last)
        for p in positions:
            pos_index[circ[p:p + k]].append((i, p))
    return dict(pos_index)


class _RefMatcher:
    """
    Shared base for the matchers: holds the refs (once, as-is) and comparisons that
    wrap from a circular ref's end back to its start, which is equivalent to comparing
    against the ref concatenated with itself.
    """

    def __init__(self, refs, circular):
        self.circular = circular
        self.ref_seqs = [seq for _, _, seq in refs]

    def _matches_at(self, ref_seq, read, start):
        """True if `read` exactly matches `ref_seq` beginning at `start`."""
        if not self.circular:
            return ref_seq.startswith(read, start)

        L, n = len(ref_seq), len(read)
        to_end = L - start
        if n <= to_end:
            return ref_seq.startswith(read, start)

        # read runs past the ref's end: compare up to the end, then wrap to the start
        if not ref_seq.startswith(read[:to_end], start):
            return False
        pos = to_end
        # full passes around the circle, for reads longer than the ref (e.g., multimers)
        while n - pos >= L:
            if not read.startswith(ref_seq, pos):
                return False
            pos += L
        return pos == n or ref_seq.startswith(read[pos:])

    def _contains(self, ref_seq, read):
        """True if `read` occurs anywhere in `ref_seq` (wrapping around if circular)."""
        if read in ref_seq:
            return True
        if not self.circular:
            return False

        L, n = len(ref_seq), len(read)
        if n <= L:
            # a match not found above must span the origin, so only the junction needs checking
            junction = ref_seq[L - (n - 1):] + ref_seq[:n - 1]
            return read in junction
        return read in ref_seq * (n // L + 2)


class ExactMatcher(_RefMatcher):
    """
    Finds which refs a read (or its reverse complement) exactly matches.

    1. the canonical k-mer prefilter excludes refs missing any sampled read k-mer
    2. for the remaining candidates, the read's first `step` k-mers (each strand) are
       looked up in the sparse positional index, giving exact offsets to compare the
       full read at; that's one direct comparison per candidate rather than a scan
       of each ref

    Reads too short for step 2 (fewer than `step` k-mers) fall back to substring scans.
    """

    def __init__(self, refs, k, n_samples, circular, step=POSITION_STEP):
        super().__init__(refs, circular)
        self.k = k
        self.n_samples = n_samples
        self.step = step
        self.index = build_kmer_index(refs, k, circular)
        self.pos_index = build_position_index(refs, k, circular, step)

    def distances(self, seq, min_frac_of_seq=0.0):
        """{ref index: edit distance} for refs `seq` matches; always 0 here (exact only)."""
        return {i: 0 for i in self.hits(seq, min_frac_of_seq)}

    def hits(self, seq, min_frac_of_seq=0.0):
        """Returns a sorted list of indices of refs that `seq` exactly matches."""
        L = len(seq)
        if L == 0:
            return []

        k = self.k
        s = seq.upper()

        if L >= k:
            cands = get_candidates(s, self.index, k, self.n_samples)
        else:
            cands = range(len(self.ref_seqs))

        cands = {i for i in cands
                 if L >= min_frac_of_seq * len(self.ref_seqs[i])
                 and (self.circular or L <= len(self.ref_seqs[i]))}
        if not cands:
            return []

        s_rc = revcomp(s)

        # short reads: plain substring checks
        if L - k + 1 < self.step:
            return sorted(i for i in cands
                          if self._contains(self.ref_seqs[i], s)
                          or self._contains(self.ref_seqs[i], s_rc))

        # longer reads: anchored comparisons from the sparse positional index
        hits = set()
        for strand in (s, s_rc):
            tried = set()
            for t in range(self.step):
                entries = self.pos_index.get(strand[t:t + k])
                if not entries:
                    continue
                for i, pos in entries:
                    if i in hits or i not in cands:
                        continue
                    start = pos - t
                    if self.circular:
                        start %= len(self.ref_seqs[i])
                    elif start < 0:
                        continue
                    if (i, start) in tried:
                        continue
                    tried.add((i, start))
                    if self._matches_at(self.ref_seqs[i], strand, start):
                        hits.add(i)
            if len(hits) == len(cands):
                break

        return sorted(hits)


def build_dense_position_index(refs, k, circular):
    """Forward k-mer -> [(ref index, position), ...] for every k-mer position in every ref."""
    pos_index = defaultdict(list)
    for i, (_, _, seq) in enumerate(refs):
        L = len(seq)
        if L < k:
            continue
        if circular:
            circ = seq + seq[:k - 1]
            n_kmers = L
        else:
            circ = seq
            n_kmers = L - k + 1
        for p in range(n_kmers):
            pos_index[circ[p:p + k]].append((i, p))
    return dict(pos_index)


class EditMatcher(_RefMatcher):
    """
    Finds the refs a read (or its reverse complement) aligns to within `max_edits` edits
    (substitutions + indels), and the edit distance to each.

    1. prefilter: the read is split into non-overlapping k-mers. Each edit can break at
       most one of them, so the number missing from a ref is a lower bound on the read's
       edit distance to it. Refs (per strand) missing more than `max_edits` are excluded
       without aligning, so no read within `max_edits` is missed.
    2. anchoring: any `max_edits + 1` of those k-mers include at least one left intact
       by the true alignment, so the first `max_edits + 1` are looked up in a dense
       positional index. An intact k-mer gives the alignment's diagonal, and the
       alignment's start in the ref must be within +/- max_edits of it.
    3. distance: an exact comparison is tried first at each diagonal. Otherwise, edlib
       aligns the read in prefix mode (start fixed, end free) from each candidate start,
       which lets edlib's banding with the `k` cutoff keep it fast. The cutoff tightens
       as better alignments are found.

    Reads with no more than `max_edits` non-overlapping k-mers fall back to full infix
    alignment against each eligible ref.
    """

    def __init__(self, refs, k, circular, max_edits):
        super().__init__(refs, circular)
        self.k = k
        self.max_edits = max_edits
        self.dense_index = build_dense_position_index(refs, k, circular)
        # k-mer -> refs containing it, with identical sets shared to keep memory low
        shared = {}
        self.kmer_refs = {}
        for km, entries in self.dense_index.items():
            ids = frozenset(i for i, _ in entries)
            self.kmer_refs[km] = shared.setdefault(ids, ids)

    def _eligible(self, i, n, min_frac_of_seq):
        L = len(self.ref_seqs[i])
        return n >= min_frac_of_seq * L and (self.circular or n <= L + self.max_edits)

    def _window(self, i, start, length):
        """Ref sequence of `length` beginning at `start` (wrapping if circular)."""
        ref = self.ref_seqs[i]
        L = len(ref)
        if not self.circular:
            return ref[start:start + length]
        start %= L
        end = start + length
        if end <= L:
            return ref[start:end]
        if end <= 2 * L:
            return ref[start:] + ref[:end - L]
        return (ref * (end // L + 1))[start:end]

    def distances(self, seq, min_frac_of_seq=0.0):
        """{ref index: edit distance} for refs `seq` aligns to within max_edits."""
        n = len(seq)
        if n == 0:
            return {}

        k, e = self.k, self.max_edits
        s = seq.upper()
        n_slots = n // k
        if n_slots - e < 1:
            return self._distances_full_alignment(s, min_frac_of_seq)

        needed = n_slots - e

        best = {}
        for strand in (s, revcomp(s)):
            kmers = [strand[j * k:(j + 1) * k] for j in range(n_slots)]
            ref_sets = [self.kmer_refs.get(km) for km in kmers]
            counts = Counter(chain.from_iterable(rs for rs in ref_sets if rs))
            passing = [i for i, c in counts.items()
                       if c >= needed and self._eligible(i, n, min_frac_of_seq)]
            if not passing:
                continue

            diagonals = defaultdict(set)
            for j in range(e + 1):
                for i, p in self.dense_index.get(kmers[j], ()):
                    diag = p - j * k
                    if self.circular:
                        diag %= len(self.ref_seqs[i])
                    diagonals[i].add(diag)

            for i in passing:
                d = self._best_distance(i, strand, diagonals[i])
                if d is not None and (i not in best or d < best[i]):
                    best[i] = d

        return best

    def _best_distance(self, i, read, diagonals):
        """Lowest edit distance (<= max_edits) of `read` to ref i near the given diagonals, or None."""
        ref = self.ref_seqs[i]
        L, n, e = len(ref), len(read), self.max_edits

        for diag in diagonals:
            if self.circular or 0 <= diag:
                if self._matches_at(ref, read, diag):
                    return 0

        starts = set()
        for diag in diagonals:
            for st in range(diag - e, diag + e + 1):
                if self.circular:
                    starts.add(st % L)
                elif 0 <= st < L:
                    starts.add(st)

        best, limit = None, e
        for st in sorted(starts):
            window = self._window(i, st, n + limit)
            d = edlib.align(read, window, mode="SHW", task="distance", k=limit)["editDistance"]
            if d != -1 and (best is None or d < best):
                best = d
                if best <= 1:
                    # 0 was already ruled out above
                    break
                limit = best - 1
        return best

    def _distances_full_alignment(self, s, min_frac_of_seq):
        n, e = len(s), self.max_edits
        s_rc = revcomp(s)
        best = {}
        for i, ref in enumerate(self.ref_seqs):
            if not self._eligible(i, n, min_frac_of_seq):
                continue
            L = len(ref)
            if self.circular:
                target = self._window(i, 0, L + n + e - 1)
            else:
                target = ref
            dists = [edlib.align(r, target, mode="HW", task="distance", k=e)["editDistance"]
                     for r in (s, s_rc)]
            dists = [d for d in dists if d != -1]
            if dists:
                best[i] = min(dists)
        return best


def collapse_to_refs(seq_dists, seq_to_ref):
    """
    {sequence index: edit distance} -> {reference index: (lowest edit distance among its
    sequences, [indices of its sequences at that distance])}
    """
    out = {}
    for si, d in seq_dists.items():
        ri = seq_to_ref[si]
        current = out.get(ri)
        if current is None or d < current[0]:
            out[ri] = (d, [si])
        elif d == current[0]:
            current[1].append(si)
    return out


def summarize_distances(dists):
    """
    From {ref index: edit distance}, returns (best-matching ref indices, best distance,
    next-best distance among the other refs or None), or None if there are no matches.
    """
    if not dists:
        return None
    best = min(dists.values())
    hits = sorted(i for i, d in dists.items() if d == best)
    others = [d for d in dists.values() if d != best]
    next_best = min(others) if others else None
    return hits, best, next_best


### reading input ###

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


def _mate_base_name(name):
    return name[:-2] if name.endswith(("/1", "/2")) else name


def iter_pairs(r1_reader, r2_reader):
    r2_iter = iter(r2_reader)
    for r1 in r1_reader:
        r2 = next(r2_iter, None)
        if r2 is None:
            print(f"\n    Read-2 file ran out before read-1 file (at read '{r1[0]}').")
            notify_premature_exit()
        if _mate_base_name(r1[0]) != _mate_base_name(r2[0]):
            print(f"\n    Read names don't match between read files: '{r1[0]}' vs '{r2[0]}'.")
            print("    The read files need to be in the same order.")
            notify_premature_exit()
        yield r1, r2
    if next(r2_iter, None) is not None:
        print("\n    Read-1 file ran out before read-2 file.")
        notify_premature_exit()


def batched_with_progress(records, reader, batch_size):
    """Yields (batch, bytes of the reader's file consumed so far)."""
    it = iter(records)
    while True:
        batch = list(islice(it, batch_size))
        if not batch:
            return
        yield batch, reader.bytes_read()


### workers ###

# per-process state, set once by init_worker()
_W = {}


def init_worker(matcher, seq_to_ref, min_frac_of_seq, paired, keep_seqs):
    _W.update(matcher=matcher, seq_to_ref=seq_to_ref, min_frac_of_seq=min_frac_of_seq,
              paired=paired, keep_seqs=keep_seqs)


def process_batch(batch_and_pos):
    """
    Returns (num records in batch, bytes position, results), where results holds
    (read_id, lengths, hit ref indices, hit seq indices, best edit distance, next-best
    edit distance or None, records-or-None) for each read/pair assigned to >= 1 ref.

    Hit seqs are the sequences of the hit refs at the best distance (for pairs, from
    either mate). Records only travel back for unique hits when reads are being written.
    """
    batch, pos = batch_and_pos
    W = _W
    matcher, seq_to_ref = W["matcher"], W["seq_to_ref"]
    min_frac = W["min_frac_of_seq"]

    results = []
    if W["paired"]:
        for r1, r2 in batch:

            c1 = collapse_to_refs(matcher.distances(r1[1], min_frac), seq_to_ref)
            if not c1:
                continue
            c2 = collapse_to_refs(matcher.distances(r2[1], min_frac), seq_to_ref)
            if not c2:
                continue
            # pairs are ranked by the mates' combined edit distance
            shared = c1.keys() & c2.keys()
            summary = summarize_distances({ri: c1[ri][0] + c2[ri][0] for ri in shared})
            if summary is None:
                continue
            hits, best, next_best = summary
            hit_seqs = [si for ri in hits for si in sorted(set(c1[ri][1]) | set(c2[ri][1]))]
            recs = (r1, r2) if (W["keep_seqs"] and len(hits) == 1) else None
            results.append((_mate_base_name(r1[0]), f"{len(r1[1])},{len(r2[1])}",
                            hits, hit_seqs, best, next_best, recs))
    else:
        for rec in batch:
            collapsed = collapse_to_refs(matcher.distances(rec[1], min_frac), seq_to_ref)
            summary = summarize_distances({ri: v[0] for ri, v in collapsed.items()})
            if summary is None:
                continue
            hits, best, next_best = summary
            hit_seqs = [si for ri in hits for si in sorted(collapsed[ri][1])]
            recs = (rec,) if (W["keep_seqs"] and len(hits) == 1) else None
            results.append((rec[0], str(len(rec[1])), hits, hit_seqs, best, next_best, recs))

    return len(batch), pos, results


### outputs ###

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


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed):
    """
    Stops if `output_dir` exists (unless forced, in which case it's replaced), then creates
    it and logs the command run there. Replacing the whole directory means outputs from an
    earlier run (e.g., reads written with --write-reads) can't be left beside new ones.
    """
    check_if_output_dir_exists(output_dir, force_overwrite)
    os.makedirs(output_dir)
    log_command_run(full_cmd_executed, output_dir, output_paths(output_dir, output_prefix)["log"])


class SortedRowWriter:
    """
    Collects (sort key, line) rows and writes them to `path` after `header`, sorted by key
    (a tuple of ints). Up to `buffer_rows` rows are held in memory; beyond that, sorted
    chunks are spilled to a temporary directory beside `path` and merged at the end, so
    memory stays bounded however many rows there are. Both the sort and the merge are
    stable, so rows with equal keys keep the order they were added in.
    """

    def __init__(self, path, header, buffer_rows=None):
        self.path = path
        self.header = header
        self.buffer_rows = SORT_BUFFER_ROWS if buffer_rows is None else buffer_rows
        self.rows = []
        self.chunk_paths = []
        self.key_len = None
        self._tmp_dir = None

    def add(self, key, line):
        if self.key_len is None:
            self.key_len = len(key)
        self.rows.append((key, line))
        if len(self.rows) >= self.buffer_rows:
            self._spill()

    def _spill(self):
        if self._tmp_dir is None:
            self._tmp_dir = tempfile.mkdtemp(prefix=".sorting-", dir=os.path.dirname(self.path) or ".")
        self.rows.sort(key=itemgetter(0))
        chunk_path = os.path.join(self._tmp_dir, f"chunk-{len(self.chunk_paths)}.tsv")
        with open(chunk_path, "w") as f:
            for key, line in self.rows:
                f.write("\t".join(map(str, key)) + "\t" + line + "\n")
        self.chunk_paths.append(chunk_path)
        self.rows = []

    def _read_chunk(self, chunk_path):
        with open(chunk_path) as f:
            for raw in f:
                parts = raw.rstrip("\n").split("\t", self.key_len)
                yield tuple(int(x) for x in parts[:self.key_len]), parts[self.key_len]

    def finish(self):
        try:
            with open(self.path, "w") as out:
                out.write("\t".join(self.header) + "\n")
                if not self.chunk_paths:
                    self.rows.sort(key=itemgetter(0))
                    rows = self.rows
                else:
                    if self.rows:
                        self._spill()
                    # chunks are passed in the order they were spilled, so ties stay in input order
                    rows = heapq.merge(*(self._read_chunk(c) for c in self.chunk_paths),
                                       key=itemgetter(0))
                for _, line in rows:
                    out.write(line + "\n")
        finally:
            self.cleanup()

    def cleanup(self):
        self.rows = []
        if self._tmp_dir is not None:
            shutil.rmtree(self._tmp_dir, ignore_errors=True)
            self._tmp_dir = None


class ReadWriter:
    """Writes uniquely assigned reads to per-reference files, opened lazily."""

    def __init__(self, reads_dir, refs, paired):
        self.reads_dir = reads_dir
        self.refs = refs
        self.paired = paired
        self.handles = {}
        os.makedirs(reads_dir, exist_ok=True)

    def write(self, ref_idx, recs):
        if ref_idx not in self.handles:
            safe = self.refs[ref_idx][0].replace("/", "_")
            ext = "fastq" if recs[0][2] is not None else "fasta"
            if self.paired:
                names = [f"{safe}_R1.{ext}", f"{safe}_R2.{ext}"]
            else:
                names = [f"{safe}.{ext}"]
            self.handles[ref_idx] = [open(os.path.join(self.reads_dir, n), "w") for n in names]

        for handle, (name, seq, qual) in zip(self.handles[ref_idx], recs):
            if qual is not None:
                handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
            else:
                handle.write(f">{name}\n{seq}\n")

    def close(self):
        for hs in self.handles.values():
            for h in hs:
                h.close()


### driver ###

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

    print()
    building = spinner("Building indexes...", "Built indexes ") if show_progress else nullcontext()
    with building:
        if max_edits > 0:
            matcher = EditMatcher(seqs, k, circular, max_edits)
        else:
            matcher = ExactMatcher(seqs, k, NUM_KMER_SAMPLES, circular)

    worker_args = (matcher, seq_to_ref, min_frac_of_seq, paired, write_reads)

    # setting up reads input
    r1_reader = FastxReader(read_1)
    r2_reader = FastxReader(read_2) if paired else None
    records = iter_pairs(r1_reader, r2_reader) if paired else iter(r1_reader)
    batches = batched_with_progress(records, r1_reader, BATCH_SIZE)

    pool = None
    if jobs > 1:
        # using 'spawn' rather than 'fork' to match gen-reads and summarize-assembly (fork can
        # deadlock from multi-threaded contexts and warns on Python 3.12+). the indexes are pickled
        # once per worker at startup; identical k-mer sets are shared objects, so pickle keeps
        # that compact even when the sequences are highly similar
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(jobs, initializer=init_worker, initargs=worker_args)
        # imap keeps results in input order
        batch_results = pool.imap(process_batch, batches)
    else:
        init_worker(*worker_args)
        batch_results = map(process_batch, batches)

    unit = "pairs" if paired else "reads"
    n_total = n_unique = n_ambig = 0
    unique_counts, ambig_counts = Counter(), Counter()
    seq_unique_counts, seq_ambig_counts = Counter(), Counter()
    writer = ReadWriter(paths["reads_dir"], refs, paired) if write_reads else None

    print(f"\n    Assigning {unit} with {plural(jobs, 'job')}...\n")
    pbar = tqdm(total=r1_reader.total_bytes, unit="B", unit_scale=True, unit_divisor=1024,
                ncols=80, disable=not show_progress, file=sys.stdout,
                bar_format="    {l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]")
    pbar.set_postfix_str(f"0 {unit}")

    header = ["read_id", "read_lengths" if paired else "read_length", "status", "num_hits",
              "matching_refs", "edit_distance", "next_best_edit_distance"]
    if not per_seq:
        header.append("matching_seqs")
    # rows are sorted by read length (combined for pairs), longest first, then by edit distance
    hits_writer = SortedRowWriter(paths["hits"], header)

    try:
        last_pos = 0
        for n_in_batch, pos, results in batch_results:
            n_total += n_in_batch
            pbar.update(pos - last_pos)
            last_pos = pos
            pbar.set_postfix_str(f"{n_total:,} {unit}", refresh=False)

            for read_id, lengths, hits, hit_seqs, best, next_best, recs in results:
                if len(hits) == 1:
                    status = "unique"
                    n_unique += 1
                    unique_counts[hits[0]] += 1
                    seq_unique_counts.update(hit_seqs)
                    if writer is not None:
                        writer.write(hits[0], recs)
                else:
                    status = "ambiguous"
                    n_ambig += 1
                    ambig_counts.update(hits)
                    seq_ambig_counts.update(hit_seqs)

                line = (f"{read_id}\t{lengths}\t{status}\t{len(hits)}\t"
                        f"{','.join(refs[i][0] for i in hits)}\t{best}\t"
                        f"{'NA' if next_best is None else next_best}")
                if not per_seq:
                    line += "\t" + ",".join(seq_label(si) for si in hit_seqs)
                total_length = sum(int(x) for x in lengths.split(","))
                hits_writer.add((-total_length, best), line)

        # finishing the bar cleanly (the buffered position can lag the file size)
        pbar.update(r1_reader.total_bytes - last_pos)

    except BaseException:
        if pool is not None:
            pool.terminate()
            pool = None
        hits_writer.cleanup()
        raise

    finally:
        pbar.close()
        if pool is not None:
            pool.close()
            pool.join()
        if writer is not None:
            writer.close()
        r1_reader.close()
        if r2_reader is not None:
            r2_reader.close()

    hits_writer.finish()

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
