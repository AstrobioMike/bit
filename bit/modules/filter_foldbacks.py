"""
Core logic for `bit filter-foldbacks`.

Whole-genome amplification by PTA (and MDA) leaves many molecules folded back: a stretch of the
template, then its own reverse complement, where the copy switched to the opposite strand a few
hundred bases away. Both arms are real sequence, but the join between them is not, and neither is
the pairing of reads from such a fragment (its two ends are the same strand of one place, so the
pair looks like an inverted rearrangement with a tiny insert).

Found from the reads alone (no reference; the k-mer work is in the compiled core, see
rust/src/foldback.rs):
  - a read with a fold in it: k-mers whose reverse complement is also in the read, agreeing on
    where the read turns. The read is cut to its longer arm, ending where that arm's paired
    stretch ends, so the join and anything unpaired beside it (often untemplated bases) go.
  - a pair whose fold is between the mates: read 2 holds stretches of read 1 in the same
    orientation. Normal mates are on opposite strands, so the only sequence they can share (when
    the fragment is shorter than the two reads) is read 1 against read 2 reverse complemented.

Pairs with neither are written as pairs. Mates of fold-back pairs are written as single reads
(cut where their own fold was, and kept if still at least `min_length` long), since each is
still real sequence even though the pairing isn't.

Long reads (`hifi=True`, e.g., PacBio HiFi): the fold finding and cutting of `shifia --pta`
(github.com/AstrobioMike/shifia), on the same sampled 21-mers, so the two cut reads at the same
places. A long read can hold both arms and the unpaired middle between them (~50-400 bp in PTA
reads), so every fold is found within the read; it is cut to its longer arm, ending where that
arm's paired stretch ends, and checked again (a read can fold more than once). As a long read is
mostly not fold, a fold must also pair a quarter of the anchors on the read's shorter side of it,
which leaves genomic inverted repeats inside long flanks alone. Keep this in step with shifia's
shifia/modules/pta.py and rust/src/pta.rs.
"""

import gzip
import os
from itertools import islice

from tqdm import tqdm # type: ignore

from bit.modules.general import check_if_output_dir_exists, log_command_run
from bit.modules.match_reads import FastxReader, iter_pairs

# short reads: 13-mers, 16 shared (at least 28 bp shared as one stretch). Tested on pairs cut
# from real PTA molecules and from reference genomes (Oct 2026): against 17-mers with 20 shared
# (36 bp), it catches more fold-backs (2x150: 27% of pairs with a fold in a read vs 22%; 2x250:
# 53% vs 48%) for more false calls at real short inverted repeats in some genomes (Pseudomonas,
# rich in short palindromic repeats: 0.28% of 2x150 and 0.88% of 2x250 pairs vs 0.07% / 0.26%;
# most others under 0.1%). Shorter still (24 bp) gained little more for twice the false calls
KMER_SIZE = 13
MIN_SHARED = 16
BATCH_SIZE = 20_000
MIN_LENGTH = 50

# long reads (hifi=True): shifia's --pta settings (shifia/modules/pta.py): 21-mers, about 1 in 4
# of them used, every fold cut however narrow its unpaired middle, pieces under 100 bp dropped
LONG_K = 21
LONG_DENSITY = 4
LONG_MIN_LENGTH = 100
LONG_BATCH_SIZE = 2_000


def output_paths(output_dir, output_prefix="", paired=True, fastq=True):
    """Output locations within `output_dir`, with `output_prefix` prepended to each name."""
    ext = "fastq.gz" if fastq else "fasta.gz"
    p = os.path.join(output_dir, output_prefix)
    out = {"summary": f"{p}summary.tsv", "log": f"{p}command-execution-info.txt"}
    if paired:
        out.update(r1=f"{p}filtered_R1.{ext}", r2=f"{p}filtered_R2.{ext}",
                   singles=f"{p}foldback-mates.{ext}")
    else:
        out.update(reads=f"{p}filtered.{ext}")
    return out


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed):
    check_if_output_dir_exists(output_dir, force_overwrite)
    os.makedirs(output_dir)
    log_command_run(full_cmd_executed, output_dir, output_paths(output_dir, output_prefix)["log"])


def _write(handle, name, seq, qual):
    if qual is not None:
        handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
    else:
        handle.write(f">{name}\n{seq}\n")


def _cut(rec, start, end):
    name, seq, qual = rec
    return name, seq[start:end], (qual[start:end] if qual is not None else None)


class Counts:
    """ what was found, for the summary """

    FIELDS = ["reads_or_pairs_in", "foldback_pairs", "pairs_with_fold_in_a_read",
              "pairs_with_fold_between_mates", "reads_cut_at_a_fold", "folds_cut",
              "mates_written_as_singles", "reads_dropped_short", "pairs_out", "reads_out",
              "bases_in", "bases_out"]

    def __init__(self):
        for f in self.FIELDS:
            setattr(self, f, 0)
        self.middles = []    # unpaired middle of each cut read's first fold (long reads)

    def middle_quartiles(self):
        """ (25th, 50th, 75th percentile) of the unpaired middles, or None """
        if not self.middles:
            return None
        m = sorted(self.middles)
        return tuple(m[min(len(m) - 1, int(q * len(m)))] for q in (0.25, 0.5, 0.75))

    def rows(self):
        return [(f, getattr(self, f)) for f in self.FIELDS]


def classify_pair(scan, min_shared=MIN_SHARED):
    """
    (fold in read 1 or 2, fold between the mates) from one scan_foldbacks row. Sequence shared in
    the same orientation calls a fold whatever else the mates share: in a fold-back fragment,
    a read crossing the fold also shares the other arm with its mate the normal way.
    """
    _, _, f1, _, _, f2, same, _ = scan
    return f1 > 0 or f2 > 0, same >= min_shared


def filter_foldbacks(read_1, read_2=None, output_dir="filter-foldbacks", output_prefix="",
                     min_length=None, k=KMER_SIZE, min_shared=MIN_SHARED, batch_size=None,
                     show_progress=True, hifi=False):
    """
    Writes the filtered reads to output_dir and returns (Counts, output paths). With hifi, the
    reads are long single reads (read_2 must be None) and k and min_shared aren't used (see
    the module's notes); min_length defaults to 50 for short reads and 100 for long ones.
    """
    if hifi:
        if read_2 is not None:
            raise ValueError("long reads (hifi) are single reads: no read 2")
        return _filter_long(read_1, output_dir, output_prefix,
                            LONG_MIN_LENGTH if min_length is None else min_length,
                            batch_size or LONG_BATCH_SIZE, show_progress)
    from bit._native import scan_foldbacks

    if min_length is None:
        min_length = MIN_LENGTH
    if batch_size is None:
        batch_size = BATCH_SIZE
    paired = read_2 is not None
    r1_reader = FastxReader(read_1)
    r2_reader = FastxReader(read_2) if paired else None
    peek = FastxReader(read_1)
    first = next(iter(peek), None)
    peek.close()
    fastq = first is not None and first[2] is not None
    paths = output_paths(output_dir, output_prefix, paired, fastq)
    counts = Counts()

    def opener(path):
        return gzip.open(path, "wt", compresslevel=1)

    records = iter_pairs(r1_reader, r2_reader) if paired else ((r, None) for r in r1_reader)
    bar = tqdm(total=r1_reader.total_bytes, unit="B", unit_scale=True, disable=not show_progress,
               desc="    Scanning reads", leave=False)

    handles = {key: opener(paths[key]) for key in (("r1", "r2", "singles") if paired else ("reads",))}
    try:
        while True:
            batch = list(islice(records, batch_size))
            if not batch:
                break
            seqs1 = [r1[1] for r1, _ in batch]
            seqs2 = [r2[1] for _, r2 in batch] if paired else None
            scans = scan_foldbacks(seqs1, seqs2, k, min_shared)
            for (r1, r2), scan in zip(batch, scans):
                counts.reads_or_pairs_in += 1
                counts.bases_in += len(r1[1]) + (len(r2[1]) if paired else 0)
                s1, e1, f1, s2, e2, f2, _, _ = scan
                counts.reads_cut_at_a_fold += (f1 > 0) + (f2 > 0)
                counts.folds_cut += f1 + f2
                if not paired:
                    if f1 == 0:
                        _write(handles["reads"], *r1)
                        counts.reads_out += 1
                        counts.bases_out += len(r1[1])
                    elif e1 - s1 >= min_length:
                        rec = _cut(r1, s1, e1)
                        _write(handles["reads"], *rec)
                        counts.reads_out += 1
                        counts.bases_out += len(rec[1])
                    else:
                        counts.reads_dropped_short += 1
                    continue
                in_read, between = classify_pair(scan, min_shared)
                if not in_read and not between:
                    _write(handles["r1"], *r1)
                    _write(handles["r2"], *r2)
                    counts.pairs_out += 1
                    counts.bases_out += len(r1[1]) + len(r2[1])
                    continue
                counts.foldback_pairs += 1
                counts.pairs_with_fold_in_a_read += in_read
                counts.pairs_with_fold_between_mates += between and not in_read
                for rec, (s, e) in ((r1, (s1, e1)), (r2, (s2, e2))):
                    rec = _cut(rec, s, e)
                    if len(rec[1]) >= min_length:
                        _write(handles["singles"], *rec)
                        counts.mates_written_as_singles += 1
                        counts.bases_out += len(rec[1])
                    else:
                        counts.reads_dropped_short += 1
            bar.update(r1_reader.bytes_read() - bar.n)
    finally:
        bar.close()
        for h in handles.values():
            h.close()
        r1_reader.close()
        if r2_reader is not None:
            r2_reader.close()

    _write_summary(paths, counts, paired)
    return counts, paths


SINGLE_FIELDS = ("reads_or_pairs_in", "reads_cut_at_a_fold", "folds_cut", "reads_dropped_short",
                 "reads_out", "bases_in", "bases_out")


def _write_summary(paths, counts, paired, long=False):
    """ summary.tsv: the pair counts too with paired reads; the unpaired middles with long ones """
    keep = SINGLE_FIELDS
    with open(paths["summary"], "w") as out:
        out.write("measure\tvalue\n")
        for name, value in counts.rows():
            if paired or name in keep:
                out.write(f"{name}\t{value}\n")
        q = counts.middle_quartiles() if long else None
        if q:
            for label, v in zip(("p25", "median", "p75"), q):
                out.write(f"unpaired_middle_{label}\t{v}\n")


def _filter_long(read_1, output_dir, output_prefix, min_length, batch_size, show_progress):
    """ long reads (hifi): each read cut at its folds as shifia --pta does (see the notes above) """
    from bit._native import trim_long_folds

    reader = FastxReader(read_1)
    peek = FastxReader(read_1)
    first = next(iter(peek), None)
    peek.close()
    fastq = first is not None and first[2] is not None
    paths = output_paths(output_dir, output_prefix, paired=False, fastq=fastq)
    counts = Counts()
    bar = tqdm(total=reader.total_bytes, unit="B", unit_scale=True, disable=not show_progress,
               desc="    Scanning reads", leave=False)
    handle = gzip.open(paths["reads"], "wt", compresslevel=1)
    records = iter(reader)
    try:
        while True:
            batch = list(islice(records, batch_size))
            if not batch:
                break
            cuts = trim_long_folds([r[1] for r in batch], k=LONG_K, density=LONG_DENSITY)
            for rec, (start, end, folds, _, middle) in zip(batch, cuts):
                counts.reads_or_pairs_in += 1
                counts.bases_in += len(rec[1])
                if folds:
                    counts.reads_cut_at_a_fold += 1
                    counts.folds_cut += folds
                    if len(counts.middles) < 1_000_000:
                        counts.middles.append(middle)
                    rec = _cut(rec, start, end)
                    if len(rec[1]) < min_length:
                        counts.reads_dropped_short += 1
                        continue
                _write(handle, *rec)
                counts.reads_out += 1
                counts.bases_out += len(rec[1])
            bar.update(reader.bytes_read() - bar.n)
    finally:
        bar.close()
        handle.close()
        reader.close()
    _write_summary(paths, counts, paired=False, long=True)
    return counts, paths
