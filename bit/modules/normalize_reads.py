"""
Core logic for `bit normalize-reads`.

Caps read depth for assembly while keeping reads (or pairs) that hold sequence still under
the target depth, including from low-abundance members and minority variants. Builds on
shifia's normalizer (github.com/AstrobioMike/shifia, shifia/modules/normalize.py); the passes
over the reads are in the compiled core (rust/src/normrun.rs).

Each read, or each pair (both mates together), is reduced to sampled canonical k-mers
("anchors": by default 25-mers, about 1 in 2 of them). Median-depth normalization (bbnorm,
khmer) keeps or drops a read by its median k-mer depth, so a read from a rare genome or a
minority variant that is mostly shared, deep sequence has a high median and is thinned like
the rest. Here a read is judged by its least-covered part instead.

Two passes:
  1. counting: every read's distinct anchors go into a count-min sketch (fixed memory,
     `--memory`), so each anchor's read count in the whole input is known before any read is
     judged. (shifia judges its long reads in one pass: a HiFi read is long and accurate enough
     to tell new sequence from its errors on its own, and kept reads are few at the depths it
     works at. A 150 bp read holding one SNP looks like an error until the count is known, and
     in a metagenome most reads are kept.)
  2. judging and writing: an anchor is solid when it's in at least 2 reads and at least 1% of
     the upper quartile of the reads of the anchors within 2k bases of it (an error's
     neighbors are the same sequence without it); others are errors and count for nothing. A
     read's depth is that of its solid anchors at the `novel_fraction` quantile (its 10% least
     covered), and it's kept with probability target / that depth (always when at or under the
     target). With 3 or more solid anchors on a minor branch (in 3+ reads, at 5-50% of the
     read's typical depth: a minority variant), it's kept with at least probability
     (target / 2) / the branch's depth. A read with anchors but none solid (at most a read or
     two's worth of anything) is kept. The coin for each read comes from a hash of its own
     sequence, so what's kept doesn't depend on read order or threads, and needs nothing held
     in memory but the sketch.
"""

import os
import sys
import threading

from tqdm import tqdm # type: ignore

from bit.modules.general import (check_if_output_dir_exists,
                                 log_command_run,
                                 notify_premature_exit)

TARGET_DEPTH = 100
KMER_SIZE = 25
DENSITY = 2
NOVEL_FRACTION = 0.1
MIN_READS = 2
SEED = 1
ERROR_WEIGHT = 1.0
ERROR_SHARE = 0.01
# the counting pass's sketch, when --memory isn't given: about as many bytes as the input on
# disk, within these bounds and at most half the computer's memory
MIN_SKETCH = 256 * 2**20
MAX_SKETCH = 8 * 2**30
# a sketch fuller than this counts rare anchors high often enough to matter
FULL_SKETCH = 0.5


def output_paths(output_dir, output_prefix="", paired=True, fastq=True):
    """Output locations within `output_dir`, with `output_prefix` prepended to each name."""
    ext = "fastq.gz" if fastq else "fasta.gz"
    p = os.path.join(output_dir, output_prefix)
    out = {"summary": f"{p}summary.tsv", "log": f"{p}command-execution-info.txt",
           "depth_table": f"{p}depths.tsv.gz"}
    if paired:
        out.update(r1=f"{p}normalized_R1.{ext}", r2=f"{p}normalized_R2.{ext}")
    else:
        out.update(reads=f"{p}normalized.{ext}")
    return out


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed):
    check_if_output_dir_exists(output_dir, force_overwrite)
    os.makedirs(output_dir)
    log_command_run(full_cmd_executed, output_dir, output_paths(output_dir, output_prefix)["log"])


def _is_fastq(path):
    import gzip
    opener = gzip.open if _is_gzipped(path) else open
    with opener(path, "rt") as f:
        for line in f:
            if line.strip():
                return line.startswith("@")
    return True


def _is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def default_sketch_bytes(read_paths):
    """ about as many bytes as the input on disk, within MIN_SKETCH-MAX_SKETCH and at most
    half the computer's memory """
    on_disk = sum(os.path.getsize(p) for p in read_paths)
    size = min(max(on_disk, MIN_SKETCH), MAX_SKETCH)
    try:
        ram = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        size = min(size, ram // 2)
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        pass
    return int(size)


def _run(desc, total, show_progress, fn, *args, **kwargs):
    """
    Runs one pass of the compiled core in a thread of its own (it lets go of the GIL) while
    this one draws bit's usual capped, kept progress bar over the input's bytes.
    """
    from bit import _native as native
    if not show_progress:
        return fn(*args, **kwargs)
    prog = native.Progress()
    result, error = [None], [None]

    def work():
        try:
            result[0] = fn(*args, progress=prog, **kwargs)
        except BaseException as e:  # passed on below
            error[0] = e

    print()
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    with tqdm(total=total, desc=desc, unit="B", unit_scale=True, unit_divisor=1024, ncols=80,
              file=sys.stdout,
              bar_format="    {l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]") as bar:
        while thread.is_alive():
            thread.join(0.2)
            bar.update(min(prog.bytes, total) - bar.n)
            bar.set_postfix_str(f"{prog.reads:,} reads", refresh=False)
        if error[0] is None:
            bar.set_postfix_str(f"{prog.reads:,} reads", refresh=False)
            bar.update(total - bar.n)
            bar.refresh()
    if error[0] is not None:
        raise error[0]
    return result[0]


class Counts:
    """ what was done, for the summary """

    FIELDS = ["reads_or_pairs_in", "bases_in", "reads_or_pairs_out", "bases_out",
              "no_anchors", "kept_without_solid_kmers", "kmer_occurrences",
              "kmer_occurrences_judged_errors", "median_depth_in", "median_depth_out",
              "sketch_bytes", "sketch_fill", "clean_weight"]

    def __init__(self, totals):
        for f in self.FIELDS:
            setattr(self, f, totals.get(f))

    @property
    def sketch_too_full(self):
        return self.sketch_fill is not None and self.sketch_fill > FULL_SKETCH

    @property
    def error_share(self):
        if not self.kmer_occurrences:
            return 0.0
        return self.kmer_occurrences_judged_errors / self.kmer_occurrences

    def rows(self):
        rows = []
        for f in self.FIELDS:
            v = getattr(self, f)
            if v is None:
                continue
            if f == "sketch_fill":
                v = round(v, 4)
            rows.append((f, v))
            if f == "kmer_occurrences_judged_errors":
                rows.append(("share_of_kmers_judged_errors", round(self.error_share, 4)))
        return rows


def normalize(read_1, read_2=None, output_dir="normalize-reads", output_prefix="",
              target_depth=TARGET_DEPTH, k=KMER_SIZE, density=DENSITY,
              novel_fraction=NOVEL_FRACTION, minor_branches=True, min_reads=MIN_READS,
              sketch_bytes=None, depth_table=False, seed=SEED, error_weight=ERROR_WEIGHT,
              error_share=ERROR_SHARE, show_progress=True):
    """
    Writes the normalized reads to output_dir and returns (Counts, output paths). With
    read_2, mates are judged and kept together.
    """
    from bit import _native as native

    paired = read_2 is not None
    paths = output_paths(output_dir, output_prefix, paired, _is_fastq(read_1))
    inputs = [read_1] + ([read_2] if paired else [])
    total = sum(os.path.getsize(p) for p in inputs)
    out_1 = paths["r1"] if paired else paths["reads"]
    out_2 = paths["r2"] if paired else None
    table = paths["depth_table"] if depth_table else None
    if not depth_table:
        del paths["depth_table"]
    minor_anchors = 3 if minor_branches else 0

    if sketch_bytes is None:
        sketch_bytes = default_sketch_bytes(inputs)
    try:
        # (the compiled core adds a sketch half this size of anchors by the strand reads hold
        # them on, for telling systematic, one-strand errors from sequence)
        counted = _run("Counting", total, show_progress, native.count_anchors,
                       read_1, read_2, k, density, int(sketch_bytes))
        totals = _run("Writing ", total, show_progress, counted.normalize,
                      read_1, read_2, out_1, out_2, target_depth,
                      novel_fraction=novel_fraction, minor_anchors=minor_anchors,
                      min_reads=min_reads, seed=seed, error_weight=error_weight,
                      error_share=error_share,
                      depth_table=table)
    except ValueError as e:
        if "Read names" in str(e) or "ran out before" in str(e):
            print(f"\n    {e}", file=sys.stderr)
            notify_premature_exit()
        raise

    counts = Counts(totals)
    _write_summary(paths, counts)
    return counts, paths


def _write_summary(paths, counts):
    with open(paths["summary"], "w") as out:
        out.write("measure\tvalue\n")
        for name, value in counts.rows():
            out.write(f"{name}\t{value}\n")
