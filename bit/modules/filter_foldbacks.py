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

import os
import sys
import threading

from bit.modules.general import (check_if_output_dir_exists, color_text, log_command_run,
                                 notify_premature_exit)
from bit.modules.match_reads import FastxReader

# short reads: 13-mers, 16 shared (at least 28 bp shared as one stretch). Tested on pairs cut
# from real PTA molecules and from reference genomes (Oct 2026): against 17-mers with 20 shared
# (36 bp), it catches more fold-backs (2x150: 27% of pairs with a fold in a read vs 22%; 2x250:
# 53% vs 48%) for more false calls at real short inverted repeats in some genomes (Pseudomonas,
# rich in short palindromic repeats: 0.28% of 2x150 and 0.88% of 2x250 pairs vs 0.07% / 0.26%;
# most others under 0.1%). Shorter still (24 bp) gained little more for twice the false calls
KMER_SIZE = 13
MIN_SHARED = 16
MIN_LENGTH = 50

# long reads (hifi=True): shifia's --pta settings (shifia/modules/pta.py): 21-mers, about 1 in 4
# of them used, every fold cut however narrow its unpaired middle, pieces under 100 bp dropped
LONG_K = 21
LONG_DENSITY = 4
LONG_MIN_LENGTH = 100


def full_prefix(output_prefix):
    """ the prefix as it goes on the file names: a '-' is added unless it already ends in '-' or '_' """
    if output_prefix and not output_prefix.endswith(("-", "_")):
        return output_prefix + "-"
    return output_prefix


def output_paths(output_dir, output_prefix="", paired=True, fastq=True):
    """Output locations within `output_dir`, with `output_prefix` (see full_prefix) prepended to each name."""
    ext = "fastq.gz" if fastq else "fasta.gz"
    p = os.path.join(output_dir, full_prefix(output_prefix))
    out = {"summary": f"{p}summary.tsv", "log": f"{p}command-execution-info.txt"}
    if paired:
        out.update(r1=f"{p}filtered_R1.{ext}", r2=f"{p}filtered_R2.{ext}",
                   singles=f"{p}foldback-mates.{ext}")
    else:
        out.update(reads=f"{p}filtered.{ext}")
    return out


def _all_output_paths(output_dir, output_prefix):
    """ every file a run with this prefix can write, whatever its input """
    paths = set()
    for paired in (True, False):
        for fastq in (True, False):
            paths.update(output_paths(output_dir, output_prefix, paired, fastq).values())
    return sorted(paths)


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed):
    """
    Without a prefix, the output directory must be new (or, with force_overwrite, it is replaced),
    as with bit's other programs. With a prefix, the directory can already exist (e.g., holding
    other samples' output); only this prefix's own files must not, or with force_overwrite they
    are replaced (nothing else in the directory is touched).
    """
    if output_prefix and os.path.isdir(output_dir):
        existing = [p for p in _all_output_paths(output_dir, output_prefix) if os.path.exists(p)]
        if existing and not force_overwrite:
            print(f"\n    {color_text('Output files with this prefix already exist:', 'yellow')}")
            for p in existing:
                print(f"      {p}")
            print("\n    Please specify a different output prefix or directory, or add the `-F/--force-overwrite` flag"
                  "\n    (which replaces only these files).")
            notify_premature_exit()
        for p in existing:
            os.remove(p)
    else:
        check_if_output_dir_exists(output_dir, force_overwrite)
        os.makedirs(output_dir)
    log_command_run(full_cmd_executed, output_dir, output_paths(output_dir, output_prefix)["log"])


class Counts:
    """ what was found, for the summary """

    FIELDS = ["reads_or_pairs_in", "foldback_pairs", "pairs_with_fold_in_a_read",
              "pairs_with_fold_between_mates", "reads_cut_at_a_fold", "folds_cut",
              "mates_written_as_singles", "reads_dropped_short", "pairs_out", "reads_out",
              "bases_in", "bases_out"]

    def __init__(self, values=None, middles=None):
        for f, v in zip(self.FIELDS, values or [0] * len(self.FIELDS)):
            setattr(self, f, v)
        self.middles = list(middles or [])    # unpaired middle of each cut read's first fold (long reads)

    def middle_quartiles(self):
        """ (25th, 50th, 75th percentile) of the unpaired middles, or None """
        if not self.middles:
            return None
        m = sorted(self.middles)
        return tuple(m[min(len(m) - 1, int(q * len(m)))] for q in (0.25, 0.5, 0.75))

    def rows(self):
        return [(f, getattr(self, f)) for f in self.FIELDS]


def _run(fn, total, unit, show_progress):
    """
    Runs one pass of the compiled core with bit's usual capped, kept bar over the bytes of the
    input and a running count. The pass runs in a thread of its own (it lets go of the GIL)
    while this one draws the bar.
    """
    from bit._native import Progress
    prog = Progress()
    if not show_progress:
        return fn(prog)
    result, error = [None], [None]

    def work():
        try:
            result[0] = fn(prog)
        except BaseException as e:  # passed on below
            error[0] = e

    from tqdm import tqdm # type: ignore
    print()
    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    with tqdm(total=total, unit="B", unit_scale=True, unit_divisor=1024, ncols=80, file=sys.stdout,
              bar_format="    {l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]") as bar:
        bar.set_postfix_str(f"0 {unit}")
        while thread.is_alive():
            thread.join(0.2)
            bar.update(min(prog.bytes, total) - bar.n)
            bar.set_postfix_str(f"{prog.reads:,} {unit}", refresh=False)
        if error[0] is None:
            bar.set_postfix_str(f"{prog.reads:,} {unit}", refresh=False)
            bar.update(total - bar.n)
            bar.refresh()
    if error[0] is not None:
        raise error[0]
    return result[0]


def filter_foldbacks(read_1, read_2=None, output_dir="filter-foldbacks", output_prefix="",
                     min_length=None, k=KMER_SIZE, min_shared=MIN_SHARED, show_progress=True,
                     hifi=False):
    """
    Writes the filtered reads to output_dir and returns (Counts, output paths). With hifi, the
    reads are long single reads (read_2 must be None) and k and min_shared aren't used (see
    the module's notes); min_length defaults to 50 for short reads and 100 for long ones.
    Reading, fold finding and writing are in the compiled core (filter_foldback_files), spread
    over the threads set with bit._native.set_threads.
    """
    from bit._native import filter_foldback_files

    if hifi and read_2 is not None:
        raise ValueError("long reads (hifi) are single reads: no read 2")
    if min_length is None:
        min_length = LONG_MIN_LENGTH if hifi else MIN_LENGTH
    paired = read_2 is not None
    peek = FastxReader(read_1)
    first = next(iter(peek), None)
    peek.close()
    fastq = first is not None and first[2] is not None
    paths = output_paths(output_dir, output_prefix, paired, fastq)
    outputs = [paths["r1"], paths["r2"], paths["singles"]] if paired else [paths["reads"]]

    def work(prog):
        return filter_foldback_files(read_1, read_2, outputs, LONG_K if hifi else k, min_shared,
                                     min_length, hifi=hifi, density=LONG_DENSITY, progress=prog)

    try:
        total = sum(os.path.getsize(p) for p in ([read_1, read_2] if paired else [read_1]))
        values, written, middles = _run(work, total, "pairs" if paired else "reads", show_progress)
    except ValueError as e:  # a problem with the input files
        print(f"\n    {e}")
        notify_premature_exit()
    counts = Counts(values, middles)
    if paired and not written[2]:
        paths["singles"] = None
    _write_summary(paths, counts, paired, long=hifi)
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
