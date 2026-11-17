import os
from dataclasses import dataclass, field
import pyfastx #type: ignore
from collections import defaultdict
from colorama import Fore, init #type: ignore
from tqdm import tqdm #type: ignore
from bit.modules.general import check_files_are_found, report_message, sort_bam_file, spinner
from bit.modules.mapped_read_stats import (get_per_contig_pid_stats, PidStats,
                                           _check_opens_with_pysam, _pid_stats_from_native)


init(autoreset=True)


def get_cov_stats(args):

    preflight_checks(args)

    refs, contigs, header_to_ref_dict = parse_refs(args.reference_fastas, args.skip_per_contig)

    # coverage and percent IDs both from one pass over the bam
    contig_pid_stats = scan_bam(args.bam, contigs, header_to_ref_dict, args.include_non_primary,
                                want_pids=not args.skip_read_pids)

    (ref_mean_pids, ref_median_pids, ref_mapped_counts,
     contig_mean_pids, contig_median_pids, contig_mapped_counts) = summarize_pid_stats(contig_pid_stats, refs, contigs,
                                                                                       args.skip_per_contig)

    generate_output(refs, args.output_prefix, ref_mean_pids, ref_median_pids, ref_mapped_counts,
                    contigs, contig_mean_pids, contig_median_pids, contig_mapped_counts, args.skip_per_contig)


def preflight_checks(args):

    check_files_are_found(list(args.reference_fastas) + [args.bam])


@dataclass
class CoverageStats:
    length: int = 0
    total_coverage_count: int = 0
    total_bases_detected_at_all: int = 0
    total_bases_detected_at_10x: int = 0
    depth_hist: dict = field(default_factory=lambda: defaultdict(int))

    def compute_metrics(self):
        """ num_bases_detected, detection (%), detection_at_10x (%), mean_coverage, median_coverage """
        if self.length == 0:
            return 0, 0.0, 0.0, 0.0, 0.0
        detection = round(100 * self.total_bases_detected_at_all / self.length, 2)
        detection_at_10x = round(100 * self.total_bases_detected_at_10x / self.length, 2)
        mean_coverage = round(self.total_coverage_count / self.length, 2)
        median_coverage = self._median_from_hist()

        return self.total_bases_detected_at_all, detection, detection_at_10x, mean_coverage, median_coverage

    def _median_from_hist(self):
        if self.length == 0:
            return 0.0

        # getting two (possible) middle indices
        i1 = (self.length - 1) // 2
        i2 = self.length // 2

        # handling potential missing bases as 0s
        total_counted = sum(self.depth_hist.values())
        needed_extra_zeros = self.length - total_counted

        cumulative_count = 0
        v1 = v2 = None

        depths = sorted(self.depth_hist)
        if 0 not in self.depth_hist:
            depths = [0] + depths

        for depth in depths:
            count = self.depth_hist[depth]

            if depth == 0:
                count += needed_extra_zeros

            next_cumulative_count = cumulative_count + count

            if v1 is None and i1 < next_cumulative_count:
                v1 = depth
            if v2 is None and i2 < next_cumulative_count:
                v2 = depth
                break

            cumulative_count = next_cumulative_count

        return round((float(v1) + float(v2)) / 2.0, 2)


@dataclass
class ContigData:
    path: str
    header: str
    stats: CoverageStats = field(default_factory=CoverageStats)


@dataclass
class RefData:
    path: str
    headers: set = field(default_factory=set)
    contig_order: list = field(default_factory=list)
    stats: CoverageStats = field(default_factory=CoverageStats)


def parse_refs(reference_fastas, skip_per_contig=False):

    print(f"\n  {Fore.YELLOW}Parsing input reference(s)...")

    refs = []
    contigs_dict = {}
    header_to_ref_dict = {}

    with spinner("", ""):
        for fasta in reference_fastas:
            ref = RefData(fasta)
            fxi_path = fasta + ".fxi"
            fxi_existed = os.path.exists(fxi_path)
            fa = pyfastx.Fasta(fasta)
            for seq in fa:
                name = seq.name
                seq_len = len(seq)
                ref.headers.add(name)
                ref.contig_order.append(name)
                ref.stats.length += seq_len
                if not skip_per_contig:
                    contigs_dict[name] = ContigData(path=fasta, header=name, stats=CoverageStats(length=seq_len))
                header_to_ref_dict[name] = ref
            refs.append(ref)
            # clean up pyfastx index file if we created it
            if not fxi_existed and os.path.exists(fxi_path):
                os.remove(fxi_path)

    return refs, contigs_dict, header_to_ref_dict


def summarize_pid_stats(contig_pid_stats, refs, contigs_dict, skip_per_contig=False):
    """ mean and median percent ID and number of mapped reads, per ref and per contig, from the
    per-contig PidStats (all None if there are none) """

    if contig_pid_stats is None:
        return None, None, None, None, None, None

    contig_mean_pids = {}
    contig_median_pids = {}
    contig_mapped_counts = {}

    if not skip_per_contig:
        for contig_name in contigs_dict.keys():
            if contig_name in contig_pid_stats:
                stats = contig_pid_stats[contig_name]
                contig_mapped_counts[contig_name] = stats.count
                contig_mean_pids[contig_name] = round(stats.mean, 2) if stats.count > 0 else None
                contig_median_pids[contig_name] = round(stats.median, 2) if stats.count > 0 else None
            else:
                contig_mapped_counts[contig_name] = 0
                contig_mean_pids[contig_name] = None
                contig_median_pids[contig_name] = None

    ref_mean_pids = {}
    ref_median_pids = {}
    ref_mapped_counts = {}

    for ref in refs:
        ref_pid_stats = PidStats()
        for contig_name in ref.headers:
            if contig_name in contig_pid_stats:
                ref_pid_stats.merge(contig_pid_stats[contig_name])

        ref_mapped_counts[ref.path] = ref_pid_stats.count
        ref_mean_pids[ref.path] = round(ref_pid_stats.mean, 2) if ref_pid_stats.count > 0 else None
        ref_median_pids[ref.path] = round(ref_pid_stats.median, 2) if ref_pid_stats.count > 0 else None

    return ref_mean_pids, ref_median_pids, ref_mapped_counts, contig_mean_pids, contig_median_pids, contig_mapped_counts


def stats_lookup(contigs, header_to_ref_dict):
    """ {contig header: (its own CoverageStats or None, its ref's CoverageStats)} """

    lookup = {}
    for header, ref in header_to_ref_dict.items():
        contig_stats = contigs[header].stats if header in contigs else None
        lookup[header] = (contig_stats, ref.stats)

    return lookup


def add_contig_totals(lookup, wanted, totals):
    """ adds each contig's coverage totals (total coverage count, bases detected at all, at 10x,
    [(depth, bases), ...]; None where there are none) to its own stats and its ref's """

    for header, t in zip(wanted, totals):
        if t is None:
            continue
        total_cov, detected, detected_10x, hist = t
        for stats in lookup[header]:
            if stats is None:
                continue
            stats.total_coverage_count += total_cov
            stats.total_bases_detected_at_all += detected
            stats.total_bases_detected_at_10x += detected_10x
            for depth, span in hist:
                stats.depth_hist[depth] += span


class NotGroupedByReference(Exception):
    """ the bam's alignments to a contig aren't all together (it isn't sorted) """


# the alignments that add no coverage, by flag: unmapped, failing QC, duplicates and
# secondary ones (1796, what mosdepth leaves out too), or with --include-non-primary all but
# the secondary ones (1540). Supplementary alignments count either way.
COVERAGE_EXCLUDE_FLAGS = 1796
COVERAGE_EXCLUDE_FLAGS_WITH_NON_PRIMARY = 1540


def _native_bam_pass(bam_path, wanted, exclude_flags, include_non_primary, want_pids):
    """
    Coverage totals per wanted contig and (if wanted) per-contig PidStats from one pass over
    the bam by the compiled core (bit._native). Returns None if anything in the bam calls
    for the Python (pysam) path.
    """
    try:
        from bit import _native as native
    except ImportError:
        return None

    import threading
    prog = native.Progress()
    result, error = [None], [None]

    def work():
        try:
            result[0] = native.bam_coverage_and_pid_stats(str(bam_path), wanted, exclude_flags, want_pids,
                                                          include_non_primary, True, progress=prog)
        except BaseException as e:  # passed on below
            error[0] = e

    thread = threading.Thread(target=work, daemon=True)
    thread.start()
    with tqdm(ncols=74, unit=" alignments", leave=True,
              bar_format="    Processed: {n:,} alignments [{elapsed}, {rate_fmt}]") as bar:
        while thread.is_alive():
            thread.join(0.2)
            bar.update(prog.reads - bar.n)
        bar.update(prog.reads - bar.n)
        if error[0] is not None:
            bar.leave = False
    if isinstance(error[0], ValueError):
        if "not grouped by reference" in str(error[0]):
            raise NotGroupedByReference()
        return None
    if error[0] is not None:
        raise error[0]

    totals, accs = result[0]
    contig_pid_stats = None if accs is None else {name: _pid_stats_from_native(acc) for name, acc in accs}

    return totals, contig_pid_stats


def _pysam_coverage(bam_path, wanted, exclude_flags):
    """
    Coverage totals per wanted contig from the bam via pysam: what the compiled core's pass
    counts, the same way. Each alignment not flagged with any of exclude_flags adds 1 from
    its start to its end on the contig; a contig's alignments have to be together.
    """
    import numpy as np # type: ignore
    import pysam # type: ignore

    slot = {header: i for i, header in enumerate(wanted)}
    totals = [None] * len(wanted)

    def add(header, total_cov, detected, detected_10x, hist):
        i = slot[header]
        if totals[i] is None:
            totals[i] = (total_cov, detected, detected_10x, hist)
        else:
            t = totals[i]
            totals[i] = (t[0] + total_cov, t[1] + detected, t[2] + detected_10x, t[3] + hist)

    def close(tid, diff):
        depths, bases = np.unique(np.cumsum(diff[:-1], dtype=np.int64), return_counts=True)
        hist = [(int(d), int(n)) for d, n in zip(depths, bases)]
        add(names[tid],
            sum(d * n for d, n in hist), sum(n for d, n in hist if d > 0), sum(n for d, n in hist if d >= 10), hist)

    with pysam.AlignmentFile(bam_path, "rb", threads=min(4, os.cpu_count() or 1)) as bam:

        names, lengths = bam.references, bam.lengths
        seen = [False] * len(names)
        current, diff = None, None

        for read in tqdm(bam.fetch(until_eof=True), ncols=74, unit=" alignments", leave=True,
                         bar_format="    Processed: {n:,} alignments [{elapsed}, {rate_fmt}]"):
            tid = read.reference_id
            if read.flag & exclude_flags or tid < 0 or read.reference_start < 0 or names[tid] not in slot:
                continue
            if tid != current:
                if current is not None:
                    close(current, diff)
                if seen[tid]:
                    raise NotGroupedByReference()
                seen[tid] = True
                current, diff = tid, np.zeros(lengths[tid] + 1, dtype=np.int32)
            start = min(read.reference_start, lengths[tid])
            # (no end without a CIGAR: such an alignment covers its one start base)
            end = min(read.reference_end or start + 1, lengths[tid])
            if start < end:
                diff[start] += 1
                diff[end] -= 1

        if current is not None:
            close(current, diff)

        # contigs with no alignments are all depth 0
        for tid, header in enumerate(names):
            if header in slot and not seen[tid] and lengths[tid] > 0:
                add(header, 0, 0, 0, [(0, lengths[tid])])

    return totals


def scan_bam(bam_path, contigs, header_to_ref_dict, include_non_primary=False, want_pids=True):
    """
    One pass over the bam: adds each contig's per-base coverage to its stats and its ref's,
    and returns {contig: PidStats} of the reads mapped to each (None if not wanted, or the
    bam can't be read for them).

    Coverage is counted as `mosdepth -x` counts it: each alignment adds 1 from its start to
    its end on the contig, deletions and skips in between included, and overlapping mates
    both count.
    """

    print(f"\n  {Fore.YELLOW}Reading the bam for coverage{' and percent IDs' if want_pids else ''}...")

    lookup = stats_lookup(contigs, header_to_ref_dict)
    wanted = list(lookup)
    exclude_flags = COVERAGE_EXCLUDE_FLAGS_WITH_NON_PRIMARY if include_non_primary else COVERAGE_EXCLUDE_FLAGS

    def one_pass(bam_path):
        _check_opens_with_pysam(bam_path)
        done = _native_bam_pass(bam_path, wanted, exclude_flags, include_non_primary, want_pids)
        if done is not None:
            return done
        # the compiled core isn't there, or left this bam to pysam
        totals = _pysam_coverage(bam_path, wanted, exclude_flags)
        contig_pid_stats = None
        if want_pids:
            try:
                contig_pid_stats = get_per_contig_pid_stats(bam_path, include_non_primary)
            except Exception:
                pass
        return totals, contig_pid_stats

    try:
        totals, contig_pid_stats = one_pass(bam_path)
    except NotGroupedByReference:
        report_message("The BAM isn't sorted, so we're sorting it for you. Why? Because it's common courtesy.",
                       color="orange", initial_indent="    ", subsequent_indent="    ")
        totals, contig_pid_stats = one_pass(sort_bam_file(bam_path))

    add_contig_totals(lookup, wanted, totals)

    return contig_pid_stats


def two_decimals(value):
    """ a table value always written with 2 decimals (87.10, 0.00), or NA when there is none """
    return "NA" if value is None else f"{value:.2f}"


def generate_output(refs, output_prefix, ref_mean_pids=None, ref_median_pids=None, ref_mapped_counts=None,
                    contigs=None, contig_mean_pids=None, contig_median_pids=None, contig_mapped_counts=None,
                    skip_per_contig=False):


    primary_out_path = f"{output_prefix}-per-ref.tsv"

    with open(primary_out_path, "w") as f:

        cols = ["ref", "length", "num_contigs", "num_bases_detected", "detection (%)", "detection_at_10x (%)", "mean_coverage", "median_coverage"]
        if ref_mean_pids is not None:
            cols.extend(["mean_pid", "median_pid", "num_mapped_reads"])
        f.write("\t".join(cols) + "\n")

        for ref in refs:
            num_bases_detected, *metrics = ref.stats.compute_metrics()
            row = [ref.path, ref.stats.length, len(ref.headers), num_bases_detected, *map(two_decimals, metrics)]

            if ref_mean_pids is not None:
                mean_pid = ref_mean_pids.get(ref.path)
                median_pid = ref_median_pids.get(ref.path)
                mapped_reads = ref_mapped_counts.get(ref.path, 0)
                row += [two_decimals(mean_pid), two_decimals(median_pid), f"{mapped_reads:,}"]
            f.write("\t".join(map(str, row)) + "\n")

    print(f"\n\n  Ref-level coverage stats written to: {Fore.YELLOW}{primary_out_path}")

    if not skip_per_contig:

        contig_out_path = f"{output_prefix}-per-contig.tsv"

        with open(contig_out_path, "w") as cf:

            cols = ["ref", "contig", "length", "num_bases_detected", "detection (%)", "detection_at_10x (%)", "mean_coverage", "median_coverage"]

            if contig_mean_pids is not None:
                cols.extend(["mean_pid", "median_pid", "num_mapped_reads"])
            cf.write("\t".join(cols) + "\n")

            written_contigs = set()
            for ref in sorted(refs, key=lambda r: r.path):
                for contig_name in ref.contig_order:
                    if contig_name not in contigs:
                        continue
                    contig_data = contigs[contig_name]
                    num_bases_detected, *metrics = contig_data.stats.compute_metrics()
                    row = [contig_data.path, contig_name, contig_data.stats.length, num_bases_detected, *map(two_decimals, metrics)]

                    if contig_mean_pids is not None:
                        mean_pid = contig_mean_pids.get(contig_name)
                        median_pid = contig_median_pids.get(contig_name)
                        mapped_reads = contig_mapped_counts.get(contig_name, 0)
                        row += [two_decimals(mean_pid), two_decimals(median_pid), f"{mapped_reads:,}"]
                    cf.write("\t".join(map(str, row)) + "\n")
                    written_contigs.add(contig_name)

        print(f"  Contig-level coverage stats written to: {Fore.YELLOW}{contig_out_path}")

        print()
