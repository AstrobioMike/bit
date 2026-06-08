#! /usr/bin/env python

import os
import math
import pysam # type: ignore
from collections import defaultdict


class PidStats:

    __slots__ = ("count", "total", "min_val", "max_val", "hist")

    def __init__(self):
        self.count = 0
        self.total = 0.0
        self.min_val = float("inf")
        self.max_val = float("-inf")
        self.hist = defaultdict(int)

    def add(self, pid):
        self.count += 1
        self.total += pid
        if pid < self.min_val:
            self.min_val = pid
        if pid > self.max_val:
            self.max_val = pid
        self.hist[round(pid, 2)] += 1

    @property
    def mean(self):
        return self.total / self.count if self.count else 0.0

    @property
    def median(self):
        if self.count == 0:
            return 0.0

        i1 = (self.count - 1) // 2
        i2 = self.count // 2

        cumulative = 0
        v1 = v2 = None

        for pid_bin in sorted(self.hist):
            cumulative += self.hist[pid_bin]

            if v1 is None and i1 < cumulative:
                v1 = pid_bin
            if v2 is None and i2 < cumulative:
                v2 = pid_bin

            if v1 is not None and v2 is not None:
                break

        return (v1 + v2) / 2

    @property
    def stdev(self):
        if self.count < 2:
            return 0.0
        mu = self.mean
        variance = sum(count * (val - mu) ** 2 for val, count in self.hist.items()) / self.count
        return math.sqrt(variance)

    def merge(self, other):
        """Merge another PidStats into this one."""
        self.count += other.count
        self.total += other.total
        if other.count > 0:
            if other.min_val < self.min_val:
                self.min_val = other.min_val
            if other.max_val > self.max_val:
                self.max_val = other.max_val
        for val, cnt in other.hist.items():
            self.hist[val] += cnt

    @classmethod
    def from_values(cls, pids):
        """Create a PidStats from an iterable of PID values."""
        stats = cls()
        for pid in pids:
            stats.add(pid)
        return stats


class ClipAggregate:
    """Primary-only running totals for clip/alignment reporting."""

    __slots__ = ("n_reads", "total_aligned", "total_soft", "total_hard",
                 "total_read_len", "clip_frac_sum", "length_stats")

    def __init__(self):
        self.n_reads = 0
        self.total_aligned = 0
        self.total_soft = 0
        self.total_hard = 0
        self.total_read_len = 0
        self.clip_frac_sum = 0.0
        self.length_stats = PidStats()

    def add(self, b):
        self.n_reads += 1
        self.total_aligned += b["aligned_cols"]
        self.total_soft += b["soft_clipped"]
        self.total_hard += b["hard_clipped"]
        self.total_read_len += b["full_read_length"]
        self.length_stats.add(b["full_read_length"])
        if b["full_read_length"]:
            self.clip_frac_sum += b["clipped"] / b["full_read_length"]

    @property
    def mean_clipped_frac(self):
        return self.clip_frac_sum / self.n_reads if self.n_reads else 0.0

    @property
    def total_clipped(self):
        return self.total_soft + self.total_hard

    @property
    def mean_aligned(self):
        return self.total_aligned / self.n_reads if self.n_reads else 0.0

    @property
    def mean_soft(self):
        return self.total_soft / self.n_reads if self.n_reads else 0.0

    @property
    def mean_hard(self):
        return self.total_hard / self.n_reads if self.n_reads else 0.0

    @property
    def mean_clipped(self):
        return self.total_clipped / self.n_reads if self.n_reads else 0.0

    @property
    def mean_clipped_pct(self):
        return self.mean_clipped_frac * 100

    def merge(self, other):
        """Merge another ClipAggregate into this one."""
        self.n_reads += other.n_reads
        self.total_aligned += other.total_aligned
        self.total_soft += other.total_soft
        self.total_hard += other.total_hard
        self.total_read_len += other.total_read_len
        self.clip_frac_sum += other.clip_frac_sum
        self.length_stats.merge(other.length_stats)


def parse_cigar_breakdown(cigartuples):
    """
    Base counts by category from a CIGAR.
    Ops: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7=(=), 8=X

    aligned_cols == full_aligned_length used as the PID denominator (M+I+D).
    full_read_length counts query bases incl. hard-clips (deletions are
    reference-side and excluded).
    """
    counts = defaultdict(int)
    for op, length in cigartuples:
        counts[op] += length

    matches_mismatches = counts[0] + counts[7] + counts[8]
    insertions = counts[1]
    deletions = counts[2]
    soft_clipped = counts[4]
    hard_clipped = counts[5]

    aligned_cols = matches_mismatches + insertions + deletions
    clipped = soft_clipped + hard_clipped
    full_read_length = matches_mismatches + insertions + soft_clipped + hard_clipped

    return {
        "aligned_cols": aligned_cols,
        "soft_clipped": soft_clipped,
        "hard_clipped": hard_clipped,
        "clipped": clipped,
        "full_read_length": full_read_length,
    }


def get_mapped_reads_pids(input_bam, include_non_primary=False, store_read_pids=True):
    """
    The calculation is alignment-centric:
        PID = (full_aligned_length - NM) / full_aligned_length * 100
        where full_aligned_length = Matches + Mismatches + Insertions + Deletions

    When store_read_pids is False, per-read data is not kept in memory
    and only summary statistics (via PidStats / ClipAggregate) are tracked.

    Clip/length metrics in ClipAggregate are gathered on primary alignments
    only, so they stay meaningful even when include_non_primary is set.
    """
    decompression_threads = min(4, os.cpu_count() or 1)

    with pysam.AlignmentFile(input_bam, "rb", threads=decompression_threads) as bam:

        ref_name_map = {i: bam.get_reference_name(i) for i in range(bam.nreferences)}

        ref_read_pids = defaultdict(list) if store_read_pids else None
        pid_stats = PidStats()
        clip_agg = ClipAggregate()

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            is_non_primary = read.is_secondary or read.is_supplementary
            if not include_non_primary and is_non_primary:
                continue
            if not read.has_tag("NM"):
                continue
            nm = read.get_tag("NM")

            b = parse_cigar_breakdown(read.cigartuples)
            full_aligned_length = b["aligned_cols"]

            pid = (full_aligned_length - nm) / full_aligned_length * 100
            pid_stats.add(pid)

            # clip/length metrics are per-read; only meaningful on primary
            if not is_non_primary:
                clip_agg.add(b)

            if store_read_pids:
                read_id = read.query_name
                if read.is_paired:
                    if read.is_read1:
                        read_id += "/1"
                    elif read.is_read2:
                        read_id += "/2"

                clipped_frac = (b["clipped"] / b["full_read_length"]
                                if b["full_read_length"] else 0.0)

                ref_read_pids[ref_name_map[read.reference_id]].append(
                    (read_id, pid, b["aligned_cols"], b["full_read_length"],
                     b["soft_clipped"], b["hard_clipped"], clipped_frac)
                )

        return ref_read_pids, pid_stats, clip_agg


def get_per_contig_pid_stats(input_bam, include_non_primary=False):
    """
    Single-pass BAM scan that builds a PidStats object per contig (reference name).
    Returns a dict of {contig_name: PidStats}. No per-read data is stored.
    Used by get_cov_stats where read IDs are not needed.
    """
    decompression_threads = min(4, os.cpu_count() or 1)

    with pysam.AlignmentFile(input_bam, "rb", threads=decompression_threads) as bam:

        ref_name_map = {i: bam.get_reference_name(i) for i in range(bam.nreferences)}
        contig_pid_stats = defaultdict(PidStats)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if not include_non_primary and (read.is_secondary or read.is_supplementary):
                continue
            if not read.has_tag("NM"):
                continue
            nm = read.get_tag("NM")

            full_aligned_length = sum(
                length for op, length in read.cigartuples
                if op in {0, 1, 2, 7, 8}
            )

            pid = (full_aligned_length - nm) / full_aligned_length * 100
            contig_pid_stats[ref_name_map[read.reference_id]].add(pid)

        return dict(contig_pid_stats)


def get_summary_stats(pid_stats, clip_agg=None):

    summary_stats = [
        ("Num mapped reads:", f"{pid_stats.count:,}"),
    ]

    if clip_agg is not None and clip_agg.n_reads:
        ls = clip_agg.length_stats
        summary_stats += [
            ("", ""),
            ("Mean read-length:", f"{ls.mean:,.2f}"),
            ("Median read-length:", f"{ls.median:,.2f}"),
            ("Min read-length:", f"{ls.min_val:,}"),
            ("Max read-length:", f"{ls.max_val:,}"),
            ("", ""),
            ("Mean aligned-bases:", f"{clip_agg.mean_aligned:,.2f}"),
            ("Mean soft-clipped:", f"{clip_agg.mean_soft:,.2f}"),
            ("Mean hard-clipped:", f"{clip_agg.mean_hard:,.2f}"),
            ("Mean clipped:", f"{clip_agg.mean_clipped:,.2f}"),
            ("Mean clipped-percent:", f"{clip_agg.mean_clipped_pct:.2f}"),
        ]

    summary_stats += [
        ("", ""),
        ("Mean percent ID:", f"{pid_stats.mean:,.2f}"),
        ("Median percent ID:", f"{pid_stats.median:,.2f}"),
        ("Min percent ID:", f"{pid_stats.min_val:,.2f}"),
        ("Max percent ID:", f"{pid_stats.max_val:,.2f}"),
        ("StDev of percent ID:", f"{pid_stats.stdev:,.2f}"),
    ]

    return summary_stats