#! /usr/bin/env python

import os
import math
import pysam
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


def get_mapped_reads_pids(input_bam, include_non_primary=False, store_read_pids=True):
    """
    The calculation is alignment-centric:
        PID = (full_aligned_length - NM) / full_aligned_length * 100
        where full_aligned_length = Matches + Mismatches + Insertions + Deletions

    When store_read_pids is False, per-read data is not kept in memory
    and only summary statistics (via PidStats) are tracked.
    """
    decompression_threads = min(4, os.cpu_count() or 1)

    with pysam.AlignmentFile(input_bam, "rb", threads=decompression_threads) as bam:

        ref_name_map = {i: bam.get_reference_name(i) for i in range(bam.nreferences)}

        ref_read_pids = defaultdict(list) if store_read_pids else None
        pid_stats = PidStats()

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
                if op in {0, 1, 2, 7, 8}   # M, I, D, =, X in CIGAR
            )

            pid = (full_aligned_length - nm) / full_aligned_length * 100
            pid_stats.add(pid)

            if store_read_pids:
                read_id = read.query_name
                if read.is_paired:
                    if read.is_read1:
                        read_id += "/1"
                    elif read.is_read2:
                        read_id += "/2"

                ref_read_pids[ref_name_map[read.reference_id]].append((read_id, pid))

        return ref_read_pids, pid_stats


def get_summary_stats(pid_stats):

    summary_stats = [
        ("Num mapped reads:", f"{pid_stats.count:,}"),
        ("Min pid:", pid_stats.min_val),
        ("Max pid:", pid_stats.max_val),
        ("Median pid:", pid_stats.median),
        ("StDev of pids:", pid_stats.stdev),
    ]

    return pid_stats.mean, summary_stats
