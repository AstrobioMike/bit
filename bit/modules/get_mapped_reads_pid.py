#! /usr/bin/env python

import pysam
from collections import defaultdict
import numpy as np


def get_mapped_reads_pids(input_bam):
    """
    The calculation is alignment-centric:
        PID = (full_aligned_length - NM) / full_aligned_length * 100
        where full_aligned_length = Matches + Mismatches + Insertions + Deletions
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam:

        ref_read_pids = defaultdict(list)
        all_pids = []

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            try:
                nm = read.get_tag("NM")
            except KeyError:
                continue

            read_id = read.query_name
            if read.is_paired:
                if read.is_read1:
                    read_id += "/1"
                elif read.is_read2:
                    read_id += "/2"

            full_aligned_length = sum(
                length for op, length in read.cigartuples
                if op in {0, 1, 2, 7, 8}   # M, I, D, =, X in CIGAR
            )

            pid = (full_aligned_length - nm) / full_aligned_length * 100

            refname = bam.get_reference_name(read.reference_id)

            ref_read_pids[refname].append((read_id, pid))
            all_pids.append(pid)

        return ref_read_pids, all_pids


def get_summary_stats(all_pids):

    pid_array = np.array(all_pids)
    mean = np.mean(pid_array)

    summary_stats = [
        ("Num mapped reads:", f"{len(pid_array):,}"),
        ("Min pid:", np.min(pid_array)),
        ("Max pid:", np.max(pid_array)),
        ("Median pid:", np.median(pid_array)),
        ("StDev of pids:", np.std(pid_array)),
    ]

    return mean, summary_stats
