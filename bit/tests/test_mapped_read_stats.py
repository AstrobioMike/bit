import pysam # type: ignore
from bit.modules.mapped_read_stats import (get_mapped_reads_pids,
                                              get_summary_stats,
                                              PidStats,
                                              ClipAggregate)


def create_test_bam(tmp_path):

    bam_path = tmp_path / "test.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'ref'}]
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:

        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 50
        a.flag = 99  # paired, read1, properly paired
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 50)]  # 50M
        a.set_tag("NM", 2)  # 2 mismatches
        outf.write(a)

        b = pysam.AlignedSegment()
        b.query_name = "read1"
        b.query_sequence = "A" * 50
        b.flag = 147  # paired, read2, properly paired, reverse strand
        b.reference_id = 0
        b.reference_start = 100
        b.mapping_quality = 60
        b.cigar = [(0, 50)]  # 50M
        b.set_tag("NM", 0)
        outf.write(b)

        # this one has no NM tag and should be skipped
        c = pysam.AlignedSegment()
        c.query_name = "read_no_nm"
        c.query_sequence = "A" * 50
        c.flag = 0
        c.reference_id = 0
        c.reference_start = 200
        c.mapping_quality = 60
        c.cigar = [(0, 50)]
        outf.write(c)

    return bam_path


def create_clip_indel_bam(tmp_path):
    """
    Single primary read with a soft-clip, an insertion event, and a deletion
    event so gap-aware and gap-compressed PID diverge and clip columns are non-zero.

    CIGAR: 5S 40M 3I 20M 2D 10M   ; NM = 6 (1 mismatch + 3 ins bases + 2 del bases)
      aligned_cols (M+I+D)        = 70 + 3 + 2 = 75
      matches                     = 75 - 6      = 69
      gap-aware PID  = 69/75*100  = 92.0
      gap_events                  = 1 ins + 1 del = 2
      gap-compr PID  = 69/(70+2)*100 = 95.8333...
      query_aligned (M+I)         = 73
      full_read_length (+soft)    = 78
    """
    bam_path = tmp_path / "clip_indel.bam"
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'ref'}]}

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        a = pysam.AlignedSegment()
        a.query_name = "read_ci"
        a.query_sequence = "A" * 78  # query bases consumed: 5S + 40M + 3I + 20M + 10M
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(4, 5), (0, 40), (1, 3), (0, 20), (2, 2), (0, 10)]
        a.set_tag("NM", 6)
        outf.write(a)

    return bam_path


def test_get_mapped_reads_pids(tmp_path):

    bam = create_test_bam(tmp_path)
    ref_read_pids, pid_stats, pid_gc_stats, clip_agg = get_mapped_reads_pids(bam)

    assert isinstance(ref_read_pids, dict)
    assert "ref" in ref_read_pids
    assert len(ref_read_pids["ref"]) == 2

    assert isinstance(pid_stats, PidStats)
    assert pid_stats.count == 2
    assert pid_stats.min_val == 96.0
    assert pid_stats.max_val == 100.0

    # no indels here, so gap-compressed == gap-aware
    assert isinstance(pid_gc_stats, PidStats)
    assert pid_gc_stats.count == 2
    assert pid_gc_stats.min_val == 96.0
    assert pid_gc_stats.max_val == 100.0

    # both reads primary, unclipped, 50 query-aligned bases each
    assert isinstance(clip_agg, ClipAggregate)
    assert clip_agg.n_reads == 2
    assert clip_agg.total_soft == 0
    assert clip_agg.total_hard == 0
    assert clip_agg.mean_aligned == 50.0

    # per-read tuple is now 10 fields:
    # read_id, pid, gc_pid, q_aln, read_len, soft, hard, clipped_frac, ref_start, ref_end
    assert all(len(row) == 10 for row in ref_read_pids["ref"])

    rows_by_id = {row[0]: row for row in ref_read_pids["ref"]}
    assert set(rows_by_id) == {"read1/1", "read1/2"}

    assert rows_by_id["read1/1"][8:] == (1, 50)
    assert rows_by_id["read1/2"][8:] == (101, 150)


def test_clip_indel_metrics(tmp_path):

    bam = create_clip_indel_bam(tmp_path)
    ref_read_pids, pid_stats, pid_gc_stats, clip_agg = get_mapped_reads_pids(bam)

    (
        read_id,
        pid,
        gc_pid,
        q_aln,
        read_len,
        soft,
        hard,
        clipped_frac,
        ref_start,
        ref_end,
    ) = ref_read_pids["ref"][0]

    assert read_id == "read_ci"
    assert q_aln == 73          # M(70) + I(3)
    assert soft == 5
    assert hard == 0
    assert read_len == 78
    assert q_aln + soft + hard == read_len

    assert ref_start == 1
    assert ref_end == 72        # M(70) + D(2), 1-based inclusive coordinate

    assert abs(pid - 92.0) < 1e-9
    assert abs(gc_pid - 95.8333333) < 1e-6
    assert abs(clipped_frac - 5 / 78) < 1e-9

    assert pid_stats.count == 1
    assert pid_gc_stats.count == 1

    assert clip_agg.n_reads == 1
    assert clip_agg.mean_aligned == 73.0
    assert clip_agg.total_soft == 5


def test_get_summary_stats():

    pid_stats = PidStats()
    pid_stats.add(96.0)
    pid_stats.add(100.0)

    summary = get_summary_stats(pid_stats)

    # get_summary_stats now returns only the list; mean lives on the stats object
    assert pid_stats.mean == 98.0
    assert any("Mean gap-aware PID:" in str(row) for row in summary)
    assert any("Median gap-aware PID:" in str(row) for row in summary)
    assert any("StDev gap-aware PID:" in str(row) for row in summary)


def test_summary_count_split_with_non_primary():

    # 3 alignments enter the PID stats, but only 2 primary reads in clip_agg
    pid_stats = PidStats()
    for v in (99.0, 98.0, 97.0):
        pid_stats.add(v)

    clip_agg = ClipAggregate()
    b = {"query_aligned": 50, "soft_clipped": 0, "hard_clipped": 0,
         "clipped": 0, "full_read_length": 50}
    clip_agg.add(b)
    clip_agg.add(b)

    summary = get_summary_stats(pid_stats, clip_agg=clip_agg, include_non_primary=True)
    text = " ".join(str(row) for row in summary)

    assert "Num mapped reads:" in text     # 2 (primary count)
    assert "Num alignments:" in text        # 3 (all alignments)
