import pysam
import numpy as np
from bit.modules.get_mapped_reads_pid import (get_mapped_reads_pids,
                                              get_summary_stats)


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
        a.set_tag("NM", 2) # 2 mismatches
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


def test_get_mapped_reads_pids(tmp_path):

    bam = create_test_bam(tmp_path)
    ref_read_pids, all_pids = get_mapped_reads_pids(bam)

    assert isinstance(ref_read_pids, dict)
    assert "ref" in ref_read_pids
    assert len(ref_read_pids["ref"]) == 2
    assert sorted(all_pids) == [96.0, 100.0]

    read_ids = [read_id for read_id, pid in ref_read_pids["ref"]]
    assert set(read_ids) == {"read1/1", "read1/2"}


def test_get_summary_stats():
    all_pids = [96.0, 100.0]
    mean, summary = get_summary_stats(all_pids)

    assert np.isclose(mean, 98.0)
    assert any("Median pid:" in str(row) for row in summary)
    assert any("StDev of pids:" in str(row) for row in summary)
