from pathlib import Path
import pytest # type: ignore
from collections import defaultdict
from bit.modules.cov_stats import CoverageStats
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli

test_fasta_path = get_package_path("tests/data/ez-screen-assembly.fasta")
test_bam_path = get_package_path("tests/data/ez-screen.bam")

def test_get_cov_stats(tmp_path):
    out_prefix = tmp_path / "cov-stats"

    cmd = [
        "bit", "cov-stats",
        "-r", str(test_fasta_path),
        "-b", str(test_bam_path),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    ref_summary_tsv = Path(f"{out_prefix}-per-ref.tsv")
    assert ref_summary_tsv.exists(), f"Coverage stats TSV not found at {ref_summary_tsv}"

    observed = ref_summary_tsv.read_text().splitlines()
    expected_header = "ref\tlength\tnum_contigs\tnum_bases_detected\tdetection (%)\tdetection_at_10x (%)\tmean_coverage\tmedian_coverage\tmean_pid\tmedian_pid\tnum_mapped_reads"
    assert observed[0] == expected_header, f"Header does not match expected:\n{observed[0]}\nExpected:\n{expected_header}"

    fields = observed[1].split("\t")
    assert len(fields) == 11, f"Unexpected number of columns:\n{fields}"

    ref, length, num_contigs, num_bases_detected, detection, detection10x, mean_cov, median_cov, mean_pid, median_pid, num_mapped_reads = fields
    assert ref.endswith("ez-screen-assembly.fasta"), f"Unexpected ref path: {ref}"
    assert int(length) == 3010
    assert int(num_contigs) == 1
    assert int(num_bases_detected) == 2623
    assert float(detection) == 87.14
    # always 2 decimals
    assert (detection, detection10x, mean_cov, median_cov, mean_pid, median_pid) == ("87.14", "0.00", "1.99", "2.00", "100.00", "100.00")
    assert float(detection10x) == pytest.approx(0.0, rel=1e-1)
    assert float(mean_cov) == pytest.approx(1.99, rel=1e-1)
    assert float(median_cov) == pytest.approx(2, rel=1e-1)
    assert float(mean_pid) == pytest.approx(100.0, rel=1e-2)
    assert float(median_pid) == pytest.approx(100.0, rel=1e-2)
    assert int(num_mapped_reads.replace(",", "")) == 40

    contig_summary_tsv = Path(f"{out_prefix}-per-contig.tsv")
    assert contig_summary_tsv.exists(), f"Contig-level coverage stats TSV not found at {contig_summary_tsv}"

    observed = contig_summary_tsv.read_text().splitlines()
    expected_header = "ref\tcontig\tlength\tnum_bases_detected\tdetection (%)\tdetection_at_10x (%)\tmean_coverage\tmedian_coverage\tmean_pid\tmedian_pid\tnum_mapped_reads"
    assert observed[0] == expected_header, f"Contig-level TSV header does not match expected:\n{observed[0]}\nExpected:\n{expected_header}"

    fields = observed[1].split("\t")
    assert len(fields) == 11, f"Unexpected number of columns in contig-level tsv:\n{fields}"

    ref, contig, length, num_bases_detected, detection, detection10x, mean_cov, median_cov, mean_pid, median_pid, num_mapped_reads = fields
    assert ref.endswith("ez-screen-assembly.fasta"), f"Unexpected ref path: {ref}"
    assert contig == "partial-NC_003131.1"
    assert int(length) == 3010
    assert int(num_bases_detected) == 2623
    assert float(detection) == 87.14
    # always 2 decimals
    assert (detection, detection10x, mean_cov, median_cov, mean_pid, median_pid) == ("87.14", "0.00", "1.99", "2.00", "100.00", "100.00")
    assert float(detection10x) == pytest.approx(0.0, rel=1e-1)
    assert float(mean_cov) == pytest.approx(1.99, rel=1e-1)
    assert float(median_cov) == pytest.approx(2, rel=1e-1)
    assert float(mean_pid) == pytest.approx(100.0, rel=1e-2)
    assert float(median_pid) == pytest.approx(100.0, rel=1e-2)
    assert int(num_mapped_reads.replace(",", "")) == 40


def _make_stats(length,  hist_dict):
    stats = CoverageStats(length)
    stats.depth_hist = defaultdict(int, hist_dict)

    stats.total_coverage_count = sum(depth * bases for depth, bases in stats.depth_hist.items())
    stats.total_bases_detected_at_all = sum(bases for depth, bases in stats.depth_hist.items() if depth > 0)
    stats.total_bases_detected_at_10x = sum(bases for depth, bases in stats.depth_hist.items() if depth >= 10)

    return stats


def test_median_from_hist_odd_length():

    stats = _make_stats(length = 5, hist_dict = {0: 1, 1: 3, 10: 1})

    # mean = (0 + 1 + 1 + 1 + 10) / 5 = 13/5 = 2.6
    assert round(stats.total_coverage_count / stats.length, 2) == 2.60

    # median = middle element = 1
    assert stats._median_from_hist() == 1.00


def test_median_from_hist_even_length_fractional():
    # 10 bases: 0,0,0,0,0,10,10,10,10,10
    stats = _make_stats(length = 10, hist_dict = {0: 5, 10: 5})

    # mean = 50/10 = 5
    assert round(stats.total_coverage_count / stats.length, 2) == 5.00

    # median indices are 4 and 5, so values 0 and 10, median should be 5
    assert stats._median_from_hist() == 5.00


def test_compute_metrics_returns_expected_mean_and_median():
    # 6 bases: 1,1,1,1,9,9
    stats = _make_stats(length = 6, hist_dict = {1: 4, 9: 2})

    num_bases_detected, detection, detection_at_10x, mean_coverage, median_coverage = stats.compute_metrics()

    # mean = (1 + 1 + 1 + 1 + 9 + 9) / 6 = 22/6 = 3.666... -> 3.67
    assert mean_coverage == 3.67

    assert median_coverage == 1.00

    assert (num_bases_detected, detection, detection_at_10x) == (6, 100.0, 0.0)


def test_compute_metrics_detection_is_a_percent_to_two_decimals():
    # 10,000 bases, 1,598 covered, 337 of those at 10x or more
    stats = _make_stats(length = 10000, hist_dict = {0: 8402, 3: 1261, 12: 337})

    num_bases_detected, detection, detection_at_10x, _, _ = stats.compute_metrics()

    assert num_bases_detected == 1598
    assert detection == 15.98
    assert detection_at_10x == 3.37


def test_update_from_bed_line_accumulates_hist_and_mean():
    stats = CoverageStats(length = 10)

    stats.update_from_bed_line(span = 3, num_reads = 5)

    stats.update_from_bed_line(span = 2, num_reads = 0)

    assert stats.depth_hist[5] == 3
    assert stats.depth_hist[0] == 2
    assert stats.total_coverage_count == 15

    # mean = 15/10 = 1.5
    assert round(stats.total_coverage_count / stats.length, 2) == 1.50


def test_per_base_bed_read_the_same_as_bgzf_or_plain_gzip(tmp_path):
    """ mosdepth's bed is BGZF (read block by block, in parallel); a plain gzip file of the
    same lines gives the same totals """
    import gzip
    import random
    import pysam # type: ignore
    from bit import _native as native

    rng = random.Random(7)
    lines, pos = [], 0
    for contig in ("c1", "c2", "c3"):
        pos = 0
        for _ in range(40000):
            end = pos + rng.randint(1, 30)
            lines.append(f"{contig}\t{pos}\t{end}\t{rng.choice([0, 0, 1, 2, 9, 10, 11, 250])}\n")
            pos = end
    text = "".join(lines)
    plain, bgzf = str(tmp_path / "plain.bed.gz"), str(tmp_path / "bgzf.bed.gz")
    with gzip.open(plain, "wt") as f:
        f.write(text)
    with pysam.BGZFile(bgzf, "wb") as f:
        f.write(text.encode())
    wanted = ["c1", "c3", "absent"]
    a = native.per_base_bed_stats(plain, wanted)
    b = native.per_base_bed_stats(bgzf, wanted)
    norm = lambda r: [None if t is None else (t[0], t[1], t[2], sorted(t[3])) for t in r]
    assert norm(a) == norm(b)
    assert a[2] is None and a[0][0] == sum(
        (int(l.split("\t")[2]) - int(l.split("\t")[1])) * int(l.split("\t")[3]) for l in lines if l.startswith("c1\t"))
