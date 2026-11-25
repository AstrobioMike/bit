from pathlib import Path
import pytest # type: ignore
from collections import defaultdict
from bit.modules.cov_stats import CoverageStats
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli, run_cli_here

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
    # coverage comes from the bam itself: no mosdepth run, so none of its files, and no index made
    assert not Path(f"{out_prefix}-mosdepth-files").exists()

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


def _write_bam(path, alignments, refs=(("c1", 100), ("c2", 50), ("c3", 30))):
    """ a bam of (contig index, start, cigar, flag[, NM]) alignments, in the order given (NM 0 if not given; None: no NM tag) """
    import pysam # type: ignore

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": length} for name, length in refs]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, (tid, start, cigar, flag, *nm) in enumerate(alignments):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = f"r{i}"
            a.flag = flag
            a.reference_id = tid
            a.reference_start = start
            a.cigarstring = cigar
            qlen = a.infer_query_length() or 0
            a.query_sequence = "A" * qlen
            a.mapping_quality = 60
            if (nm or [0])[0] is not None:
                a.set_tag("NM", (nm or [0])[0])
            bam.write(a)
    return str(path)


def _norm(totals):
    return [None if t is None else (t[0], t[1], t[2], sorted(t[3])) for t in totals]


# c1: 20M at 10, 10M5D10M at 20 (spans 20-45, the deletion included), a supplementary 10M
# at 80, a secondary 10M at 0, a duplicate and a QC-fail 10M at 50; c2: none; c3: 10M at 25
# (runs 5 bases off the end of its 30)
_ALIGNMENTS = [(0, 0, "10M", 0x100), (0, 10, "20M", 0), (0, 20, "10M5D10M", 16), (0, 50, "10M", 0x400),
               (0, 50, "10M", 0x200), (0, 80, "10M", 0x800), (2, 25, "10M", 0)]


def test_bam_coverage_counts_alignment_spans_as_mosdepth_x_does(tmp_path):
    from bit import _native as native
    from bit.modules.cov_stats import _pysam_coverage

    bam = _write_bam(tmp_path / "t.bam", _ALIGNMENTS)
    wanted = ["c1", "c2", "c3", "absent"]

    # primary and supplementary alignments: depth 1 at 10-20, 2 at 20-30, 1 at 30-45 and 80-90
    # (55 in all over 45 bases)
    expected = [(55, 45, 0, [(0, 55), (1, 35), (2, 10)]), (0, 0, 0, [(0, 50)]), (5, 5, 0, [(0, 25), (1, 5)]), None]
    totals, pids = native.bam_coverage_and_pid_stats(bam, wanted, 1796, False, False, True)
    assert _norm(totals) == expected and pids is None
    assert _norm(_pysam_coverage(bam, wanted, 1796)) == expected

    # with the secondary one too: depth 1 at 0-10 as well
    expected[0] = (65, 55, 0, [(0, 45), (1, 45), (2, 10)])
    totals, _ = native.bam_coverage_and_pid_stats(bam, wanted, 1540, False, True, True)
    assert _norm(totals) == expected
    assert _norm(_pysam_coverage(bam, wanted, 1540)) == expected


def test_bam_coverage_and_pids_from_the_one_pass_match_the_pid_pass(tmp_path):
    from bit import _native as native

    for include_non_primary in (False, True):
        _, pids = native.bam_coverage_and_pid_stats(str(test_bam_path), [], 1796, True, include_non_primary, True)
        assert pids == native.bam_per_contig_pid_stats(str(test_bam_path), include_non_primary, True)


def test_bam_with_a_contigs_alignments_apart_is_not_read_for_coverage(tmp_path):
    from bit import _native as native
    from bit.modules.cov_stats import _pysam_coverage, NotGroupedByReference

    bam = _write_bam(tmp_path / "apart.bam", [(0, 10, "20M", 0), (2, 5, "10M", 0), (0, 40, "20M", 0)])
    with pytest.raises(ValueError, match="not grouped by reference"):
        native.bam_coverage_and_pid_stats(bam, ["c1", "c3"], 1796, False, False, True)
    with pytest.raises(NotGroupedByReference):
        _pysam_coverage(bam, ["c1", "c3"], 1796)
    # (alignments of one contig in any order are fine, and so are contigs that aren't wanted)
    totals, _ = native.bam_coverage_and_pid_stats(bam, ["c3"], 1796, False, False, True)
    assert _norm(totals) == [(10, 10, 0, [(0, 20), (1, 10)])]


def test_get_cov_stats_without_the_compiled_core_gives_the_same_tables(tmp_path, monkeypatch):
    """ the pysam path, taken when the compiled core isn't there or leaves a bam to it """
    import argparse
    import bit.modules.cov_stats as cov_stats

    def tables(prefix):
        args = argparse.Namespace(reference_fastas=[str(test_fasta_path)], bam=str(test_bam_path), bed=None,
                                  output_prefix=str(prefix), skip_per_contig=False, include_non_primary=False,
                                  skip_read_pids=False, min_read_pid=0)
        cov_stats.get_cov_stats(args)
        return Path(f"{prefix}-per-ref.tsv").read_text(), Path(f"{prefix}-per-contig.tsv").read_text()

    native = tables(tmp_path / "native")
    monkeypatch.setattr(cov_stats, "_native_bam_pass", lambda *a, **k: None)
    assert tables(tmp_path / "pysam") == native


def test_alignments_under_the_minimum_percent_id_count_for_nothing(tmp_path):
    """ --min-read-pid: in the compiled pass and the pysam one alike, alignments under it (or without NM) add to
    neither the coverage nor the percent IDs """
    from bit import _native as native
    from bit.modules.cov_stats import _pysam_coverage
    from bit.modules.mapped_read_stats import get_per_contig_pid_stats

    # on c1, 20M at 10 with 0 mismatches (100%), 20M at 40 with 2 (90%), 20M at 70 with no NM tag
    bam = _write_bam(tmp_path / "pid.bam", [(0, 10, "20M", 0, 0), (0, 40, "20M", 0, 2), (0, 70, "20M", 0, None)])

    everything = (60, 60, 0, [(0, 40), (1, 60)])
    kept = (20, 20, 0, [(0, 80), (1, 20)])
    for min_pid, expected, pids in ((0, everything, [100.0, 90.0]), (95, kept, [100.0]), (90, None, [100.0, 90.0])):
        totals, accs = native.bam_coverage_and_pid_stats(bam, ["c1"], 1796, True, False, True, min_pid=min_pid)
        python = _pysam_coverage(bam, ["c1"], 1796, min_pid)
        if expected:
            assert _norm(totals) == [expected] == _norm(python)
        else:
            assert _norm(totals) == _norm(python)
        assert get_per_contig_pid_stats(bam, min_pid=min_pid)["c1"].count == len(pids)
        assert dict(accs)["c1"][0] == len(pids)


def test_existing_outputs_need_force_which_replaces_only_them(tmp_path):
    cmd = ["bit", "cov-stats", "-r", str(test_fasta_path), "-b", str(test_bam_path)]
    (tmp_path / "keep.txt").write_text("mine")
    # the default prefix, in the current directory
    run_cli_here(cmd + ["--skip-per-contig"], cwd=tmp_path)
    (tmp_path / "cov-stats-per-contig.tsv").write_text("from an older run")
    again = run_cli_here(cmd, cwd=tmp_path, check=False)
    assert again.returncode != 0 and "already exist" in again.stdout
    assert (tmp_path / "cov-stats-per-contig.tsv").read_text() == "from an older run"
    run_cli_here(cmd + ["--skip-per-contig", "-F"], cwd=tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["cov-stats-per-ref.tsv", "keep.txt"]
    assert (tmp_path / "keep.txt").read_text() == "mine"
