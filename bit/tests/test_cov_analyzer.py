import os
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import pytest
from bit.modules import cov_analyzer as ca


def _write_fasta(path: Path, name: str, length: int, base: str = "A"):
    seq = (base * length).encode()
    with open(path, "wb") as fh:
        fh.write(f">{name}\n".encode())
        # wrapping at 60 so faidx is happy
        for i in range(0, length, 60):
            fh.write(seq[i:i+60] + b"\n")
    # indexing
    pysam.faidx(str(path))


def _make_bam(path: Path, contig: str, length: int):
    """
    Creates a coordinate-sorted BAM with reads only in 100-200 on the contig.
    Nothing maps around 900-1000+, so zero-coverage is detectable there.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": length}],
        "PG": [{"ID": "pytest", "PN": "pytest"}],
    }

    # writing unsorted BAM first
    unsorted = str(path.with_suffix(".unsorted.bam"))
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        # 20 reads of length 100 covering 100-200-ish
        for i, start in enumerate(range(100, 200, 5)):
            a = pysam.AlignedSegment()
            a.query_name = f"read{i}"
            a.query_sequence = "A" * 100
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = ((0, 100),)
            a.next_reference_id = -1
            a.next_reference_start = -1
            a.template_length = 0
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(a)

    # sorting and indexing
    sorted_bam = str(path)
    pysam.sort("-o", sorted_bam, unsorted)
    os.remove(unsorted)
    pysam.index(sorted_bam)


@pytest.mark.parametrize("per_contig", [False, True])
def test_run_cov_analyzer_end_to_end(tmp_path, monkeypatch, per_contig):
    """
    Runs run_cov_analyzer with:
      - sliding windows of 50bp, step 50bp
      - synthetic mosdepth regions:
          * most windows: cov=1
          * windows 100-200: cov=10 (two adjacent windows -> one merged HIGH region)
          * windows 900-1000: cov=0 (two adjacent windows -> one merged LOW region)
      - BAM has reads only in 100-200; the 900-1000 region has zero coverage.

    Asserts:
      - output files exist
      - one high and one low merged region
      - low region reports zero_cov_bases equal to region length
    """
    outdir = tmp_path / "out"
    (outdir / "mosdepth-files").mkdir(parents=True)

    # building reference and BAM
    ref = tmp_path / "ref.fa"
    _write_fasta(ref, "chr1", 2000, "A")
    bam = tmp_path / "reads.bam"
    _make_bam(bam, "chr1", 2000)

    # monkey patching noisy/loggy helpers to be no-ops for test cleanliness
    monkeypatch.setattr(ca, "tee", lambda *a, **k: None)
    monkeypatch.setattr(ca, "log_command_run", lambda *a, **k: None)

    # monkey patching run_mosdepth
    def fake_run_mosdepth(bam_file, window_bed_path, output_dir):
        bam_stem = Path(bam_file).stem
        mosdepth_out_prefix = Path(output_dir) / "mosdepth-files" / bam_stem
        regions_path = mosdepth_out_prefix.with_suffix(".regions.bed.gz")

        regions_path.parent.mkdir(parents=True, exist_ok=True)

        bed = pd.read_csv(
            window_bed_path,
            sep="\t",
            header=None,
            names=["contig", "start", "end"]
        )

        cov = np.ones(len(bed), dtype=float)
        # two high-cov windows: 100–150 and 150–200
        high_mask = (bed["start"] >= 100) & (bed["end"] <= 200)
        cov[high_mask] = 10.0

        # two low windows: 900–950 and 950–1000
        low_mask = (bed["start"] >= 900) & (bed["end"] <= 1000)
        cov[low_mask] = 0.0

        df = bed.copy()
        df["cov"] = cov

        df.to_csv(
            regions_path,
            sep="\t",
            header=False,
            index=False,
            compression="gzip",
        )
        return str(regions_path)

    monkeypatch.setattr(ca, "run_mosdepth", fake_run_mosdepth)

    ca.run_cov_analyzer(
        reference_fasta=str(ref),
        bam_file=str(bam),
        output_dir=str(outdir),
        high_threshold=2.0,
        low_threshold=2.0,
        min_region_length=0,
        exclude_contigs=[],
        per_contig=per_contig,
        sliding_window_size=50,
        step_size=50,
        allowed_gap=0,
        buffer=10,
        no_window_stats=False,
        log_file=None,
        full_cmd_executed="pytest-run",
    )

    # expected outputs
    expected_files = [
        outdir / "window-coverage-overview.txt",
        outdir / "window-coverage-overview.tsv",
        outdir / "window-coverage-stats.tsv.gz",
        outdir / "window-coverage-histogram.png",
        outdir / "high-coverage-regions-of-interest.tsv",
        outdir / "high-coverage-regions-of-interest.fasta",
        outdir / "low-coverage-regions-of-interest.tsv",
        outdir / "low-coverage-regions-of-interest.fasta",
    ]
    for p in expected_files:
        assert p.is_file(), f"Missing expected output: {p}"

    # validating merged high region (should have merge 100-150 and 150-200)
    high_df = pd.read_csv(outdir / "high-coverage-regions-of-interest.tsv", sep="\t")
    # should be one merged region
    assert len(high_df) == 1
    r = high_df.iloc[0]
    assert r["contig"] == "chr1"
    assert int(r["start"]) == 100 and int(r["end"]) == 200
    assert int(r["length"]) == 100
    # should be a high mean coverage (close to 10)
    assert r["cov"] >= 9.0

    # validating merged low region (should merge 900-950 and 950-1000)
    low_df = pd.read_csv(outdir / "low-coverage-regions-of-interest.tsv", sep="\t")
    assert len(low_df) == 1
    r = low_df.iloc[0]
    assert r["contig"] == "chr1"
    assert int(r["start"]) == 900 and int(r["end"]) == 1000
    assert int(r["length"]) == 100
    assert r["cov"] == 0.0
    assert int(r["zero_cov_bases"]) == int(r["length"])

    # just checking histogram & overview tsv aren't empty
    assert (outdir / "window-coverage-histogram.png").stat().st_size > 0
    overview = pd.read_csv(outdir / "window-coverage-overview.tsv", sep="\t")
    assert {"name", "mean", "median", "std", "num_windows"}.issubset(overview.columns)
