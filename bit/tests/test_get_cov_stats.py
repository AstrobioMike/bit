from pathlib import Path
import pytest
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli

test_fasta_path = get_package_path("tests/data/ez-screen-assembly.fasta")
test_bam_path = get_package_path("tests/data/ez-screen.bam")

def test_get_cov_stats(tmp_path):
    out_prefix = tmp_path / "cov-stats"

    cmd = [
        "bit-get-cov-stats",
        "-r", str(test_fasta_path),
        "-b", str(test_bam_path),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}.tsv")
    assert summary_tsv.exists(), f"Coverage stats TSV not found at {summary_tsv}"

    observed = summary_tsv.read_text().splitlines()
    expected_header = "ref\tdetection\tdetection_at_10x\tmean_coverage\tmean_pid\tnum_mapped_reads"
    assert observed[0] == expected_header, f"Header does not match expected:\n{observed[0]}\nExpected:\n{expected_header}"

    fields = observed[1].split("\t")
    assert len(fields) == 6, f"Unexpected number of columns:\n{fields}"

    ref, detection, detection10x, mean_cov, mean_pid, num_mapped_reads = fields
    assert ref.endswith("ez-screen-assembly.fasta"), f"Unexpected ref path: {ref}"
    assert float(detection) == pytest.approx(0.8714, rel=1e-3)
    assert float(detection10x) == pytest.approx(0.0, rel=1e-1)
    assert float(mean_cov) == pytest.approx(1.9934, rel=1e-1)
    assert float(mean_pid) == pytest.approx(100.0, rel=1e-2)
    assert int(num_mapped_reads.replace(",", "")) == 40
