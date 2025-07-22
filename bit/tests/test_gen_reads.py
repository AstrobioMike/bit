import shutil
import gzip
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli

test_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")

def test_gen_reads(tmp_path):
    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit-gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "perfect-reads"),
        "-n", "1",
        "-r", "10",
        "-s", "1"
    ]

    run_cli(cmd)

    R1_path = tmp_path / "perfect-reads_R1.fastq.gz"
    R2_path = tmp_path / "perfect-reads_R2.fastq.gz"
    assert R1_path.exists(), f"R1 file not found at {R1_path}"
    assert R2_path.exists(), f"R2 file not found at {R2_path}"

    with gzip.open(R1_path, 'rt') as f:
        expected = f.read().splitlines()

    observed = [
        "@partial-NC_003131.1_550/1",
        "GTGGACGACT",
        "+",
        "IIIIIIIIII"
    ]
    assert expected == observed, f"R1 content does not match expected:\n{observed}\nExpected:\n{expected}"
