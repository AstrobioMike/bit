import shutil
import gzip
import pytest
from bit.modules.gen_reads import parse_proportions_file
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


def test_parse_proportions_file_logic(tmp_path):

    prop_file = tmp_path / "proportions.tsv"

    # test case if the user provides abundances instead of proportions summing to 1
    prop_file.write_text("file1.fasta\t3\nfile2.fasta\t1\n")

    # total_proportion will be 4.0; file1 should become 0.75, file2 should become 0.25
    observed = parse_proportions_file(str(prop_file), ["file1.fasta", "file2.fasta"])

    assert observed["file1.fasta"] == 0.75
    assert observed["file2.fasta"] == 0.25
    assert len(observed) == 2


def test_parse_proportions_equal_split():

    input_fastas = ["a.fasta", "b.fasta", "c.fasta", "d.fasta"]

    observed = parse_proportions_file(None, input_fastas)

    assert observed["a.fasta"] == 0.25
    assert observed["b.fasta"] == 0.25
    assert observed["c.fasta"] == 0.25
    assert observed["d.fasta"] == 0.25
    assert sum(observed.values()) == pytest.approx(1.0)
