import shutil
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli

test_R1 = get_package_path("tests/data/mock-R1.fastq.gz")
test_R2 = get_package_path("tests/data/mock-R2.fastq.gz")

def test_assemble(tmp_path):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    shutil.copy(test_R1, reads_dir)
    shutil.copy(test_R2, reads_dir)

    out_dir = tmp_path / "out"

    cmd = [
        "bit-assemble",
        "-r", str(reads_dir),
        "-o", str(out_dir),
        "--run-fastp",
        # "--run-bbnorm",
        "-j", "1",
        "-t", "1",
        "-m", "1e8",
    ]

    res = run_cli(cmd)

    assembly = out_dir / "mock" / "mock-assembly.fasta"
    assert assembly.exists(), f"Assembly file not found at {assembly}"

    summary_tsv = out_dir / "assembly-summaries.tsv"
    assert summary_tsv.exists(), f"Summary TSV not found at {summary_tsv}"
    text = summary_tsv.read_text().splitlines()
    assert text[0].startswith("Assembly\t")
