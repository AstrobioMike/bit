from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
from pathlib import Path
import shutil

test_targets_fasta = get_package_path("tests/data/ez-screen-targets.fasta")
test_assembly_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")

def test_ez_screen_assembly(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    cmd = [
        "bit-ez-screen",
        "assembly",
        "-a", str(test_assembly_fasta),
        "-t", str(test_targets_fasta),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-assembly-summary.tsv")

    observed = summary_tsv.read_text().splitlines()
    expected = [
        "input-assembly\tyopE\tyopK",
        "ez-screen-assembly\t1\t0"
    ]

    assert observed == expected, f"Summary TSV content does not match expected:\n{observed}\nExpected:\n{expected}"

def test_ez_screen_reads(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    R1 = get_package_path("tests/data/ez-screen-R1.fastq.gz")
    R2 = get_package_path("tests/data/ez-screen-R2.fastq.gz")

    shutil.copy(R1, reads_dir)
    shutil.copy(R2, reads_dir)

    cmd = [
        "bit-ez-screen",
        "reads",
        "-t", str(test_targets_fasta),
        "-r", str(reads_dir),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-reads-summary.tsv")

    observed = summary_tsv.read_text().splitlines()
    expected = [
        "sample\ttarget\tnum_reads_recruited\tdetection\tmean_perc_id",
        "ez-screen\tyopE\t20\t0.9\t100.0"
    ]

    assert observed == expected, f"Summary TSV content does not match expected:\n{observed}\nExpected:\n{expected}"
