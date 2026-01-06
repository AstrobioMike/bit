import shutil
import pytest
from unittest.mock import patch, MagicMock
from bit.modules.assemble import check_conda_env, RunConfig
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


@patch("bit.modules.assemble.os.path.exists")
@patch("bit.modules.assemble.run")
@patch("bit.modules.assemble.get_package_path")
@patch("bit.modules.assemble.sys.platform", "linux")
def test_check_conda_env_creation_logic(mock_pkg_path, mock_run, mock_exists):

    config = RunConfig()

    # when the environment does NOT exist
    mock_exists.return_value = False

    mock_pkg_path.return_value = "/fake/path/assemble.yaml"

    mock_run.return_value = MagicMock(returncode=0)

    updated_config = check_conda_env(config)

    called_path = mock_exists.call_args[0][0]
    assert "bit-assemble" in called_path

    # checking the conda create command was made correctly
    mock_run.assert_called_once()
    args_passed = mock_run.call_args[0][0]
    assert args_passed[0] == "conda"
    assert "env" in args_passed
    assert "create" in args_passed
    assert "--file" in args_passed
    assert "/fake/path/assemble.yaml" in args_passed

    assert updated_config.conda_env == "bit-assemble"


@patch("bit.modules.assemble.os.path.exists")
@patch("bit.modules.assemble.run")
@patch("bit.modules.assemble.sys.platform", "linux")
def test_check_conda_env_creation_failure(mock_run, mock_exists):

    config = RunConfig()

    # when the conda env does NOT exist
    mock_exists.return_value = False

    mock_run.side_effect = Exception("Conda process died")

    with pytest.raises(SystemExit) as e:
        check_conda_env(config)

    assert e.value.code == 1
