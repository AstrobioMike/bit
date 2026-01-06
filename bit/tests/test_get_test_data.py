import pytest
import subprocess
from unittest.mock import patch, MagicMock, call
from argparse import Namespace

from bit.modules import get_test_data as gtd

@pytest.fixture
def mock_args():
    return Namespace(datatype="metagenomics")

@patch('bit.modules.get_test_data.subprocess.run')
@patch('bit.modules.get_test_data.os.remove')
@patch('bit.modules.get_test_data.time.sleep')
def test_dl_test_data_success_first_try(mock_sleep, mock_remove, mock_run, mock_args):

    mock_run.return_value = MagicMock(returncode=0)

    gtd.dl_test_data(mock_args)

    assert mock_run.call_args_list[0] == call([
        "curl", "-L", "--connect-timeout", "30", "-o", "test-metagenomics-reads.zip",
        "https://figshare.com/ndownloader/files/46096083"
    ])

    assert mock_run.call_args_list[1] == call(["unzip", "-qo", "test-metagenomics-reads.zip"], check=True)

    mock_remove.assert_called_once_with("test-metagenomics-reads.zip")

@patch('bit.modules.get_test_data.subprocess.run')
@patch('bit.modules.get_test_data.time.sleep')
def test_dl_test_data_retry_logic(mock_sleep, mock_run, mock_args):

    mock_run.side_effect = [
        MagicMock(returncode=1), # try 1
        MagicMock(returncode=1), # try 2
        MagicMock(returncode=0), # try 3
        MagicMock(returncode=0)  # unzip call
    ]

    with patch('bit.modules.get_test_data.os.remove'):
        gtd.dl_test_data(mock_args)

    assert mock_run.call_count == 4


@patch('bit.modules.get_test_data.subprocess.run')
@patch('bit.modules.get_test_data.time.sleep')
def test_dl_test_data_all_retries_fail(mock_sleep, mock_run, mock_args):

    mock_run.return_value = MagicMock(returncode=1)

    with pytest.raises(SystemExit) as e:
        gtd.dl_test_data(mock_args)

    assert e.value.code == 1
    assert mock_run.call_count == 3


@patch('bit.modules.get_test_data.subprocess.run')
@patch('bit.modules.get_test_data.os.remove')
def test_dl_test_data_unzip_fail(mock_remove, mock_run, mock_args):

    mock_run.side_effect = [
        MagicMock(returncode=0),
        subprocess.CalledProcessError(returncode=1, cmd="unzip")
    ]

    with pytest.raises(SystemExit) as e:
        gtd.dl_test_data(mock_args)

    assert e.value.code == 1
