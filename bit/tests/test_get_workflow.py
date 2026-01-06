import pytest
from argparse import Namespace
from unittest.mock import patch, MagicMock
from bit.modules import get_workflow as gwf


@patch('bit.modules.get_workflow.requests.get')
def test_get_all_releases_pagination(mock_get):

    mock_resp_page1 = MagicMock()
    mock_resp_page1.json.return_value = [{"tag_name": "v1.0.0"}, {"tag_name": "v1.1.0"}]
    mock_resp_page1.raise_for_status.return_value = None

    mock_resp_page2 = MagicMock()
    mock_resp_page2.json.return_value = []  # Empty list to stop the loop
    mock_resp_page2.raise_for_status.return_value = None

    mock_get.side_effect = [mock_resp_page1, mock_resp_page2]

    all_releases = gwf.get_all_releases()

    assert len(all_releases) == 2
    assert all_releases[0]["tag_name"] == "v1.0.0"

    assert mock_get.call_count == 2

    first_call_params = mock_get.call_args_list[0][1]['params']
    second_call_params = mock_get.call_args_list[1][1]['params']

    assert first_call_params["page"] == 1
    assert second_call_params["page"] == 2


@patch('bit.modules.get_workflow.get_all_releases')
def test_get_versions_available(mock_get_all):

    mock_get_all.return_value = [
        {'html_url': 'https://github.com/.../releases/tag/metagenomics-wf-v1.0.0'},
        {'html_url': 'https://github.com/.../releases/tag/metagenomics-wf-v1.1.0'},
        {'html_url': 'https://github.com/.../releases/tag/other-wf-v2.0.0'}
    ]

    versions, mapping = gwf.get_versions_available("metagenomics-wf")

    assert "1.0.0" in versions
    assert "1.1.0" in versions
    assert "2.0.0" not in versions
    assert mapping["1.1.0"] == "metagenomics-wf-v1.1.0"


def test_check_version_available_fails():

    available = ["1.0.0", "1.1.0"]
    with pytest.raises(SystemExit) as e:
        gwf.check_version_available("metagenomics", "2.0.0", available)
    assert e.value.code == 1


def test_check_version_available_passes():

    available = ["1.0.0", "1.1.0"]
    gwf.check_version_available("metagenomics", "1.0.0", available)


@patch('os.path.exists')
def test_check_if_dir_exists_exits(mock_exists):

    mock_exists.return_value = True
    with pytest.raises(SystemExit) as e:
        gwf.check_if_dir_already_exists("some_dir")
    assert e.value.code == 1


@patch('os.path.exists')
def test_check_if_dir_not_exists_passes(mock_exists):

    mock_exists.return_value = False
    gwf.check_if_dir_already_exists("some_dir")


@patch('requests.get')
@patch('os.system')
@patch('os.remove')
def test_download_and_unzip(mock_remove, mock_os_system, mock_requests_get):

    mock_resp = MagicMock()
    mock_resp.content = b"fake_zip_binary_data"
    mock_requests_get.return_value = mock_resp

    with patch("builtins.open", MagicMock()):
        gwf.download_and_unzip("metagenomics", "https://fake.link/test.zip")

    mock_os_system.assert_called_with("unzip -q test.zip")
    mock_remove.assert_called_with("test.zip")


@patch('bit.modules.get_workflow.get_versions_available')
def test_dl_wf_list_only(mock_get_versions):

    mock_get_versions.return_value = (["1.0.0"], {"1.0.0": "meta-v1.0.0"})
    args = Namespace(workflow="metagenomics", wanted_version=None, list_available_versions=True)

    with pytest.raises(SystemExit) as e:
        gwf.dl_wf(args)
    assert e.value.code == 0


@patch('bit.modules.get_workflow.get_versions_available')
@patch('bit.modules.get_workflow.check_if_dir_already_exists')
@patch('bit.modules.get_workflow.download_and_unzip')
def test_dl_wf_full_flow(mock_download, mock_check_dir, mock_get_versions):

    mock_get_versions.return_value = (
        ["1.0.0", "1.1.0"],
        {"1.0.0": "metagenomics-wf-v1.0.0", "1.1.0": "metagenomics-wf-v1.1.0"}
    )

    # simulating input args (e.g., bit-get-wf metagenomics)
    args = Namespace(
        workflow="metagenomics",
        wanted_version=None,
        list_available_versions=False
    )

    gwf.dl_wf(args)

    expected_full_name = "metagenomics-wf-v1.1.0"
    expected_link = (
        "https://github.com/astrobiomike/bit/releases/download/"
        f"{expected_full_name}/bit-{expected_full_name}.zip"
    )

    mock_check_dir.assert_called_once_with("bit-" + expected_full_name)
    mock_download.assert_called_once_with("metagenomics", expected_link)
