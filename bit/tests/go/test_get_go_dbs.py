import pytest
from pathlib import Path
from unittest.mock import patch, call

from bit.modules.go.get_go_dbs import (
    check_go_data_location_var_is_set,
    check_if_data_present,
    download_go_data,
    get_go_data,
)


def _mock_download(link, label, dest):
    """Side effect for download_with_tqdm that writes fake obo content to dest."""
    Path(dest).write_text(f"[Term]\nid: GO:0000001\n")


################################################################################
# check_go_data_location_var_is_set
################################################################################

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("GO_DB_DIR", str(tmp_path))
    result = check_go_data_location_var_is_set()
    assert result == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("GO_DB_DIR", raising=False)
    with pytest.raises(SystemExit):
        check_go_data_location_var_is_set()


################################################################################
# check_if_data_present
################################################################################

def test_check_if_data_present_both_files_nonempty(tmp_path):
    (tmp_path / "go-basic.obo").write_text("data")
    (tmp_path / "goslim_metagenomics.obo").write_text("data")
    assert check_if_data_present(str(tmp_path)) is True


def test_check_if_data_present_go_basic_missing(tmp_path):
    (tmp_path / "goslim_metagenomics.obo").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "goslim_metagenomics.obo").exists()


def test_check_if_data_present_goslim_missing(tmp_path):
    (tmp_path / "go-basic.obo").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "go-basic.obo").exists()


def test_check_if_data_present_both_missing(tmp_path):
    assert check_if_data_present(str(tmp_path)) is False


def test_check_if_data_present_empty_files_removed(tmp_path):
    (tmp_path / "go-basic.obo").write_text("")
    (tmp_path / "goslim_metagenomics.obo").write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "go-basic.obo").exists()
    assert not (tmp_path / "goslim_metagenomics.obo").exists()


################################################################################
# download_go_data
################################################################################

def test_download_go_data_writes_both_files(tmp_path):
    with patch("bit.modules.go.get_go_dbs.download_with_tqdm", side_effect=_mock_download):
        download_go_data(str(tmp_path))
    assert (tmp_path / "go-basic.obo").exists()
    assert (tmp_path / "goslim_metagenomics.obo").exists()


def test_download_go_data_downloads_correct_urls(tmp_path):
    with patch("bit.modules.go.get_go_dbs.download_with_tqdm", side_effect=_mock_download) as mock_dl:
        download_go_data(str(tmp_path))
    urls = [c.args[0] for c in mock_dl.call_args_list]
    assert any("go-basic.obo" in u for u in urls)
    assert any("goslim_metagenomics" in u for u in urls)


def test_download_go_data_failure_calls_notify_premature_exit(tmp_path):
    with patch("bit.modules.go.get_go_dbs.download_with_tqdm",
               side_effect=Exception("network error")), \
         patch("bit.modules.go.get_go_dbs.notify_premature_exit") as mock_exit:
        mock_exit.side_effect = SystemExit(1)
        with pytest.raises(SystemExit):
            download_go_data(str(tmp_path))
    mock_exit.assert_called_once()


################################################################################
# get_go_data
################################################################################

def test_get_go_data_already_present_skips_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GO_DB_DIR", str(tmp_path))
    (tmp_path / "go-basic.obo").write_text("data")
    (tmp_path / "goslim_metagenomics.obo").write_text("data")
    with patch("bit.modules.go.get_go_dbs.download_go_data") as mock_dl:
        get_go_data(force_update=False, quiet=True)
    mock_dl.assert_not_called()


def test_get_go_data_force_update_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GO_DB_DIR", str(tmp_path))
    (tmp_path / "go-basic.obo").write_text("data")
    (tmp_path / "goslim_metagenomics.obo").write_text("data")
    with patch("bit.modules.go.get_go_dbs.download_go_data") as mock_dl:
        get_go_data(force_update=True)
    mock_dl.assert_called_once_with(str(tmp_path))


def test_get_go_data_missing_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GO_DB_DIR", str(tmp_path))
    with patch("bit.modules.go.get_go_dbs.download_go_data") as mock_dl:
        get_go_data()
    mock_dl.assert_called_once_with(str(tmp_path))
