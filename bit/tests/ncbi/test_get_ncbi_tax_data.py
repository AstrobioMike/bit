import pytest # type: ignore
from unittest.mock import patch, MagicMock

from bit.modules.ncbi.get_ncbi_tax_data import (
    check_tax_location_var_is_set,
    check_if_data_present,
    download_ncbi_tax_data,
    get_ncbi_tax_data,
)


def test_check_tax_location_var_is_set_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("TAXONKIT_DB", str(tmp_path))
    result = check_tax_location_var_is_set()
    assert result == str(tmp_path)


def test_check_tax_location_var_is_set_exits_if_missing(monkeypatch):
    monkeypatch.delenv("TAXONKIT_DB", raising=False)
    with pytest.raises(SystemExit):
        check_tax_location_var_is_set()


def test_check_if_data_present_both_files_nonempty(tmp_path):
    (tmp_path / "names.dmp").write_text("data")
    (tmp_path / "nodes.dmp").write_text("data")
    assert check_if_data_present(str(tmp_path)) is True


def test_check_if_data_present_names_missing(tmp_path):
    (tmp_path / "nodes.dmp").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "nodes.dmp").exists()


def test_check_if_data_present_nodes_missing(tmp_path):
    (tmp_path / "names.dmp").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "names.dmp").exists()


def test_check_if_data_present_both_missing(tmp_path):
    assert check_if_data_present(str(tmp_path)) is False


def test_check_if_data_present_empty_files_removed(tmp_path):
    (tmp_path / "names.dmp").write_text("")
    (tmp_path / "nodes.dmp").write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "names.dmp").exists()
    assert not (tmp_path / "nodes.dmp").exists()


def test_check_if_data_present_one_empty_cleans_both(tmp_path):
    (tmp_path / "names.dmp").write_text("data")
    (tmp_path / "nodes.dmp").write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "names.dmp").exists()


def test_download_ncbi_tax_data_success(tmp_path):
    mock_tarball = MagicMock()
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_with_tqdm") as mock_dl, \
         patch("bit.modules.ncbi.get_ncbi_tax_data.tarfile.open", return_value=mock_tarball):
        mock_tarball.__enter__ = lambda s: s
        mock_tarball.__exit__ = MagicMock(return_value=False)
        # create the file so os.remove doesn't fail
        (tmp_path / "taxdump.tar.gz").write_bytes(b"")
        download_ncbi_tax_data(str(tmp_path))
    mock_dl.assert_called_once()
    mock_tarball.extractall.assert_called_once_with(str(tmp_path))


def test_download_ncbi_tax_data_download_failure(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_with_tqdm",
               side_effect=Exception("network error")), \
         patch("bit.modules.ncbi.get_ncbi_tax_data.notify_premature_exit") as mock_exit:
        mock_exit.side_effect = SystemExit(1)
        with pytest.raises(SystemExit):
            download_ncbi_tax_data(str(tmp_path))
    mock_exit.assert_called_once()


def test_get_ncbi_tax_data_already_present_not_quiet(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("TAXONKIT_DB", str(tmp_path))
    (tmp_path / "names.dmp").write_text("data")
    (tmp_path / "nodes.dmp").write_text("data")
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_ncbi_tax_data") as mock_dl:
        get_ncbi_tax_data(force_update=False, quiet=False)
    mock_dl.assert_not_called()


def test_get_ncbi_tax_data_already_present_quiet(tmp_path, monkeypatch):
    monkeypatch.setenv("TAXONKIT_DB", str(tmp_path))
    (tmp_path / "names.dmp").write_text("data")
    (tmp_path / "nodes.dmp").write_text("data")
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_ncbi_tax_data") as mock_dl:
        get_ncbi_tax_data(force_update=False, quiet=True)
    mock_dl.assert_not_called()


def test_get_ncbi_tax_data_force_update_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("TAXONKIT_DB", str(tmp_path))
    (tmp_path / "names.dmp").write_text("data")
    (tmp_path / "nodes.dmp").write_text("data")
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_ncbi_tax_data") as mock_dl:
        get_ncbi_tax_data(force_update=True)
    mock_dl.assert_called_once_with(str(tmp_path))


def test_get_ncbi_tax_data_missing_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("TAXONKIT_DB", str(tmp_path))
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_ncbi_tax_data") as mock_dl:
        get_ncbi_tax_data()
    mock_dl.assert_called_once_with(str(tmp_path))
