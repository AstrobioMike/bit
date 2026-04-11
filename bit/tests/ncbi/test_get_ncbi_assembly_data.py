import pytest # type: ignore
from datetime import date
from pathlib import Path
from unittest.mock import patch, MagicMock

from bit.modules.ncbi.get_ncbi_assembly_data import (
    check_ncbi_assembly_info_location_var_is_set,
    check_if_data_present,
    download_ncbi_assembly_summary_data,
    get_ncbi_assembly_data,
)


def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    result = check_ncbi_assembly_info_location_var_is_set()
    assert result == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("NCBI_assembly_data_dir", raising=False)
    with pytest.raises(SystemExit):
        check_ncbi_assembly_info_location_var_is_set()


def test_check_if_data_present_both_files_nonempty(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    assert check_if_data_present(str(tmp_path)) is True


def test_check_if_data_present_table_missing(tmp_path):
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "date-retrieved.txt").exists()


def test_check_if_data_present_date_file_missing(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "ncbi-assembly-info.tsv").exists()


def test_check_if_data_present_both_missing(tmp_path):
    assert check_if_data_present(str(tmp_path)) is False


def test_check_if_data_present_empty_files_removed(tmp_path):
    (tmp_path / "ncbi-assembly-info.tsv").write_text("")
    (tmp_path / "date-retrieved.txt").write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / "ncbi-assembly-info.tsv").exists()
    assert not (tmp_path / "date-retrieved.txt").exists()


def _mock_download(link, label, dest):
    """Side effect for download_with_tqdm that writes fake content to dest."""
    if "genbank" in link:
        Path(dest).write_text("genbank_line1\ngenbank_line2\n")
    else:
        Path(dest).write_text("refseq_line1\nrefseq_line2\n")


def test_download_combines_genbank_and_refseq(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               side_effect=_mock_download):
        download_ncbi_assembly_summary_data(str(tmp_path))
    combined = (tmp_path / "ncbi-assembly-info.tsv").read_text()
    assert "genbank_line1" in combined
    assert "refseq_line1" in combined


def test_download_removes_temp_file(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               side_effect=_mock_download):
        download_ncbi_assembly_summary_data(str(tmp_path))
    assert not (tmp_path / "refseq-assembly-info.tmp").exists()


def test_download_writes_date_file(tmp_path):
    fixed_date = date(2024, 6, 15)
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               side_effect=_mock_download), \
         patch("bit.modules.ncbi.get_ncbi_assembly_data.date") as mock_date:
        mock_date.today.return_value = fixed_date
        download_ncbi_assembly_summary_data(str(tmp_path))
    date_text = (tmp_path / "date-retrieved.txt").read_text().strip()
    assert date_text == "2024,06,15"


def test_download_failure_calls_notify_premature_exit(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               side_effect=Exception("network error")), \
         patch("bit.modules.ncbi.get_ncbi_assembly_data.notify_premature_exit") as mock_exit:
        mock_exit.side_effect = SystemExit(1)
        with pytest.raises(SystemExit):
            download_ncbi_assembly_summary_data(str(tmp_path))
    mock_exit.assert_called_once()


def test_get_data_already_present_skips_download(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_ncbi_assembly_summary_data") as mock_dl:
        get_ncbi_assembly_data(force_update=False, quiet=True)
    mock_dl.assert_not_called()


def test_get_data_force_update_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / "ncbi-assembly-info.tsv").write_text("data")
    (tmp_path / "date-retrieved.txt").write_text("2024,01,01")
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_ncbi_assembly_summary_data") as mock_dl:
        get_ncbi_assembly_data(force_update=True)
    mock_dl.assert_called_once_with(str(tmp_path))


def test_get_data_missing_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_ncbi_assembly_summary_data") as mock_dl:
        get_ncbi_assembly_data()
    mock_dl.assert_called_once_with(str(tmp_path))
