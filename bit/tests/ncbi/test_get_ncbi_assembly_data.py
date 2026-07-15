import socket

import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest  # type: ignore
from unittest.mock import patch
import bit.modules.ncbi.get_ncbi_assembly_data as mod
from bit.modules.ncbi.get_ncbi_assembly_data import (
    NCBI_DATA_URL,
    NCBI_DATE_URL,
    check_ncbi_assembly_info_location_var_is_set,
    ncbi_data_table_path,
    check_if_data_present,
    get_slim_ncbi_assembly_data,
    get_ncbi_assembly_data,
)
from bit.modules.ncbi.build_ncbi_data_parquet import PARQUET_FILENAME, DATE_FILENAME

# --- helpers --------------------------------------------------------------

def _valid_parquet_bytes(path):
    pq.write_table(pa.table({"assembly_accession": pa.array(["GCF_1"])}), str(path))


def _fake_downloader(date_contents="2026,01,05"):
    """
    Stand-in for download_with_tqdm that serves BOTH assets by URL:
      - the parquet URL  -> a valid one-row parquet
      - the date URL     -> `date_contents`
    """
    def _dl(url, target, filename=None, **kw):
        if url == NCBI_DATA_URL:
            _valid_parquet_bytes(filename)
        elif url == NCBI_DATE_URL:
            with open(filename, "w") as fh:
                fh.write(date_contents + ("\n" if not date_contents.endswith("\n") else ""))
        else:
            raise AssertionError(f"unexpected download URL: {url}")
    return _dl


# --- location var ---------------------------------------------------------

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    assert check_ncbi_assembly_info_location_var_is_set() == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("NCBI_assembly_data_dir", raising=False)
    with pytest.raises(SystemExit):
        check_ncbi_assembly_info_location_var_is_set()


def test_table_path_derives_from_the_single_filename_constant(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    assert ncbi_data_table_path() == str(tmp_path / PARQUET_FILENAME)
    assert ncbi_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"


# --- check_if_data_present ------------------------------------------------

def test_present_when_both_files_nonempty(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text("2026,01,01")
    assert check_if_data_present(str(tmp_path)) is True


def test_absent_when_table_missing(tmp_path):
    (tmp_path / DATE_FILENAME).write_text("2026,01,01")
    assert check_if_data_present(str(tmp_path)) is False


def test_absent_when_date_missing(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_data_present(str(tmp_path)) is False


def test_a_half_present_pair_is_cleaned_up(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()


def test_empty_files_count_as_absent_and_are_removed(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("")
    (tmp_path / DATE_FILENAME).write_text("")
    assert check_if_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / DATE_FILENAME).exists()


# --- download path --------------------------------------------------------

def test_download_writes_table_and_date(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               _fake_downloader()):
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert (tmp_path / PARQUET_FILENAME).exists()
    assert (tmp_path / DATE_FILENAME).exists()


def test_date_file_comes_FROM_THE_ASSET_not_todays_date(tmp_path):
    """
    The recorded date must be the asset's build-time stamp (when the SNAPSHOT was
    taken), NOT the day the user downloaded it. A user pulling a 3-week-old asset must
    see the 3-week-old date, or 'date-retrieved' silently means the wrong thing.
    """
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               _fake_downloader(date_contents="2025,12,15")):
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert (tmp_path / DATE_FILENAME).read_text().strip() == "2025,12,15"


def test_no_leftover_part_file_after_success(tmp_path):
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               _fake_downloader()):
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert not (tmp_path / (DATE_FILENAME + ".part")).exists()


def test_a_malformed_date_stamp_is_rejected(tmp_path):
    """A date file that isn't 'YYYY,MM,DD' must fail rather than be trusted."""
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               _fake_downloader(date_contents="garbage")), \
         patch("bit.modules.ncbi.get_ncbi_assembly_data.notify_premature_exit") as mock_exit:
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    # table and date both cleaned up, and no stray .part left behind
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / DATE_FILENAME).exists()
    assert not (tmp_path / (DATE_FILENAME + ".part")).exists()
    mock_exit.assert_called_once()


def test_a_truncated_parquet_is_rejected_and_cleaned_up(tmp_path):
    def _bad(url, target, filename=None, **kw):
        if url == NCBI_DATA_URL:
            with open(filename, "wb") as fh:
                fh.write(b"not a parquet file")
        else:
            with open(filename, "w") as fh:
                fh.write("2026,01,01\n")

    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm", _bad), \
         patch("bit.modules.ncbi.get_ncbi_assembly_data.notify_premature_exit") as mock_exit:
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / DATE_FILENAME).exists()
    mock_exit.assert_called_once()


def test_download_failure_exits_without_a_local_rebuild(tmp_path):
    def boom(*a, **k):
        raise socket.timeout("slow")
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm", boom), \
         patch("bit.modules.ncbi.get_ncbi_assembly_data.notify_premature_exit") as mock_exit:
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert not (tmp_path / PARQUET_FILENAME).exists()
    mock_exit.assert_called_once()


def test_socket_timeout_is_restored_after_download(tmp_path):
    before = socket.getdefaulttimeout()
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.download_with_tqdm",
               _fake_downloader()):
        get_slim_ncbi_assembly_data(str(tmp_path), quiet=True)
    assert socket.getdefaulttimeout() == before


# --- routing + the return-path fix ----------------------------------------

def test_present_data_skips_download_but_still_returns_path(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text("2026,01,01")
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.get_slim_ncbi_assembly_data") as mock_dl:
        result = get_ncbi_assembly_data(quiet=True)
    mock_dl.assert_not_called()
    assert result == str(tmp_path / PARQUET_FILENAME)


def test_missing_data_triggers_download_then_returns_path(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.get_slim_ncbi_assembly_data") as mock_dl:
        result = get_ncbi_assembly_data(quiet=True)
    mock_dl.assert_called_once()
    assert result == str(tmp_path / PARQUET_FILENAME)


def test_force_update_downloads_even_when_present(tmp_path, monkeypatch):
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text("2026,01,01")
    with patch("bit.modules.ncbi.get_ncbi_assembly_data.get_slim_ncbi_assembly_data") as mock_dl:
        result = get_ncbi_assembly_data(force_update=True, quiet=True)
    mock_dl.assert_called_once()
    assert result == str(tmp_path / PARQUET_FILENAME)
