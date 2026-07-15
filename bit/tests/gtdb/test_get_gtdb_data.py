import socket
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
import pytest  # type: ignore
from unittest.mock import patch
from bit.modules.gtdb.get_gtdb_data import (
    PARQUET_FILENAME,
    VERSION_FILENAME,
    GTDB_DATA_URL,
    GTDB_VERSION_URL,
    check_gtdb_location_var_is_set,
    gtdb_data_table_path,
    check_if_gtdb_data_present,
    get_slim_gtdb_tab,
    get_gtdb_data,
    report_gtdb_version_info,
)


# --- helpers --------------------------------------------------------------

def _valid_parquet(path):
    pq.write_table(pa.table({"accession": pa.array(["GB_GCA_1"])}), str(path))


# GTDB's real VERSION.txt shape: version line, blank, "Released ..." line
_VERSION_BODY = "v232\n\nReleased Apr 15, 2026\n"


def _fake_downloader(version_body=_VERSION_BODY):
    """Serves both assets by URL: the parquet URL -> a valid parquet, the version
    URL -> `version_body`."""
    def _dl(url, target, filename=None, **kw):
        if url == GTDB_DATA_URL:
            _valid_parquet(filename)
        elif url == GTDB_VERSION_URL:
            with open(filename, "w") as fh:
                fh.write(version_body)
        else:
            raise AssertionError(f"unexpected download URL: {url}")
    return _dl


# --- location var ---------------------------------------------------------

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    assert check_gtdb_location_var_is_set() == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("GTDB_DIR", raising=False)
    with pytest.raises(SystemExit):
        check_gtdb_location_var_is_set()


def test_table_path_derives_from_the_single_filename_constant(monkeypatch, tmp_path):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    assert gtdb_data_table_path() == str(tmp_path / PARQUET_FILENAME)
    assert gtdb_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"


# --- check_if_gtdb_data_present -------------------------------------------

def test_present_when_both_files_nonempty(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / VERSION_FILENAME).write_text("v232\nReleased\n")
    assert check_if_gtdb_data_present(str(tmp_path)) is True


def test_absent_when_table_missing(tmp_path):
    (tmp_path / VERSION_FILENAME).write_text("v232\nReleased\n")
    assert check_if_gtdb_data_present(str(tmp_path)) is False


def test_absent_when_version_missing(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_gtdb_data_present(str(tmp_path)) is False


def test_a_half_present_pair_is_cleaned_up(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()


def test_empty_files_count_as_absent_and_are_removed(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("")
    (tmp_path / VERSION_FILENAME).write_text("")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / VERSION_FILENAME).exists()


# --- download path --------------------------------------------------------

def test_download_writes_table_and_version(tmp_path):
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", _fake_downloader()):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert (tmp_path / PARQUET_FILENAME).exists()
    assert (tmp_path / VERSION_FILENAME).exists()


def test_version_file_comes_FROM_THE_ASSET(tmp_path):
    """The recorded version must be the asset's actual GTDB release, fetched not
    fabricated. A downloaded v999 body must land as v999."""
    body = "v999\n\nReleased Jan 1, 2099\n"
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", _fake_downloader(body)):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    version, date = report_gtdb_version_info(str(tmp_path))
    assert version == "v999"
    assert date == "Released Jan 1, 2099"


def test_version_asset_lands_under_its_asset_name(tmp_path):
    """The version file keeps its asset name (VERSION.txt) on disk -- no rename --
    and no stray .part is left behind."""
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", _fake_downloader()):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert (tmp_path / "VERSION.txt").exists()
    assert not (tmp_path / (VERSION_FILENAME + ".part")).exists()


def test_a_truncated_parquet_is_rejected_and_cleaned_up(tmp_path):
    def _bad(url, target, filename=None, **kw):
        if url == GTDB_DATA_URL:
            with open(filename, "wb") as fh:
                fh.write(b"not a parquet")
        else:
            with open(filename, "w") as fh:
                fh.write(_VERSION_BODY)

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", _bad), \
         patch("bit.modules.gtdb.get_gtdb_data.notify_premature_exit"):
        with pytest.raises(SystemExit):
            get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / VERSION_FILENAME).exists()


def test_a_malformed_version_file_is_rejected(tmp_path):
    """A version file without the two expected lines must fail, not be trusted."""
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               _fake_downloader(version_body="oneline\n")), \
         patch("bit.modules.gtdb.get_gtdb_data.notify_premature_exit"):
        with pytest.raises(SystemExit):
            get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert not (tmp_path / PARQUET_FILENAME).exists()
    assert not (tmp_path / VERSION_FILENAME).exists()
    assert not (tmp_path / (VERSION_FILENAME + ".part")).exists()


def test_download_failure_exits_without_a_local_rebuild(tmp_path):
    def boom(*a, **k):
        raise socket.timeout("slow")
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", boom), \
         patch("bit.modules.gtdb.get_gtdb_data.notify_premature_exit"):
        with pytest.raises(SystemExit):
            get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert not (tmp_path / PARQUET_FILENAME).exists()


def test_socket_timeout_is_restored_after_download(tmp_path):
    before = socket.getdefaulttimeout()
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", _fake_downloader()):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    assert socket.getdefaulttimeout() == before


# --- get_gtdb_data (routing + the return-path contract) -------------------

def test_present_data_skips_download_but_still_returns_path(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / VERSION_FILENAME).write_text("v232\nReleased\n")
    with patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_dl:
        result = get_gtdb_data(quiet=True)
    mock_dl.assert_not_called()
    assert result == str(tmp_path / PARQUET_FILENAME)


def test_missing_data_triggers_download_then_returns_path(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    with patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_dl:
        result = get_gtdb_data(quiet=True)
    mock_dl.assert_called_once()
    assert result == str(tmp_path / PARQUET_FILENAME)


def test_force_update_downloads_even_when_present(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / VERSION_FILENAME).write_text("v232\nReleased\n")
    with patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_dl:
        result = get_gtdb_data(force_update=True, quiet=True)
    mock_dl.assert_called_once()
    assert result == str(tmp_path / PARQUET_FILENAME)
