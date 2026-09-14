"""
Tests for the NCBI-specific half of the assembly-info setup.

The download, verification and cleanup all live in
bit/modules/hosted_parquet_asset.py now and are tested there, against both this asset
and the GTDB one. What's left here is what's genuinely NCBI's: that the asset is
configured with the right variable, filenames and URLs, that the wrappers the rest of
the codebase imports still delegate to it, the present/absent/force routing (and
bit's `quiet` flag, which lives in the wrapper), and read_date_retrieved(), which has
no GTDB counterpart.

Mirrors gtotree/tests/utils/ncbi/test_get_ncbi_assembly_data.py.
"""

import pytest  # type: ignore
from unittest.mock import patch

from bit.modules.ncbi.get_ncbi_assembly_data import (
    DATE_FILENAME,
    NCBI_ASSET,
    NCBI_DATA_URL,
    NCBI_DATE_URL,
    PARQUET_FILENAME,
    check_if_data_present,
    check_ncbi_assembly_info_location_var_is_set,
    get_ncbi_assembly_data,
    ncbi_data_table_path,
    read_date_retrieved,
)

MODPATH = "bit.modules.ncbi.get_ncbi_assembly_data"
ENV_VAR = "NCBI_assembly_data_dir"

# the asset's date stamp shape: a single 'YYYY,MM,DD' line
_DATE_BODY = "2026,01,05\n"


# --- the asset spec --------------------------------------------------------

class TestAssetSpec:

    def test_it_points_at_the_ncbi_variable_and_files(self):
        assert NCBI_ASSET.env_var == ENV_VAR
        assert NCBI_ASSET.parquet_filename == PARQUET_FILENAME == "ncbi-data.parquet"
        assert NCBI_ASSET.sidecar_filename == DATE_FILENAME == "date-retrieved.txt"

    def test_the_module_urls_match_the_asset(self):
        # these two names are imported by tests elsewhere, so they have to stay in step
        assert NCBI_DATA_URL == NCBI_ASSET.data_url
        assert NCBI_DATE_URL == NCBI_ASSET.sidecar_url

    def test_the_asset_is_pulled_from_the_rolling_release(self):
        assert NCBI_DATA_URL.startswith("https://github.com/AstrobioMike/bit/releases")
        assert NCBI_DATA_URL.endswith(PARQUET_FILENAME)


# --- the wrappers delegate -------------------------------------------------

class TestWrappersDelegateToTheAsset:

    def test_location_var_returns_the_path(self, monkeypatch, tmp_path):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        assert check_ncbi_assembly_info_location_var_is_set() == str(tmp_path)

    def test_location_var_exits_nonzero_if_missing(self, monkeypatch):
        # this used to exit 0, so a wrapper script read a broken install as success
        monkeypatch.delenv(ENV_VAR, raising=False)
        with pytest.raises(SystemExit) as excinfo:
            check_ncbi_assembly_info_location_var_is_set()
        assert excinfo.value.code == 1

    def test_table_path_derives_from_the_filename_constant(self, monkeypatch,
                                                           tmp_path):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        assert ncbi_data_table_path() == str(tmp_path / PARQUET_FILENAME)
        assert ncbi_data_table_path("/somewhere") == f"/somewhere/{PARQUET_FILENAME}"

    def test_presence_check_delegates(self, tmp_path):
        (tmp_path / PARQUET_FILENAME).write_text("x")
        (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
        assert check_if_data_present(str(tmp_path)) is True


# --- routing ---------------------------------------------------------------

def _seed(tmp_path):
    (tmp_path / PARQUET_FILENAME).write_text("x")
    (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)


class TestRouting:
    """
    get_ncbi_assembly_data() returns the table path -- get_accessions_from_gtdb caches
    it as the liveness-screen path, so the return value is load-bearing.
    """

    def test_a_present_asset_is_not_re_downloaded_and_the_table_path_comes_back(
            self, monkeypatch, tmp_path):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        _seed(tmp_path)

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            result = get_ncbi_assembly_data(force_update=False, quiet=True)

        assert result == str(tmp_path / PARQUET_FILENAME)
        mock_dl.assert_not_called()

    def test_an_absent_asset_is_downloaded_and_the_table_path_still_comes_back(
            self, monkeypatch, tmp_path):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            result = get_ncbi_assembly_data(quiet=True)

        assert result == str(tmp_path / PARQUET_FILENAME)
        mock_dl.assert_called_once()

    def test_force_update_downloads_even_if_present(self, monkeypatch, tmp_path):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        _seed(tmp_path)

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            get_ncbi_assembly_data(force_update=True, quiet=True)

        mock_dl.assert_called_once()


# --- the quiet flag --------------------------------------------------------

class TestQuiet:
    """
    `quiet` silences the "already present" note and nothing else. Everything outside
    `bit data get ncbi-assembly-data` passes quiet=True, so this note is what would
    otherwise land in the middle of a gen-mg, dl-ncbi-assemblies or get-accs run.
    """

    def test_the_present_note_is_printed_by_default(self, monkeypatch, tmp_path,
                                                    capsys):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        _seed(tmp_path)

        get_ncbi_assembly_data()

        out = capsys.readouterr().out
        assert "Assembly data already present at:" in out
        assert "bit data get ncbi-assembly-data -f" in out

    def test_quiet_silences_the_present_note(self, monkeypatch, tmp_path, capsys):
        monkeypatch.setenv(ENV_VAR, str(tmp_path))
        _seed(tmp_path)

        get_ncbi_assembly_data(quiet=True)

        assert capsys.readouterr().out == ""

    def test_quiet_does_not_reach_the_download(self, monkeypatch, tmp_path):
        """
        The download path takes no `quiet`: a failure there is fatal, so it always
        explains itself regardless of what the caller asked for.
        """
        monkeypatch.setenv(ENV_VAR, str(tmp_path))

        with patch(f"{MODPATH}.get_slim_ncbi_assembly_data") as mock_dl:
            get_ncbi_assembly_data(quiet=True)

        mock_dl.assert_called_once_with(str(tmp_path))


# --- read_date_retrieved ---------------------------------------------------

class TestReadDateRetrieved:

    def test_a_well_formed_stamp_is_formatted(self, tmp_path):
        (tmp_path / DATE_FILENAME).write_text(_DATE_BODY)
        assert read_date_retrieved(str(tmp_path)) == "Jan 05, 2026"

    def test_an_unparseable_stamp_comes_back_raw(self, tmp_path):
        (tmp_path / DATE_FILENAME).write_text("not-a-date\n")
        assert read_date_retrieved(str(tmp_path)) == "not-a-date"
