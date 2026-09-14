"""
Unit tests for bit/modules/hosted_parquet_asset.py.

The NCBI assembly-info table and the GTDB metadata table run through one
implementation now, so the behaviour that used to be tested twice (once per source,
with the two copies slowly drifting) is parametrized over both assets here. That's
the point of the file: whatever the shared code does, it has to do for both.

What stays in the per-source test modules is the source-specific part -- that the
asset is configured with the right env var, filenames and URLs, the present/absent/
force routing (including bit's `quiet` flag, which lives in the wrappers), and the
sidecar readers (read_date_retrieved, report_gtdb_version_info).

Mirrors gtotree/tests/utils/misc/test_hosted_parquet_asset.py.
"""

import socket
import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore
from unittest.mock import patch

from bit.modules import hosted_parquet_asset as hpa
from bit.modules.gtdb.get_gtdb_data import GTDB_ASSET
from bit.modules.ncbi.get_ncbi_assembly_data import NCBI_ASSET

MODPATH = "bit.modules.hosted_parquet_asset"

#: a well-formed sidecar body for each asset, matching its validator.
#: GTDB's real VERSION.txt shape is: version line, blank, "Released ..." line
SIDECAR_BODIES = {
    "NCBI_assembly_data_dir": "2026,01,05\n",
    "GTDB_DIR": "v232\n\nReleased Apr 15, 2026\n",
}

#: a sidecar body each asset's validator must reject
BAD_SIDECAR_BODIES = {
    "NCBI_assembly_data_dir": "not-a-date\n",
    "GTDB_DIR": "v232\n",  # version line but no release-date line
}


@pytest.fixture(params=[NCBI_ASSET, GTDB_ASSET], ids=["ncbi", "gtdb"])
def asset(request):
    return request.param


def _sidecar_body(asset):
    return SIDECAR_BODIES[asset.env_var]


def _valid_parquet(path):
    pq.write_table(pa.table({"assembly_accession": pa.array(["GCA_1"])}), str(path))


def _fake_downloader(asset, sidecar_body=None, break_parquet=False):
    """
    Stand-in for download_with_tqdm(url, target, filename, ...): serves the parquet URL
    to a real (small) Parquet file and the sidecar URL to `sidecar_body`.
    """
    body = _sidecar_body(asset) if sidecar_body is None else sidecar_body

    def _dl(url, target, filename=None, **kw):
        if url == asset.data_url:
            if break_parquet:
                with open(filename, "w") as fh:
                    fh.write("not a parquet")
            else:
                _valid_parquet(filename)
        elif url == asset.sidecar_url:
            with open(filename, "w") as fh:
                fh.write(body)
        else:
            raise AssertionError(f"unexpected download URL: {url}")
    return _dl


def _seed(location, asset, table="x", sidecar=None):
    """Put a present-looking pair of files in `location`."""
    if table is not None:
        (location / asset.parquet_filename).write_text(table)
    if sidecar is not False:
        body = _sidecar_body(asset) if sidecar is None else sidecar
        (location / asset.sidecar_filename).write_text(body)


# ---------------------------------------------------------------------------
# the sidecar validators
# ---------------------------------------------------------------------------

class TestSidecarValidators:

    @pytest.mark.parametrize("body", ["2026,01,05\n", "2026,1,5\ntrailing junk\n"])
    def test_a_well_formed_date_stamp_passes(self, tmp_path, body):
        path = tmp_path / "date-retrieved.txt"
        path.write_text(body)
        hpa.validate_date_stamp(str(path))  # no raise

    @pytest.mark.parametrize("body", ["not-a-date\n", "2026,01\n", "2026,01,05,06\n",
                                      "2026,Jan,05\n", "\n"])
    def test_a_malformed_date_stamp_is_rejected(self, tmp_path, body):
        path = tmp_path / "date-retrieved.txt"
        path.write_text(body)
        with pytest.raises(ValueError):
            hpa.validate_date_stamp(str(path))

    def test_two_version_lines_pass(self, tmp_path):
        path = tmp_path / "VERSION.txt"
        path.write_text("v232\nReleased Apr 15, 2026\n")
        hpa.validate_version_lines(str(path))  # no raise

    def test_blank_lines_do_not_count_toward_the_two(self, tmp_path):
        path = tmp_path / "VERSION.txt"
        path.write_text("v232\n\n\n")
        with pytest.raises(ValueError):
            hpa.validate_version_lines(str(path))

    @pytest.mark.parametrize("body", ["", "v232\n"])
    def test_too_few_version_lines_is_rejected(self, tmp_path, body):
        path = tmp_path / "VERSION.txt"
        path.write_text(body)
        with pytest.raises(ValueError):
            hpa.validate_version_lines(str(path))


# ---------------------------------------------------------------------------
# locating the asset
# ---------------------------------------------------------------------------

class TestLocation:

    def test_the_env_var_gives_the_directory(self, asset, tmp_path, monkeypatch):
        monkeypatch.setenv(asset.env_var, str(tmp_path))
        assert asset.location() == str(tmp_path)

    def test_an_unset_env_var_exits_nonzero(self, asset, monkeypatch, capsys):
        """
        This diverged once already: NCBI exited 0 here while GTDB exited 1, so a
        wrapper script couldn't tell a broken install from a successful one.
        """
        monkeypatch.delenv(asset.env_var, raising=False)

        with pytest.raises(SystemExit) as excinfo:
            asset.location()

        assert excinfo.value.code == 1
        out = capsys.readouterr().out
        assert asset.env_var in out
        assert "bit data locations check" in out

    def test_table_path_uses_the_env_var_when_not_given_one(self, asset, tmp_path,
                                                            monkeypatch):
        monkeypatch.setenv(asset.env_var, str(tmp_path))
        assert asset.table_path() == str(tmp_path / asset.parquet_filename)

    def test_table_path_honours_an_explicit_location(self, asset):
        assert asset.table_path("/somewhere") == \
            f"/somewhere/{asset.parquet_filename}"

    def test_sidecar_path_sits_beside_the_table(self, asset):
        assert asset.sidecar_path("/somewhere") == \
            f"/somewhere/{asset.sidecar_filename}"


# ---------------------------------------------------------------------------
# is_present -- and its cleanup side effect
# ---------------------------------------------------------------------------

class TestIsPresent:

    def test_both_files_nonempty_is_present(self, asset, tmp_path):
        _seed(tmp_path, asset)
        assert asset.is_present(str(tmp_path)) is True

    def test_a_missing_table_is_absent(self, asset, tmp_path):
        _seed(tmp_path, asset, table=None)
        assert asset.is_present(str(tmp_path)) is False

    def test_a_missing_sidecar_is_absent(self, asset, tmp_path):
        _seed(tmp_path, asset, sidecar=False)
        assert asset.is_present(str(tmp_path)) is False

    def test_a_half_present_pair_is_cleaned_up(self, asset, tmp_path):
        """
        The table alone would be read and trusted on the next run, so the leftover has
        to go rather than just being reported absent.
        """
        _seed(tmp_path, asset, sidecar=False)

        assert asset.is_present(str(tmp_path)) is False
        assert not (tmp_path / asset.parquet_filename).exists()

    def test_empty_files_count_as_absent_and_are_removed(self, asset, tmp_path):
        _seed(tmp_path, asset, table="", sidecar="")

        assert asset.is_present(str(tmp_path)) is False
        assert not (tmp_path / asset.parquet_filename).exists()
        assert not (tmp_path / asset.sidecar_filename).exists()


# ---------------------------------------------------------------------------
# the download path
# ---------------------------------------------------------------------------

class TestDownload:

    def test_both_files_land_and_the_table_is_readable(self, asset, tmp_path):
        with patch(f"{MODPATH}.download_with_tqdm",
                   side_effect=_fake_downloader(asset)):
            asset.download(str(tmp_path))

        table = tmp_path / asset.parquet_filename
        assert pq.ParquetFile(str(table)).metadata.num_rows == 1
        assert (tmp_path / asset.sidecar_filename).read_text() == _sidecar_body(asset)

    def test_a_corrupt_table_is_caught_and_both_files_removed(self, asset, tmp_path):
        with patch(f"{MODPATH}.download_with_tqdm",
                   side_effect=_fake_downloader(asset, break_parquet=True)), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert not (tmp_path / asset.parquet_filename).exists()
        assert not (tmp_path / asset.sidecar_filename).exists()

    def test_a_network_error_is_translated_rather_than_raised(self, asset, tmp_path):
        def _boom(*a, **k):
            raise TimeoutError("no route")

        with patch(f"{MODPATH}.download_with_tqdm", side_effect=_boom), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert not (tmp_path / asset.parquet_filename).exists()

    def test_a_malformed_sidecar_is_rejected_and_cleaned_up(self, asset, tmp_path):
        bad = BAD_SIDECAR_BODIES[asset.env_var]

        with patch(f"{MODPATH}.download_with_tqdm",
                   side_effect=_fake_downloader(asset, sidecar_body=bad)), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert not (tmp_path / asset.sidecar_filename).exists()

    def test_no_part_file_is_left_behind_on_failure(self, asset, tmp_path):
        bad = BAD_SIDECAR_BODIES[asset.env_var]

        with patch(f"{MODPATH}.download_with_tqdm",
                   side_effect=_fake_downloader(asset, sidecar_body=bad)), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert list(tmp_path.glob("*.part")) == []

    def test_the_socket_timeout_is_restored_afterwards(self, asset, tmp_path):
        before = socket.getdefaulttimeout()
        with patch(f"{MODPATH}.download_with_tqdm",
                   side_effect=_fake_downloader(asset)):
            asset.download(str(tmp_path))

        assert socket.getdefaulttimeout() == before

    def test_the_socket_timeout_is_restored_even_on_failure(self, asset, tmp_path):
        def _boom(*a, **k):
            raise TimeoutError("no route")

        before = socket.getdefaulttimeout()
        with patch(f"{MODPATH}.download_with_tqdm", side_effect=_boom), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert socket.getdefaulttimeout() == before


class TestReportUnavailable:

    def test_the_manual_route_names_both_urls_and_both_filenames(self, asset, capsys):
        with patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.report_unavailable(TimeoutError("no route"))

        out = capsys.readouterr().out
        assert asset.data_url in out
        assert asset.sidecar_url in out
        assert asset.parquet_filename in out
        assert asset.sidecar_filename in out
        assert "no route" in out

    def test_the_explanation_is_printed_even_though_callers_pass_quiet(self, asset,
                                                                       tmp_path,
                                                                       capsys):
        """
        Every caller outside `bit data get ...` passes quiet=True, and the failure is
        always fatal -- so if `quiet` reached this path the user would get nothing but
        "Exiting for now :(" to act on. `quiet` deliberately stops at the wrappers.
        """
        def _boom(*a, **k):
            raise TimeoutError("no route")

        with patch(f"{MODPATH}.download_with_tqdm", side_effect=_boom), \
             patch(f"{MODPATH}.notify_premature_exit", side_effect=SystemExit(1)):
            with pytest.raises(SystemExit):
                asset.download(str(tmp_path))

        assert "Couldn't download" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# the two assets are configured distinctly
# ---------------------------------------------------------------------------

class TestAssetConfiguration:
    """
    One implementation behind two assets means a copy/paste slip in either wrapper
    would point both at the same files. These pin the pairs apart.
    """

    def test_the_assets_do_not_share_an_env_var(self):
        assert NCBI_ASSET.env_var != GTDB_ASSET.env_var

    def test_the_assets_do_not_share_a_table(self):
        assert NCBI_ASSET.parquet_filename != GTDB_ASSET.parquet_filename
        assert NCBI_ASSET.data_url != GTDB_ASSET.data_url

    def test_the_assets_do_not_share_a_sidecar(self):
        assert NCBI_ASSET.sidecar_filename != GTDB_ASSET.sidecar_filename
        assert NCBI_ASSET.sidecar_url != GTDB_ASSET.sidecar_url

    def test_each_asset_gets_the_validator_for_its_sidecar_shape(self):
        assert NCBI_ASSET.sidecar_validator is hpa.validate_date_stamp
        assert GTDB_ASSET.sidecar_validator is hpa.validate_version_lines

    def test_the_urls_are_built_from_the_release_base_and_filenames(self, asset):
        assert asset.data_url == f"{asset.release_base}/{asset.parquet_filename}"
        assert asset.sidecar_url == f"{asset.release_base}/{asset.sidecar_filename}"
