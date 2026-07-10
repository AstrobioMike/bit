import gzip
import io
import os
import shutil
import socket
import tarfile
import urllib.error
import pytest
from pathlib import Path
from unittest.mock import patch

import bit.modules.gtdb.get_gtdb_data as gtdb_mod
from bit.modules.gtdb.get_gtdb_data import (
    GTDB_KEPT_COLUMNS,
    check_gtdb_location_var_is_set,
    check_if_gtdb_data_present,
    gen_gtdb_tab,
    get_gtdb_data,
    get_slim_gtdb_tab,
)


ARC_TAXONOMY = "d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;f__Haloarculaceae;g__Haloarcula;s__Haloarcula hispanica"
BAC_TAXONOMY = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"

ARC_TSV = f"accession\tgtdb_taxonomy\nRS_GCF_000002125.1\t{ARC_TAXONOMY}\n"
BAC_TSV = f"accession\tgtdb_taxonomy\nGB_GCA_000001405.29\t{BAC_TAXONOMY}\n"


def _to_gz(content: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content.encode())
    return buf.getvalue()


def _make_mock_download(arc_tsv=ARC_TSV, bac_tsv=BAC_TSV):
    """Returns a side_effect for download_with_tqdm that writes fake gzipped TSVs."""
    def _download(link, label, dest, **kwargs):
        if "ar53" in link:
            Path(dest).write_bytes(_to_gz(arc_tsv))
        else:
            Path(dest).write_bytes(_to_gz(bac_tsv))
        return dest
    return _download


def _build_slim_tarball(dest, *, include_metadata=True, include_version=True,
                        extra_member=None, metadata_under_dir=False):
    """
    Build a .tar.gz mimicking the hosted slim bundle. By default it contains the
    two expected files at the archive root. Flags let tests omit a file, add a
    stray member that must be ignored, or nest the metadata under a directory.
    """
    meta = "accession\tdomain\tgenus\nRS_GCF_000002125.1\tArchaea\tHaloarcula\n"
    ver = "v220\n2024-04-24\n"
    with tarfile.open(dest, "w:gz") as tar:
        def _add(name, content):
            data = content.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        if include_metadata:
            name = "some-dir/GTDB-arc-and-bac-metadata.tsv" if metadata_under_dir \
                else "GTDB-arc-and-bac-metadata.tsv"
            _add(name, meta)
        if include_version:
            _add("GTDB-version-info.txt", ver)
        if extra_member:
            _add(extra_member, "ignore me\n")


def _mock_tarball_download(src_tarball):
    """download_with_tqdm side_effect that copies a prebuilt tarball to dest."""
    def _download(url, label, dest, **kwargs):
        shutil.copy(src_tarball, dest)
        return dest
    return _download


################################################################################
# check_gtdb_location_var_is_set
################################################################################

def test_location_var_returns_path(monkeypatch, tmp_path):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    result = check_gtdb_location_var_is_set()
    assert result == str(tmp_path)


def test_location_var_exits_if_missing(monkeypatch):
    monkeypatch.delenv("GTDB_DIR", raising=False)
    with pytest.raises(SystemExit):
        check_gtdb_location_var_is_set()


################################################################################
# check_if_gtdb_data_present
################################################################################

def test_check_if_data_present_both_files_nonempty(tmp_path):
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    assert check_if_gtdb_data_present(str(tmp_path)) is True


def test_check_if_data_present_metadata_missing(tmp_path):
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / "GTDB-version-info.txt").exists()


def test_check_if_data_present_version_missing(tmp_path):
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()


def test_check_if_data_present_both_missing(tmp_path):
    assert check_if_gtdb_data_present(str(tmp_path)) is False


def test_check_if_data_present_empty_files_removed(tmp_path):
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("")
    (tmp_path / "GTDB-version-info.txt").write_text("")
    assert check_if_gtdb_data_present(str(tmp_path)) is False
    assert not (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert not (tmp_path / "GTDB-version-info.txt").exists()


################################################################################
# gen_gtdb_tab
################################################################################

def test_gen_gtdb_tab_writes_metadata_file(tmp_path):
    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download()), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        gen_gtdb_tab(str(tmp_path))
    assert (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert (tmp_path / "GTDB-arc-and-bac-metadata.tsv").stat().st_size > 0


def test_gen_gtdb_tab_writes_version_file(tmp_path):
    def fake_urlretrieve(url, dest):
        Path(dest).write_text("v220\n")

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download()), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve",
               side_effect=fake_urlretrieve):
        gen_gtdb_tab(str(tmp_path))
    assert (tmp_path / "GTDB-version-info.txt").read_text() == "v220\n"


def test_gen_gtdb_tab_taxonomy_columns_split(tmp_path):
    import pandas as pd

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download()), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        gen_gtdb_tab(str(tmp_path))

    result = pd.read_csv(tmp_path / "GTDB-arc-and-bac-metadata.tsv", sep="\t")
    for col in ("domain", "phylum", "class", "order", "family", "genus", "species"):
        assert col in result.columns


def test_gen_gtdb_tab_taxonomy_values_correct(tmp_path):
    import pandas as pd

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download()), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        gen_gtdb_tab(str(tmp_path))

    result = pd.read_csv(tmp_path / "GTDB-arc-and-bac-metadata.tsv", sep="\t")
    arc_row = result[result["accession"] == "RS_GCF_000002125.1"].iloc[0]
    assert arc_row["domain"] == "Archaea"
    assert arc_row["genus"] == "Haloarcula"
    bac_row = result[result["accession"] == "GB_GCA_000001405.29"].iloc[0]
    assert bac_row["domain"] == "Bacteria"
    assert bac_row["species"] == "Escherichia coli"


def test_gen_gtdb_tab_combines_arc_and_bac(tmp_path):
    import pandas as pd

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download()), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        gen_gtdb_tab(str(tmp_path))

    result = pd.read_csv(tmp_path / "GTDB-arc-and-bac-metadata.tsv", sep="\t")
    assert len(result) == 2


def test_gen_gtdb_tab_bad_taxonomy_exits(tmp_path):
    bad_tsv = "accession\tgtdb_taxonomy\nRS_GCF_000002125.1\td__Archaea;p__Halobacteriota\n"

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download(arc_tsv=bad_tsv)), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        with pytest.raises(SystemExit):
            gen_gtdb_tab(str(tmp_path))


def test_gen_gtdb_tab_slims_to_kept_columns(tmp_path):
    """gen_gtdb_tab should drop everything outside GTDB_KEPT_COLUMNS, keeping
    the split rank columns and never leaking wide upstream-only columns."""
    import pandas as pd

    cols = ("accession\tgtdb_taxonomy\tncbi_genbank_assembly_accession\t"
            "genome_size\tncbi_organism_name\tssu_silva_taxonomy")
    arc = (cols + f"\nRS_GCF_000002125.1\t{ARC_TAXONOMY}\tGCA_000002125.1\t"
           "3890005\tHaloarcula sp\tsome_tax\n")
    bac = (cols + f"\nGB_GCA_000001405.29\t{BAC_TAXONOMY}\tGCA_000001405.29\t"
           "4600000\tE coli\tother_tax\n")

    with patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_make_mock_download(arc_tsv=arc, bac_tsv=bac)), \
         patch("bit.modules.gtdb.get_gtdb_data.urllib.request.urlretrieve"):
        gen_gtdb_tab(str(tmp_path))

    result = pd.read_csv(tmp_path / "GTDB-arc-and-bac-metadata.tsv", sep="\t")
    # only kept columns, in kept order, intersected with what was present
    expected = [c for c in GTDB_KEPT_COLUMNS if c in result.columns]
    assert list(result.columns) == expected
    # upstream-only wide columns and the raw taxonomy string are gone
    for dropped in ("gtdb_taxonomy", "ncbi_organism_name", "ssu_silva_taxonomy"):
        assert dropped not in result.columns
    # rank columns and present kept columns survive
    for kept in ("domain", "genus", "ncbi_genbank_assembly_accession", "genome_size"):
        assert kept in result.columns


################################################################################
# get_slim_gtdb_tab
################################################################################

def test_get_slim_gtdb_tab_extracts_both_files(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_slim_tarball(src, extra_member="README-junk.txt")
    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_mock_tarball_download(src)):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)

    assert (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert (tmp_path / "GTDB-version-info.txt").read_text().startswith("v220")
    # stray archive member is ignored, and the temp tarball is cleaned up
    assert not (tmp_path / "README-junk.txt").exists()
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


def test_get_slim_gtdb_tab_flattens_nested_metadata(tmp_path):
    """A metadata file nested under a directory in the archive is still written
    to the data dir root by basename (guards against path-prefixed archives)."""
    src = tmp_path / "src.tar.gz"
    _build_slim_tarball(src, metadata_under_dir=True)
    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_mock_tarball_download(src)):
        get_slim_gtdb_tab(str(tmp_path), quiet=True)

    assert (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert not (tmp_path / "some-dir").exists()


def test_get_slim_gtdb_tab_no_url_falls_back_to_rebuild(tmp_path):
    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", ""), \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    mock_gen.assert_called_once_with(str(tmp_path))


def test_get_slim_gtdb_tab_download_error_falls_back(tmp_path):
    def boom(url, label, dest, **kwargs):
        raise ConnectionError("server down")

    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", side_effect=boom), \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    mock_gen.assert_called_once_with(str(tmp_path))
    assert not (tmp_path / "GTDB-slim.tar.gz").exists()


def test_get_slim_gtdb_tab_missing_file_in_archive_falls_back(tmp_path):
    src = tmp_path / "src.tar.gz"
    _build_slim_tarball(src, include_version=False)  # metadata only, no version
    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_mock_tarball_download(src)), \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    mock_gen.assert_called_once_with(str(tmp_path))


################################################################################
# get_slim_gtdb_tab download success / fallback
################################################################################

def _copying_download(src_tarball):
    """download_with_tqdm side effect that copies the prebuilt tarball to dest and
    returns the dest path (retries now live inside download_with_tqdm itself, so
    from get_slim_gtdb_tab's view a call either lands the file or raises)."""
    def _download(url, label, dest, **kwargs):
        shutil.copy(src_tarball, dest)
        return dest
    return _download


def test_get_slim_gtdb_tab_download_success_uses_fast_path(tmp_path):
    """A successful download lands the tarball and completes via the fast path
    (no upstream rebuild)."""
    src = tmp_path / "src.tar.gz"
    _build_slim_tarball(src)
    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm",
               side_effect=_copying_download(src)), \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    mock_gen.assert_not_called()
    assert (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert (tmp_path / "GTDB-version-info.txt").exists()


def test_get_slim_gtdb_tab_download_failure_falls_back(tmp_path):
    """When the download raises (retries already exhausted inside
    download_with_tqdm), get_slim_gtdb_tab falls back to the upstream rebuild."""
    def always_fail(url, label, dest, **kwargs):
        raise socket.timeout("handshake timed out")

    with patch.object(gtdb_mod, "GTDB_SLIM_TARBALL_URL", "https://example.test/x.tar.gz"), \
         patch("bit.modules.gtdb.get_gtdb_data.download_with_tqdm", side_effect=always_fail), \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_slim_gtdb_tab(str(tmp_path), quiet=True)
    mock_gen.assert_called_once_with(str(tmp_path))


################################################################################
# get_gtdb_data
################################################################################

def test_get_gtdb_data_already_present_skips_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    with patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen, \
         patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_slim:
        get_gtdb_data(force_update=False, quiet=True)
    mock_gen.assert_not_called()
    mock_slim.assert_not_called()


def test_get_gtdb_data_default_uses_slim_path(tmp_path, monkeypatch):
    """Default (no -f) on missing data should use the fast slim path, not the
    upstream rebuild."""
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    with patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_slim, \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_gtdb_data(quiet=True)
    mock_slim.assert_called_once_with(str(tmp_path), quiet=True)
    mock_gen.assert_not_called()


def test_get_gtdb_data_force_update_forces_slim_fetch(tmp_path, monkeypatch):
    """-f re-fetches the hosted asset even when local data already exists; it
    does not rebuild from upstream directly (that's only the fallback)."""
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    with patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen, \
         patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_slim:
        get_gtdb_data(force_update=True, quiet=True)
    mock_slim.assert_called_once_with(str(tmp_path), quiet=True)
    mock_gen.assert_not_called()


def test_get_gtdb_data_missing_triggers_slim_path(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    with patch("bit.modules.gtdb.get_gtdb_data.get_slim_gtdb_tab") as mock_slim, \
         patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab"):
        get_gtdb_data()
    mock_slim.assert_called_once_with(str(tmp_path), quiet=False)


def test_get_gtdb_data_returns_gtdb_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    result = get_gtdb_data(force_update=False, quiet=True)
    assert result == str(tmp_path)


################################################################################
# main() bundle builder + fetch_upstream_version (for the refresh Action)
################################################################################

def test_main_builds_bundle_and_leaves_version_loose(tmp_path):
    def fake_gen(location):
        Path(location, "GTDB-arc-and-bac-metadata.tsv").write_text(
            "accession\tdomain\tgenus\nRS_1\tBacteria\tEscherichia\n")
        Path(location, "GTDB-version-info.txt").write_text(
            "Genome Taxonomy Database (GTDB) v232\nReleased April 2025\n")
    with patch.object(gtdb_mod, "gen_gtdb_tab", side_effect=fake_gen):
        rc = gtdb_mod.main(["-o", str(tmp_path),
                            "--archive-name", "GTDB-arc-and-bac-metadata.tar.gz"])
    assert rc == 0
    import tarfile
    tgz = tmp_path / "GTDB-arc-and-bac-metadata.tar.gz"
    assert tgz.exists()
    # large uncompressed table removed; loose version file kept for separate upload
    assert not (tmp_path / "GTDB-arc-and-bac-metadata.tsv").exists()
    assert (tmp_path / "GTDB-version-info.txt").exists()
    with tarfile.open(tgz, "r:gz") as t:
        members = sorted(m.name for m in t.getmembers())
    assert members == ["GTDB-arc-and-bac-metadata.tsv", "GTDB-version-info.txt"]


def test_fetch_upstream_version_returns_body():
    with patch("urllib.request.urlopen") as mock_open:
        mock_open.return_value.__enter__.return_value.read.return_value = \
            b"Genome Taxonomy Database (GTDB) v232\nReleased April 2025\n"
        v = gtdb_mod.fetch_upstream_version()
    assert "v232" in v
    assert "Released" in v
