import gzip
import io
import pytest
from pathlib import Path
from unittest.mock import patch

from bit.modules.gtdb.get_gtdb_data import (
    check_gtdb_location_var_is_set,
    check_if_gtdb_data_present,
    gen_gtdb_tab,
    get_gtdb_data,
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
    def _download(link, label, dest):
        if "ar53" in link:
            Path(dest).write_bytes(_to_gz(arc_tsv))
        else:
            Path(dest).write_bytes(_to_gz(bac_tsv))
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


################################################################################
# get_gtdb_data
################################################################################

def test_get_gtdb_data_already_present_skips_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    with patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_gtdb_data(force_update=False, quiet=True)
    mock_gen.assert_not_called()


def test_get_gtdb_data_force_update_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    with patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_gtdb_data(force_update=True)
    mock_gen.assert_called_once_with(str(tmp_path))


def test_get_gtdb_data_missing_triggers_download(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    with patch("bit.modules.gtdb.get_gtdb_data.gen_gtdb_tab") as mock_gen:
        get_gtdb_data()
    mock_gen.assert_called_once_with(str(tmp_path))


def test_get_gtdb_data_returns_gtdb_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("GTDB_DIR", str(tmp_path))
    (tmp_path / "GTDB-arc-and-bac-metadata.tsv").write_text("data")
    (tmp_path / "GTDB-version-info.txt").write_text("v220")
    result = get_gtdb_data(force_update=False, quiet=True)
    assert result == str(tmp_path)
