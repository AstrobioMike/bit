"""
Tests for bit.modules.ref_dbs: where the reference databases are, whether they're set up, and the setup steps (with
the downloading itself replaced, see test_downloads.py for that).
"""

import json
from pathlib import Path
from unittest.mock import patch

import pytest # type: ignore

from bit.modules import ref_dbs
from bit.modules.downloads import DownloadError
from bit.modules.ref_dbs import REF_DBS, RefDbsDirNotSet, missing, ref_dbs_root, setup


def test_the_directory_given_wins_over_the_variable(monkeypatch, tmp_path):
    monkeypatch.setenv("BIT_REF_DBS_DIR", str(tmp_path / "var"))
    assert ref_dbs_root() == tmp_path / "var"
    assert ref_dbs_root(str(tmp_path / "given")) == tmp_path / "given"
    monkeypatch.delenv("BIT_REF_DBS_DIR")
    with pytest.raises(RefDbsDirNotSet, match="bit data locations set"):
        ref_dbs_root()


def test_an_unpacked_tax_db_without_its_marker_isnt_set_up(tmp_path):
    db = REF_DBS["tax-db"]
    db.path(tmp_path).mkdir()
    for name in db.ready[:-1]:
        (db.path(tmp_path) / name).touch()
    assert missing(tmp_path, ["tax-db"]) == ["tax-db"]      # e.g. an unpack cut short in its last file


def test_a_database_is_set_up_when_all_its_marker_files_are_there(tmp_path):
    db = REF_DBS["checkm2-db"]
    assert missing(tmp_path, ["checkm2-db", "kofamscan-db"]) == ["checkm2-db", "kofamscan-db"]
    db.path(tmp_path).mkdir()
    (db.path(tmp_path) / "uniref100.KO.1.dmnd").touch()
    assert missing(tmp_path, ["checkm2-db"]) == ["checkm2-db"]         # downloaded, but not marked as finished
    (db.path(tmp_path) / "CHECKM2_DB_SETUP").touch()
    assert missing(tmp_path, ["checkm2-db", "kofamscan-db"]) == ["kofamscan-db"]


def fake_amplicon_fetch(downloads, label, **kwargs):
    for url, path, md5 in downloads:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(url)
    return ["fetched"]


def test_setup_marks_it_done_and_skips_it_after(tmp_path, capsys):
    with patch("bit.modules.downloads.fetch", side_effect=fake_amplicon_fetch) as fetch:
        setup("amplicon-refs", tmp_path)
        assert fetch.call_count == 1
        assert missing(tmp_path, ["amplicon-refs"]) == []
        setup("amplicon-refs", tmp_path)
        assert fetch.call_count == 1                                     # already set up
        assert "already set up" in capsys.readouterr().out
        setup("amplicon-refs", tmp_path, force=True)
        assert fetch.call_count == 2


def test_a_failed_setup_exits_with_a_message_and_isnt_marked(tmp_path, capsys):
    with patch("bit.modules.downloads.fetch", side_effect=DownloadError("Getting x failed: 503")):
        with pytest.raises(SystemExit):
            setup("kofamscan-db", tmp_path)
    out = " ".join(capsys.readouterr().out.split())
    assert "Getting x failed: 503" in out and "carries on from where it stopped" in out
    assert missing(tmp_path, ["kofamscan-db"]) == ["kofamscan-db"]


def test_tax_db_parts_come_from_the_md5_list_and_are_unpacked_in_order(tmp_path):
    md5_list = b"aa  bit-tax-db-2026-11.tar.zst.part-000\nbb  bit-tax-db-2026-11.tar.zst.part-001\n"
    calls = {}

    def fake_fetch(downloads, label, **kwargs):
        calls["fetch"] = downloads
        for _, path, _ in downloads:
            Path(path).write_text("x")
        return []

    def fake_unpack(paths, into, label, decompress, strip_components=0):
        calls["unpack"] = (paths, decompress, strip_components)
        for name in REF_DBS["tax-db"].ready[:-1]:     # the data files; the marker is setup's to write, last
            (Path(into) / name).touch()

    with patch("bit.modules.downloads.read_url", return_value=md5_list), \
         patch("bit.modules.downloads.fetch", side_effect=fake_fetch), \
         patch("bit.modules.downloads.unpack", side_effect=fake_unpack):
        setup("tax-db", tmp_path)

    db_dir = tmp_path / "bit-tax-db-2026-11"
    assert calls["fetch"] == [(ref_dbs.TAX_DB_URL + f"bit-tax-db-2026-11.tar.zst.part-00{i}",
                               str(db_dir / f"bit-tax-db-2026-11.tar.zst.part-00{i}"), md5) for i, md5 in ((0, "aa"), (1, "bb"))]
    assert calls["unpack"][1:] == (["zstd", "-dc"], 1)
    assert not list(db_dir.glob("*.part-*"))                            # the parts are removed once unpacked
    assert missing(tmp_path, ["tax-db"]) == [] and (db_dir / "TAX_DB_SETUP").exists()


def test_checkm2_and_gtdbtk_md5s_come_from_their_hosts(tmp_path):
    record = json.dumps({"files": [{"key": "checkm2_database.tar.gz", "checksum": "md5:abc123"}]}).encode()
    gtdb_list = (b"111  ./auxillary_files/gtdbtk_package/split_package/gtdbtk_r232_data.tar.gz.part_aa\n"
                 b"222  ./auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz\n")
    seen = []

    def fake_fetch(downloads, label, **kwargs):
        seen.append(downloads[0][2])
        Path(downloads[0][1]).write_text("x")
        return []

    def fake_unpack(paths, into, *args, **kwargs):
        (Path(into) / "uniref100.KO.1.dmnd").touch()

    with patch("bit.modules.downloads.read_url", side_effect=[record, gtdb_list]), \
         patch("bit.modules.downloads.fetch", side_effect=fake_fetch), \
         patch("bit.modules.downloads.unpack", side_effect=fake_unpack), \
         patch.object(ref_dbs, "GTDBTK_RUNS_HERE", True):
        setup("checkm2-db", tmp_path)
        setup("gtdbtk-db", tmp_path)

    assert seen == ["abc123", "222"]
    assert missing(tmp_path, ["checkm2-db", "gtdbtk-db"]) == []
    assert "zenodo" in (tmp_path / "checkm2-ref-db-1.1.0" / "DB-INFO.txt").read_text()
