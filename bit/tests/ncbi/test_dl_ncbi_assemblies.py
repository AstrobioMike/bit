import argparse
import pytest # type: ignore
import requests # type: ignore
from pathlib import Path
from unittest.mock import patch, MagicMock

from bit.modules.ncbi.dl_ncbi_assemblies import (
    RunData,
    setup,
    summarize_search,
    download_one,
    download_assemblies,
    report_finish,
)


def test_rundata_defaults():
    rd = RunData()
    assert rd.num_jobs == 10
    assert rd.num_wanted == 0
    assert rd.num_found == 0
    assert rd.num_not_found == 0
    assert rd.num_downloaded == 0
    assert rd.num_skipped == 0
    assert rd.num_not_downloaded == 0


def test_setup_reads_accessions(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCA_001\nGCA_002\nGCA_003\n")
    args = argparse.Namespace(
        wanted_accessions=str(accs),
        format="fasta",
        jobs=4,
        output_dir=str(tmp_path),
    )
    rd = setup(args)
    assert rd.wanted_accs == ["GCA_001", "GCA_002", "GCA_003"]
    assert rd.num_wanted == 3
    assert rd.wanted_format == "fasta"


def test_setup_ignores_blank_lines(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCA_001\n\nGCA_002\n\n")
    args = argparse.Namespace(
        wanted_accessions=str(accs),
        format="fasta",
        jobs=4,
        output_dir=str(tmp_path),
    )
    rd = setup(args)
    assert rd.wanted_accs == ["GCA_001", "GCA_002"]
    assert rd.num_wanted == 2


def test_setup_output_paths(tmp_path):
    accs = tmp_path / "accs.txt"
    accs.write_text("GCA_001\n")
    args = argparse.Namespace(
        wanted_accessions=str(accs),
        format="fasta",
        jobs=4,
        output_dir=str(tmp_path),
    )
    rd = setup(args)
    assert rd.ncbi_sub_table_path == tmp_path / "wanted-ncbi-accessions-info.tsv"
    assert rd.not_found_path == tmp_path / "ncbi-accessions-not-found.txt"
    assert rd.not_downloaded_path == tmp_path / "ncbi-accessions-not-downloaded.tsv"


def test_summarize_search_all_found():
    rd = RunData(num_wanted=3, num_found=3)
    summarize_search(rd)  # should complete without error or exit


def test_summarize_search_partial_found(tmp_path, capsys):
    not_found = tmp_path / "not-found.txt"
    rd = RunData(num_wanted=3, num_found=2, num_not_found=1, not_found_path=not_found)
    summarize_search(rd)
    assert "1 accession(s) not found" in capsys.readouterr().out


def test_summarize_search_none_found_exits(tmp_path):
    not_found = tmp_path / "not-found.txt"
    not_found.write_text("")  # must exist so os.remove doesn't raise
    rd = RunData(num_wanted=3, num_found=0, not_found_path=not_found)
    with pytest.raises(SystemExit):
        summarize_search(rd)


def test_download_one_success(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.raise_for_status.return_value = None
    mock_resp.headers = {"Content-Type": "application/octet-stream"}
    mock_resp.iter_content.return_value = [b"somedata"]
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest)
    assert err is None
    assert path == dest
    assert status == "downloaded"
    assert Path(dest).exists()


def test_download_one_skips_existing_valid_file(tmp_path):
    # a non-gz file that already exists and is non-empty should be skipped
    # without any network call (valid_gzip returns True for non-.gz paths)
    dest = tmp_path / "file.fasta"
    dest.write_text("already here")
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get") as mock_get:
        path, err, status = download_one("http://fake/file.fasta", str(dest))
    assert status == "skipped"
    assert err is None
    mock_get.assert_not_called()


def test_download_one_404(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 404
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest, retries=1)
    assert status == "failed"
    assert "404" in err
    assert not Path(dest).exists()


def test_download_one_html_error_page(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.raise_for_status.return_value = None
    mock_resp.headers = {"Content-Type": "text/html"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest, retries=1)
    assert status == "failed"
    assert "error page" in err


def test_download_one_xml_error_page(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.raise_for_status.return_value = None
    mock_resp.headers = {"Content-Type": "application/xml"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest, retries=1)
    assert status == "failed"
    assert "error page" in err


def test_download_one_empty_file(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.raise_for_status.return_value = None
    mock_resp.headers = {"Content-Type": "application/octet-stream"}
    mock_resp.iter_content.return_value = [b""]  # writes 0 bytes
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest, retries=1)
    assert status == "failed"
    assert err == "Downloaded file was empty"
    assert not Path(dest).exists()


def test_download_one_request_exception(tmp_path):
    dest = str(tmp_path / "file.gz")
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get",
               side_effect=requests.RequestException("timeout")):
        path, err, status = download_one("http://fake/file.gz", dest, retries=2)
    assert status == "failed"
    assert "timeout" in err
    assert not Path(dest).exists()


def _make_run_data(tmp_path, links_and_dests):
    tsv = tmp_path / "wanted-ncbi-accessions-info.tsv"
    rows = "\n".join(f"{link}\t{dest}" for link, dest in links_and_dests)
    tsv.write_text(f"target_link\tlocal_destination\n{rows}\n")
    return RunData(
        ncbi_sub_table_path=tsv,
        not_downloaded_path=tmp_path / "failed.tsv",
        num_jobs=1,
    )


def test_download_assemblies_all_succeed(tmp_path):
    dest = str(tmp_path / "a.gz")
    rd = _make_run_data(tmp_path, [("http://fake/a.gz", dest)])
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.download_one",
               return_value=(dest, None, "downloaded")):
        rd = download_assemblies(rd)
    assert rd.num_downloaded == 1
    assert rd.num_not_downloaded == 0
    assert rd.num_skipped == 0


def test_download_assemblies_all_fail(tmp_path):
    dest = str(tmp_path / "a.gz")
    rd = _make_run_data(tmp_path, [("http://fake/a.gz", dest)])
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.download_one",
               return_value=(dest, "timeout", "failed")):
        rd = download_assemblies(rd)
    assert rd.num_not_downloaded == 1
    assert rd.num_downloaded == 0
    assert (tmp_path / "failed.tsv").exists()


def test_download_assemblies_partial_failure(tmp_path):
    dest_a = str(tmp_path / "a.gz")
    dest_b = str(tmp_path / "b.gz")
    rd = _make_run_data(tmp_path, [
        ("http://fake/a.gz", dest_a),
        ("http://fake/b.gz", dest_b),
    ])
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.download_one",
               side_effect=[(dest_a, None, "downloaded"), (dest_b, "404", "failed")]):
        rd = download_assemblies(rd)
    assert rd.num_downloaded == 1
    assert rd.num_not_downloaded == 1


def test_download_assemblies_counts_skipped_as_downloaded(tmp_path):
    # skipped files are already present, so they count toward num_downloaded
    # and are tracked separately in num_skipped
    dest_a = str(tmp_path / "a.gz")
    dest_b = str(tmp_path / "b.gz")
    rd = _make_run_data(tmp_path, [
        ("http://fake/a.gz", dest_a),
        ("http://fake/b.gz", dest_b),
    ])
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.download_one",
               side_effect=[(dest_a, None, "downloaded"), (dest_b, None, "skipped")]):
        rd = download_assemblies(rd)
    assert rd.num_downloaded == 2
    assert rd.num_skipped == 1
    assert rd.num_not_downloaded == 0


def test_download_assemblies_removes_stale_failed_file(tmp_path):
    # a not-downloaded TSV left over from a prior run should be removed when the
    # current run has no failures (recovered-on-rerun behavior)
    dest = str(tmp_path / "a.gz")
    rd = _make_run_data(tmp_path, [("http://fake/a.gz", dest)])
    stale = tmp_path / "failed.tsv"
    stale.write_text("accession\terror\nGCA_999\told failure\n")
    assert stale.exists()
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.download_one",
               return_value=(dest, None, "downloaded")):
        rd = download_assemblies(rd)
    assert not stale.exists()
    assert rd.num_not_downloaded == 0


def test_report_finish_all_downloaded(capsys):
    rd = RunData(num_wanted=3, num_found=3, num_downloaded=3)
    report_finish(rd)
    assert "All 3" in capsys.readouterr().out


def test_report_finish_all_found_downloaded(capsys):
    rd = RunData(num_wanted=5, num_found=3, num_downloaded=3)
    report_finish(rd)
    assert "All 3 found" in capsys.readouterr().out


def test_report_finish_some_failed(capsys, tmp_path):
    rd = RunData(
        num_wanted=3, num_found=3,
        num_downloaded=2, num_not_downloaded=1,
        not_downloaded_path=tmp_path / "failed.tsv",
    )
    report_finish(rd)
    out = capsys.readouterr().out
    assert "1 file(s) failed" in out
    assert "2" in out


def test_report_finish_skipped_note(capsys):
    rd = RunData(num_wanted=3, num_found=3, num_downloaded=3, num_skipped=2)
    report_finish(rd)
    out = capsys.readouterr().out
    assert "already present" in out
