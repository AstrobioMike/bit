import argparse
import pytest # type: ignore
import requests # type: ignore
import gzip
import time
from pathlib import Path
from unittest.mock import patch, MagicMock

from bit.modules.ncbi.dl_ncbi_assemblies import (
    RunData,
    setup,
    summarize_search,
    download_one,
    download_assemblies,
    report_finish,
    valid_gzip,
    sleep_backoff,
    download_one
)

def _write_gzip(path, payload):
    with open(path, "wb") as raw:
        pass
    with gzip.open(path, "wb") as fh:
        fh.write(payload)
 
 
def _truncate(path, fraction=0.5):
    data = Path(path).read_bytes()
    Path(path).write_bytes(data[: int(len(data) * fraction)])
 
 
def _mock_response(status=200, content_type="application/octet-stream",
                   body=b"data", headers=None):
    resp = MagicMock()
    resp.status_code = status
    resp.headers = {"Content-Type": content_type}
    if headers:
        resp.headers.update(headers)
    resp.raise_for_status.return_value = None
    resp.iter_content.return_value = [body]
    return resp


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


def test_valid_gzip_accepts_intact_file(tmp_path):
    p = tmp_path / "good.fasta.gz"
    _write_gzip(p, b"ACGT" * 100_000)  # ~400 KB payload, multi-block
    assert valid_gzip(str(p)) is True
 
 
def test_valid_gzip_rejects_small_truncation(tmp_path):
    # a small payload truncated mid-stream raises EOFError internally; the function
    # must treat that as invalid (not crash). EOFError is not an OSError subclass,
    # which is exactly the case a first-block-only read would miss.
    p = tmp_path / "small.fasta.gz"
    _write_gzip(p, b"ACGT" * 100)  # ~400 B payload
    _truncate(p)
    assert valid_gzip(str(p)) is False
 
 
def test_valid_gzip_rejects_large_truncation(tmp_path):
    # a large payload truncated mid-stream: reading only the first block would
    # succeed (the cut is far past block 1) and wrongly pass. The full read reaches
    # the missing trailer and correctly rejects it.
    p = tmp_path / "large.fasta.gz"
    _write_gzip(p, b"ACGT" * 1_000_000)  # ~4 MB payload
    _truncate(p)
    assert valid_gzip(str(p)) is False
 
 
def test_valid_gzip_rejects_corrupted_body(tmp_path):
    # valid gzip structure but a flipped mid-stream byte -> CRC32 mismatch at EOS
    p = tmp_path / "corrupt.fasta.gz"
    _write_gzip(p, b"ACGT" * 50_000)
    data = bytearray(p.read_bytes())
    data[len(data) // 2] ^= 0xFF
    p.write_bytes(bytes(data))
    assert valid_gzip(str(p)) is False
 
 
def test_valid_gzip_rejects_non_gzip_content(tmp_path):
    p = tmp_path / "junk.fasta.gz"
    p.write_bytes(b"this is not gzip data at all")
    assert valid_gzip(str(p)) is False
 
 
def test_valid_gzip_passes_through_non_gz_paths(tmp_path):
    # non-.gz paths are not decompressed; the function returns True without reading
    p = tmp_path / "report.txt"
    p.write_text("some plain text")
    assert valid_gzip(str(p)) is True
 
 
# ---------------------------------------------------------------------------
# Tier 2a - sleep_backoff
# ---------------------------------------------------------------------------
 
def test_sleep_backoff_honors_retry_after_header():
    resp = MagicMock()
    resp.headers = {"Retry-After": "0"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep:
        sleep_backoff(1, resp=resp)
    # slept exactly once, using the header value (0.0)
    mock_sleep.assert_called_once_with(0.0)
 
 
def test_sleep_backoff_invalid_retry_after_falls_back_to_exponential():
    resp = MagicMock()
    resp.headers = {"Retry-After": "not-a-number"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep, \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        sleep_backoff(3, resp=resp)
    # 2 ** (3 - 1) + 0.0 == 4.0
    mock_sleep.assert_called_once_with(4.0)
 
 
def test_sleep_backoff_no_response_uses_exponential():
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep, \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        sleep_backoff(1)
    # 2 ** 0 + 0.0 == 1.0
    mock_sleep.assert_called_once_with(1.0)
 
 
# ---------------------------------------------------------------------------
# Tier 2b - download_one retry branches
# ---------------------------------------------------------------------------
 
def test_download_one_transient_then_success(tmp_path):
    dest = str(tmp_path / "a.gz")
    seq = [_mock_response(status=503), _mock_response(status=200)]
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", side_effect=seq), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.sleep_backoff") as mock_backoff:
        path, err, status = download_one("http://x/a.gz", dest, retries=3)
    assert status == "downloaded"
    assert err is None
    # backed off once between the 503 and the 200
    assert mock_backoff.call_count == 1
 
 
def test_download_one_transient_exhausts_retries(tmp_path):
    dest = str(tmp_path / "b.gz")
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get",
               side_effect=[_mock_response(status=503)] * 3), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.sleep_backoff"):
        path, err, status = download_one("http://x/b.gz", dest, retries=3)
    assert status == "failed"
    assert "503" in err
    assert "3 attempts" in err
 
 
def test_download_one_error_page_then_success(tmp_path):
    dest = str(tmp_path / "c.gz")
    seq = [_mock_response(content_type="text/html"),
           _mock_response(content_type="application/octet-stream", body=b"data")]
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", side_effect=seq), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.sleep_backoff") as mock_backoff:
        path, err, status = download_one("http://x/c.gz", dest, retries=3)
    assert status == "downloaded"
    assert mock_backoff.call_count == 1
 
 
def test_download_one_error_page_exhausts_retries(tmp_path):
    dest = str(tmp_path / "d.gz")
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get",
               side_effect=[_mock_response(content_type="application/xml")] * 2), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.sleep_backoff"):
        path, err, status = download_one("http://x/d.gz", dest, retries=2)
    assert status == "failed"
    assert "error page" in err
 
 
def test_download_one_empty_then_success(tmp_path):
    dest = str(tmp_path / "e.gz")
    seq = [_mock_response(body=b""), _mock_response(body=b"data")]
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", side_effect=seq), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.sleep_backoff") as mock_backoff:
        path, err, status = download_one("http://x/e.gz", dest, retries=3)
    assert status == "downloaded"
    # the first (empty) attempt removes the zero-byte file and backs off once
    assert mock_backoff.call_count == 1
