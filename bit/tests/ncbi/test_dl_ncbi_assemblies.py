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
    sleep_backoff,
    download_one,
    MAX_BACKOFF,
    MAX_RETRY_AFTER
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
    assert status == "failed_transient"
    assert "error page" in err


def test_download_one_xml_error_page(tmp_path):
    dest = str(tmp_path / "file.gz")
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.raise_for_status.return_value = None
    mock_resp.headers = {"Content-Type": "application/xml"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get", return_value=mock_resp):
        path, err, status = download_one("http://fake/file.gz", dest, retries=1)
    assert status == "failed_transient"
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
    assert status == "failed_transient"
    assert err == "Downloaded file was empty"
    assert not Path(dest).exists()


def test_download_one_request_exception(tmp_path):
    dest = str(tmp_path / "file.gz")
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.requests.get",
               side_effect=requests.RequestException("timeout")):
        path, err, status = download_one("http://fake/file.gz", dest, retries=2)
    assert status == "failed_transient"
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


# ---------------------------------------------------------------------------
# Tier 2a - sleep_backoff
# ---------------------------------------------------------------------------

def test_sleep_backoff_honors_retry_after_header():
    resp = MagicMock()
    resp.headers = {"Retry-After": "0"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep:
        sleep_backoff(1, resp=resp, throttled=True)
    # slept exactly once, using the header value (0.0)
    mock_sleep.assert_called_once_with(0.0)


def test_sleep_backoff_ignores_retry_after_when_not_throttled():
    # a plain 5xx/connection blip takes the sawtooth path, which does not consult
    # Retry-After at all (only an explicit throttle does)
    resp = MagicMock()
    resp.headers = {"Retry-After": "600"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep, \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        sleep_backoff(1, resp=resp)
    mock_sleep.assert_called_once_with(1.0)


def test_sleep_backoff_sawtooth_resets_each_cycle():
    # non-throttle sleeps cycle 1, 2, 4, 8, 16 then start over, so a straggler never
    # parks a pool thread for hours
    seen = []
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep", side_effect=seen.append), \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        for attempt in range(1, 12):
            sleep_backoff(attempt)
    assert seen == [1.0, 2.0, 4.0, 8.0, 16.0, 1.0, 2.0, 4.0, 8.0, 16.0, 1.0]


def test_sleep_backoff_throttled_is_capped():
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep, \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        sleep_backoff(20, throttled=True)
    mock_sleep.assert_called_once_with(MAX_BACKOFF)


def test_sleep_backoff_caps_absurd_retry_after():
    resp = MagicMock()
    resp.headers = {"Retry-After": "99999"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep:
        sleep_backoff(1, resp=resp, throttled=True)
    mock_sleep.assert_called_once_with(MAX_RETRY_AFTER)


def test_sleep_backoff_honors_retry_after_beyond_max_backoff():
    # Retry-After gets the larger ceiling: a server saying "come back in 120s" is
    # obeyed rather than clamped to MAX_BACKOFF, which would just burn retries on
    # requests we already know will be refused
    resp = MagicMock()
    resp.headers = {"Retry-After": "120"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep:
        sleep_backoff(1, resp=resp, throttled=True)
    mock_sleep.assert_called_once_with(120.0)


def test_sleep_backoff_invalid_retry_after_falls_back_to_exponential():
    resp = MagicMock()
    resp.headers = {"Retry-After": "not-a-number"}
    with patch("bit.modules.ncbi.dl_ncbi_assemblies.time.sleep") as mock_sleep, \
         patch("bit.modules.ncbi.dl_ncbi_assemblies.random.uniform", return_value=0.0):
        sleep_backoff(3, resp=resp, throttled=True)
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
    assert status == "failed_transient"
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
    assert status == "failed_transient"
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


# ─── taxon-based targeting (`-t`), deduplication, and --dry-run ───────────────

from types import SimpleNamespace
from unittest.mock import patch as _patch

from bit.modules.ncbi.dl_ncbi_assemblies import (
    resolve_targets, report_selection, dl_ncbi_assemblies, TaxonSelection)
from bit.modules.taxonomy.tax_select import CrossDomainTaxon, AmbiguousTaxon


def _dl_args(**kwargs):
    defaults = dict(
        wanted_accessions=None, target_taxon=None, source="ncbi",
        ncbi_section="both", derep_rank="auto", target_rank=None,
        target_domain=None, assembly_level=None, min_completeness=None,
        max_contamination=None, representatives_only=False,
        dry_run=False, format="fasta", jobs=10, output_dir=".", quiet=False,
    )
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


def _fake_selection(canonical, accessions, rank="phylum", derep="family"):
    return SimpleNamespace(canonical=canonical, resolved_rank=rank,
                           effective_derep_rank=derep, accessions=list(accessions),
                           rows=[], warnings=[])


class TestResolveTargets:

    def test_accessions_file_alone_still_works(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCA_001\nGCA_002\n")
        got, selections, note = resolve_targets(_dl_args(wanted_accessions=str(accs)))
        assert got == ["GCA_001", "GCA_002"]
        assert selections == [] and note is None

    def test_duplicate_lines_in_the_file_are_collapsed(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCA_001\nGCA_002\nGCA_001\n")
        got, _, _ = resolve_targets(_dl_args(wanted_accessions=str(accs)))
        assert got == ["GCA_001", "GCA_002"]

    def test_a_single_taxon_resolves_to_its_accessions(self):
        sel = _fake_selection("Nitrospirota", ["GCF_1", "GCF_2"])
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    return_value=(sel.accessions, sel)), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["Nitrospirota"], [])):
            got, selections, _ = resolve_targets(_dl_args(target_taxon=["Nitrospirota"]))
        assert got == ["GCF_1", "GCF_2"]
        assert selections[0].num_selected == 2
        assert selections[0].num_new == 2

    def test_repeated_taxa_dedupe_and_report_the_overlap(self):
        """
        Two overlapping taxa must contribute each shared genome once, and the per-taxon
        counts must still show what each SELECTED -- that gap is the point of the report.
        """
        first = _fake_selection("Bacteria", ["GCF_1", "GCF_2"])
        second = _fake_selection("Escherichia", ["GCF_2", "GCF_3"])
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    side_effect=[(first.accessions, first), (second.accessions, second)]), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["Bacteria", "Escherichia"], [])):
            got, selections, _ = resolve_targets(
                _dl_args(target_taxon=["Bacteria", "Escherichia"]))

        assert got == ["GCF_1", "GCF_2", "GCF_3"]
        assert (selections[0].num_selected, selections[0].num_new) == (2, 2)
        assert (selections[1].num_selected, selections[1].num_new) == (2, 1)
        # the total is the union, NOT the sum of the per-taxon counts
        assert len(got) == 3 != sum(s.num_selected for s in selections)

    def test_taxon_accessions_dedupe_against_the_accessions_file(self, tmp_path):
        accs = tmp_path / "accs.txt"
        accs.write_text("GCF_1\n")
        sel = _fake_selection("Bacteria", ["GCF_1", "GCF_2"])
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    return_value=(sel.accessions, sel)), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["Bacteria"], [])):
            got, selections, _ = resolve_targets(
                _dl_args(wanted_accessions=str(accs), target_taxon=["Bacteria"]))
        assert got == ["GCF_1", "GCF_2"]
        assert selections[0].num_new == 1

    def test_a_cross_domain_name_exits_pointing_at_target_domain(self, capsys):
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["Bacillus"], [])), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    side_effect=CrossDomainTaxon("Bacillus", ["Bacteria", "Eukaryota"])):
            with pytest.raises(SystemExit):
                resolve_targets(_dl_args(target_taxon=["Bacillus"]))
        assert "--target-domain" in capsys.readouterr().out

    def test_an_ambiguous_rank_exits_pointing_at_dash_r(self, capsys):
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["X"], [])), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    side_effect=AmbiguousTaxon("X", ["order", "family"])):
            with pytest.raises(SystemExit):
                resolve_targets(_dl_args(target_taxon=["X"]))
        assert "`-r`" in capsys.readouterr().out


class TestSelectionKwargsPlumbing:

    def _capture(self, **args):
        sel = _fake_selection("T", ["GCF_1"])
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    return_value=(sel.accessions, sel)) as m, \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["T"], [])):
            resolve_targets(_dl_args(target_taxon=["T"], **args))
        return m.call_args.kwargs

    def test_reps_only_is_off_by_default_for_both_sources(self):
        """
        NOT the source's own default (which would make GTDB representatives-only).
        Derep already controls volume here, and get-accs-from-gtdb defaults to all
        genomes too, so both sources start unrestricted.
        """
        assert self._capture()["reps_only"] is False
        assert self._capture(source="gtdb")["reps_only"] is False

    def test_representatives_only_restricts_either_source(self):
        assert self._capture(representatives_only=True)["reps_only"] is True
        assert self._capture(source="gtdb", representatives_only=True)["reps_only"] is True

    def test_quality_floor_and_levels_are_passed_through(self):
        kw = self._capture(min_completeness=90.0, max_contamination=5.0,
                           assembly_level=["Complete Genome"])
        assert kw["min_completeness"] == 90.0
        assert kw["max_contamination"] == 5.0
        assert kw["assembly_levels"] == ["Complete Genome"]

    def test_section_and_derep_are_passed_through(self):
        kw = self._capture(ncbi_section="genbank", derep_rank="genus")
        assert kw["ncbi_section"] == "genbank"
        assert kw["derep_rank"] == "genus"


class TestDryRun:

    def _run_dry(self, tmp_path, **args):
        first = _fake_selection("Bacteria", ["GCF_1", "GCF_2"])
        second = _fake_selection("Escherichia", ["GCF_2", "GCF_3"])
        a = _dl_args(target_taxon=["Bacteria", "Escherichia"], dry_run=True,
                     output_dir=str(tmp_path / "out"), **args)
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies.resolve_target_taxon_accessions",
                    side_effect=[(first.accessions, first), (second.accessions, second)]), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["Bacteria", "Escherichia"], [])), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.get_ncbi_assembly_data"), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.ensure_source_data"), \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.download_assemblies") as dl:
            dl_ncbi_assemblies(a)
        return dl, a

    def test_dry_run_downloads_nothing(self, tmp_path, capsys):
        dl, _ = self._run_dry(tmp_path)
        capsys.readouterr()
        dl.assert_not_called()

    def test_dry_run_writes_no_files_and_no_output_dir(self, tmp_path, capsys):
        _, a = self._run_dry(tmp_path)
        capsys.readouterr()
        assert not Path(a.output_dir).exists()

    def test_dry_run_reports_per_taxon_counts_and_the_deduped_total(self, tmp_path, capsys):
        self._run_dry(tmp_path)
        out = capsys.readouterr().out
        assert "Bacteria" in out and "Escherichia" in out
        # union of {1,2} and {2,3} is 3, not 4
        assert "3" in out
        assert "already counted" in out


class TestReportSelection:

    def test_a_clean_run_reports_no_overlap_wording(self, capsys):
        sels = [TaxonSelection(taxon="A", canonical="A", resolved_rank="phylum",
                               effective_derep_rank="family", num_selected=5,
                               num_new=5, warnings=[])]
        report_selection(["a"] * 5, sels, (None, 0), _dl_args())
        out = capsys.readouterr().out
        assert "5 genome(s) selected" in out
        assert "already counted" not in out

    def test_derep_off_is_labelled(self, capsys):
        sels = [TaxonSelection(taxon="A", canonical="A", resolved_rank="species",
                               effective_derep_rank=None, num_selected=2,
                               num_new=2, warnings=[])]
        report_selection(["a", "b"], sels, (None, 0), _dl_args())
        assert "dereplication off" in capsys.readouterr().out

    def test_the_all_expansion_note_is_surfaced(self, capsys):
        sels = [TaxonSelection(taxon="Bacteria", canonical="Bacteria",
                               resolved_rank="domain", effective_derep_rank="class",
                               num_selected=1, num_new=1, warnings=[])]
        report_selection(["a"], sels, ("`-t all` was expanded to: Bacteria, Archaea", 0),
                         _dl_args())
        assert "expanded to" in capsys.readouterr().out

    def test_quiet_suppresses_the_report(self, capsys):
        sels = [TaxonSelection(taxon="A", canonical="A", resolved_rank="phylum",
                               effective_derep_rank="family", num_selected=5,
                               num_new=5, warnings=[])]
        report_selection(["a"] * 5, sels, (None, 0), _dl_args(quiet=True))
        assert capsys.readouterr().out == ""


class TestNotFoundMessagingBranchesOnInputMode:

    def test_accession_input_says_inputs_may_be_invalid(self):
        rd = RunData(from_taxon=False)
        assert "may be invalid" in rd.not_found_reason
        assert "assembly accessions?" in rd.none_found_hint

    def test_taxon_input_does_not_blame_the_user_or_suggest_a_refresh(self):
        """
        A `-t` user never typed an accession, so "may be invalid" is wrong, and there
        is nothing for them to go re-download. State it plainly, like GToTree does.
        """
        rd = RunData(from_taxon=True)
        assert rd.not_found_reason == ""
        assert "bit data get" not in rd.none_found_hint
        assert "invalid" not in rd.none_found_hint


class TestDryRunCountMatchesTheRealSelection:
    """
    The dry run's whole value is that its number is the number. It must come from the
    same selection call the real run uses -- not a cheaper distinct-group count, which
    doesn't see liveness screening or the quality floor (see
    test_domain_aware_derep.TestCheapCountPathWouldOverreport).
    """

    def test_dry_run_and_real_run_resolve_identically(self, tmp_path):
        sel = _fake_selection("Bacteria", ["GCF_1", "GCF_2", "GCF_3"])

        def run(dry):
            a = _dl_args(target_taxon=["Bacteria"], dry_run=dry,
                         output_dir=str(tmp_path / ("dry" if dry else "real")))
            with _patch(f"bit.modules.ncbi.dl_ncbi_assemblies."
                        f"resolve_target_taxon_accessions",
                        return_value=(sel.accessions, sel)), \
                 _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                        return_value=(["Bacteria"], [])):
                return resolve_targets(a)

        dry_accs, dry_sels, _ = run(True)
        real_accs, real_sels, _ = run(False)
        assert dry_accs == real_accs
        assert dry_sels[0].num_selected == real_sels[0].num_selected

    def test_the_selection_call_skips_metadata_rows(self):
        """
        dl-ncbi-assemblies downloads files rather than writing a metadata TSV, so it
        asks for accessions only -- one less filtered read of the asset per taxon.
        """
        sel = _fake_selection("T", ["GCF_1"])
        with _patch("bit.modules.ncbi.dl_ncbi_assemblies."
                    "resolve_target_taxon_accessions",
                    return_value=(sel.accessions, sel)) as m, \
             _patch("bit.modules.ncbi.dl_ncbi_assemblies.expand_target_taxa",
                    return_value=(["T"], [])):
            resolve_targets(_dl_args(target_taxon=["T"]))
        assert m.call_args.kwargs["include_rows"] is False
