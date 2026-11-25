"""
Tests for the check made before anything is downloaded (bit.modules.connectivity), and for where it is and isn't
made. The network is replaced throughout: `probe` (whether one site answers) is set per test, apart from the tests of
probe itself, which ask a server started here.
"""

import sys
import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest # type: ignore

from bit.cli import run_wf
from bit.modules import connectivity, ref_dbs
from bit.modules import snakemake as smk
from bit.modules import workflows as wf
from bit.modules.connectivity import GENERAL_SITES, NO_CONNECTION, require_internet
from bit.tests.utils import require_local_connections

REAL_PROBE = connectivity.probe      # the autouse fixture (conftest.py) replaces it for every test

TRY_AGAIN = ("You can try again later, or you can get the db from another machine that can access the site(s) and "
             "then move the data here.")


def answering(monkeypatch, *down, nothing=False):
    """ sets which sites answer (all but those of `down`, or none at all); returns the list of sites asked """
    asked = []

    def probe(url, timeout=None):
        asked.append(url)
        return not nothing and connectivity.host(url) not in down

    monkeypatch.setattr(connectivity, "probe", probe)
    return asked


def said(capsys):
    captured = capsys.readouterr()
    return " ".join((captured.out + captured.err).split())


def run_cli(monkeypatch, *argv):
    if argv[0] in wf.WORKFLOW_COMMANDS:
        monkeypatch.setattr(sys, "argv", [f"bit {argv[0]}", *argv[1:]])
        run_wf.command_main(argv[0])
    else:
        monkeypatch.setattr(sys, "argv", ["bit run-wf", *argv])
        run_wf.main()


@pytest.fixture
def genomes_dir(tmp_path):
    d = tmp_path / "genomes"
    d.mkdir()
    (d / "g1.fasta").write_text(">c1\nACGT\n")
    return d


def set_up_dbs(root, names):
    for name in names:
        db = ref_dbs.REF_DBS[name]
        db.path(root).mkdir(parents=True, exist_ok=True)
        for f in db.ready:
            (db.path(root) / f).touch()


############################################################
######################## the check #########################
############################################################

def test_no_connection_at_all_stops_with_that_message(monkeypatch, capsys):
    asked = answering(monkeypatch, nothing=True)
    with pytest.raises(SystemExit) as e:
        require_internet(["https://ftp.ncbi.nlm.nih.gov/genomes/all/x.fna.gz"])
    assert e.value.code == 1
    out = said(capsys)
    assert NO_CONNECTION in out
    assert NO_CONNECTION == ("This would require downloading something from the internet, but it seems there is no "
                             "connection right now. So we can't proceed :(")
    assert "ftp.ncbi.nlm.nih.gov" not in out and "try again later" not in out
    # the site needed and the always-up ones were all asked
    assert set(asked) == {"https://ftp.ncbi.nlm.nih.gov", *GENERAL_SITES}


def test_connected_but_one_site_down_names_it(monkeypatch, capsys):
    answering(monkeypatch, "ftp.ncbi.nlm.nih.gov")
    with pytest.raises(SystemExit):
        require_internet(["https://ftp.ncbi.nlm.nih.gov/a", "https://github.com/b"], what="the db")
    out = said(capsys)
    assert ("The internet is reachable, but ftp.ncbi.nlm.nih.gov isn't answering right now; it may be down, or "
            "blocked from this network.") in out
    assert out.count("github.com") == 0
    assert TRY_AGAIN in out
    assert NO_CONNECTION not in out


def test_several_sites_down_are_all_named_and_the_wording_fits_what_is_downloaded(monkeypatch, capsys):
    answering(monkeypatch, "zenodo.org", "ftp.ebi.ac.uk")
    with pytest.raises(SystemExit):
        require_internet(["https://zenodo.org/a", "https://ftp.ebi.ac.uk/b", "https://github.com/c"], what="the data")
    out = said(capsys)
    assert "but zenodo.org and ftp.ebi.ac.uk aren't answering right now; they may be down" in out
    assert "you can get the data from another machine that can access the site(s) and then move the data here." in out


def test_the_sites_needed_answering_is_enough(monkeypatch, capsys):
    # e.g. a network that lets through only what it has to: the always-up sites not answering doesn't matter
    answering(monkeypatch, *[connectivity.host(s) for s in GENERAL_SITES])
    require_internet(["https://ftp.ncbi.nlm.nih.gov/a"])
    assert said(capsys) == ""


def test_a_slow_site_is_given_longer_before_it_is_said_to_be_down(monkeypatch, capsys):
    asked = []

    def probe(url, timeout=None):
        asked.append((url, timeout))
        return connectivity.host(url) != "ftp.ebi.ac.uk" or timeout == connectivity.SLOW_TIMEOUT

    monkeypatch.setattr(connectivity, "probe", probe)
    require_internet(["https://ftp.ebi.ac.uk/a", "https://github.com/b"])
    assert said(capsys) == ""
    assert [t for url, t in asked if url == "https://ftp.ebi.ac.uk"] == [connectivity.TIMEOUT, connectivity.SLOW_TIMEOUT]
    assert [t for url, t in asked if url == "https://github.com"] == [connectivity.TIMEOUT]

    # with no connection at all there is nothing to wait longer for
    asked.clear()
    monkeypatch.setattr(connectivity, "probe", lambda url, timeout=None: asked.append(timeout))
    with pytest.raises(SystemExit):
        require_internet(["https://ftp.ebi.ac.uk/a"])
    assert set(asked) == {connectivity.TIMEOUT}


def test_each_site_is_asked_once_whatever_the_number_of_urls(monkeypatch):
    asked = answering(monkeypatch)
    require_internet(["https://zenodo.org/api/records/1", "https://zenodo.org/api/records/1/files/x/content"])
    assert asked.count("https://zenodo.org") == 1


def test_nothing_to_download_means_nothing_is_asked(monkeypatch, capsys):
    asked = answering(monkeypatch, nothing=True)
    require_internet([])
    assert asked == [] and said(capsys) == ""


def test_site_and_host_come_from_the_url():
    assert connectivity.site(ref_dbs.EUKCC_DL_LINK) == "https://ftp.ebi.ac.uk"
    assert connectivity.site("http://purl.obolibrary.org/obo/go/go-basic.obo") == "http://purl.obolibrary.org"
    assert connectivity.host("http://localhost:8080/x") == "localhost"
    assert connectivity.site("github.com") == "https://github.com"


############################################################
################# whether a site answers ###################
############################################################

def serve(statuses, drop_head=False):
    """ a local server answering each method with the status given for it; returns (its URL, the methods asked) """
    require_local_connections()
    asked = []

    class Handler(BaseHTTPRequestHandler):
        def answer(self):
            asked.append(self.command)
            if drop_head and self.command == "HEAD":
                # no response at all to a HEAD request (shut down first, so the client sees it end at once on
                # every system: a close alone can leave the socket open while the handler's files refer to it)
                self.close_connection = True
                try:
                    self.connection.shutdown(socket.SHUT_RDWR)
                except OSError:
                    pass
                self.connection.close()
                return
            self.send_response(statuses[self.command])
            self.send_header("Content-Length", "0")
            self.end_headers()

        do_HEAD = do_GET = answer

        def log_message(self, *args):
            pass

    server = HTTPServer(("127.0.0.1", 0), Handler)
    # (shutdown() waits for the serving loop's next look, half a second apart unless told otherwise)
    threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_port}/some/file", asked


@pytest.fixture
def no_proxy(monkeypatch):
    for name in ("http_proxy", "HTTP_PROXY", "https_proxy", "HTTPS_PROXY", "all_proxy", "ALL_PROXY"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")


@pytest.mark.parametrize("status", [200, 301, 403, 404, 405, 500, 503])
def test_any_http_status_counts_as_reachable(status, no_proxy):
    server, url, asked = serve({"HEAD": status, "GET": status})
    try:
        assert REAL_PROBE(url, timeout=5) is True
        assert asked == ["HEAD"]                 # an error status is an answer: nothing more is asked
    finally:
        server.shutdown()
        server.server_close()


def test_get_is_tried_when_head_gets_no_response(no_proxy):
    server, url, asked = serve({"GET": 404}, drop_head=True)
    try:
        assert REAL_PROBE(url, timeout=5) is True
        assert asked == ["HEAD", "GET"]
    finally:
        server.shutdown()
        server.server_close()


def test_no_answer_means_unreachable(no_proxy):
    # nothing listens on this port any more (connection refused), and this name doesn't resolve
    server, url, _ = serve({"HEAD": 200, "GET": 200})
    server.shutdown()
    server.server_close()
    assert REAL_PROBE(url, timeout=5) is False
    assert REAL_PROBE("https://no-such-host.invalid/x", timeout=5) is False


def test_the_real_probe_drives_the_messages(no_proxy, monkeypatch, capsys):
    """ end to end with real requests: one local site answering 404, one refusing connections """
    up, up_url, _ = serve({"HEAD": 404, "GET": 404})
    gone, gone_url, _ = serve({"HEAD": 200, "GET": 200})
    gone.shutdown()
    gone.server_close()
    monkeypatch.setattr(connectivity, "probe", REAL_PROBE)
    monkeypatch.setattr(connectivity, "GENERAL_SITES", (up_url,))
    try:
        require_internet([up_url])                                 # a 404 is an answer
        with pytest.raises(SystemExit):
            require_internet([up_url, gone_url.replace("127.0.0.1", "localhost")], what="the data")
        assert "The internet is reachable, but localhost isn't answering right now" in said(capsys)
    finally:
        up.shutdown()
        up.server_close()


############################################################
############### missing reference databases ################
############################################################

def test_missing_dbs_are_listed_with_their_commands_in_both_messages(monkeypatch, tmp_path, capsys):
    names = ["tax-db", "checkm2-db"]
    urls = ref_dbs.download_urls(names)
    for nothing, down in ((True, ()), (False, ("zenodo.org",))):
        answering(monkeypatch, *down, nothing=nothing)
        with pytest.raises(SystemExit):
            require_internet(urls, what="the db", missing_dbs=names, ref_dbs_dir=tmp_path)
        out = said(capsys)
        assert f"This run needs these reference databases, which aren't set up yet in {tmp_path}/:" in out
        for name in names:
            assert f"{name} {ref_dbs.REF_DBS[name].size}" in out
            assert f"bit data get {name} -d <dir>" in out
        assert "They can be downloaded on another machine with internet access and then copied here" in out
        assert "then point BIT_REF_DBS_DIR (or -d/--ref-dbs-dir) at that directory." in out
        assert "gtdbtk-db" not in out
        if nothing:
            assert NO_CONNECTION in out and TRY_AGAIN not in out
        else:
            assert "but zenodo.org isn't answering right now" in out
            assert out.index("bit data get tax-db") < out.index(TRY_AGAIN)      # the message ends with it


def test_every_database_says_where_it_is_downloaded_from():
    assert set(ref_dbs.DOWNLOAD_URLS) == set(ref_dbs.REF_DBS)
    hosts = {name: {connectivity.host(u) for u in urls} for name, urls in ref_dbs.DOWNLOAD_URLS.items()}
    assert hosts == {
        "tax-db": {"github.com"}, "checkm2-db": {"zenodo.org"}, "gtdbtk-db": {"data.ace.uq.edu.au"},
        "kofamscan-db": {"github.com"}, "eukcc-db": {"ftp.ebi.ac.uk"}, "amplicon-refs": {"github.com"},
    }


def test_data_get_checks_only_when_it_would_download(monkeypatch, tmp_path, capsys):
    asked = answering(monkeypatch, nothing=True)

    # already set up: nothing to download, nothing asked
    set_up_dbs(tmp_path, ["eukcc-db"])
    ref_dbs.setup("eukcc-db", tmp_path)
    assert asked == []

    # not set up: stops before anything is fetched, saying how to get it elsewhere
    with patch("bit.modules.downloads.fetch") as fetch, pytest.raises(SystemExit):
        ref_dbs.setup("checkm2-db", tmp_path)
    fetch.assert_not_called()
    out = said(capsys)
    assert NO_CONNECTION in out
    assert f"This reference database isn't set up yet in {tmp_path}/:" in out
    assert "bit data get checkm2-db -d <dir>" in out
    assert "https://zenodo.org" in asked and "https://ftp.ebi.ac.uk" not in asked

    # forced: it would be downloaded again, so it is checked, and what is there is left alone when that can't be done
    with pytest.raises(SystemExit):
        ref_dbs.setup("eukcc-db", tmp_path, force=True)
    assert ref_dbs.missing(tmp_path, ["eukcc-db"]) == []


def test_set_up_missing_checks_all_of_them_before_the_first_download(monkeypatch, tmp_path, capsys):
    answering(monkeypatch, "ftp.ebi.ac.uk")
    with patch("bit.modules.downloads.fetch") as fetch, pytest.raises(SystemExit):
        ref_dbs.set_up_missing(tmp_path, ["amplicon-refs", "eukcc-db"])
    fetch.assert_not_called()
    out = said(capsys)
    assert "but ftp.ebi.ac.uk isn't answering right now" in out
    assert "bit data get amplicon-refs -d <dir>" in out and "bit data get eukcc-db -d <dir>" in out


############################################################
###################### the workflows #######################
############################################################

ENV_LINE = "Conda environment {} will be created.\n"


def planned(*env_files, rules=()):
    """ stands in for the dry run made before a run starts: the environments snakemake says it still has to create,
    and the steps it has to do """
    return patch("bit.modules.workflows.dry_run_report",
                 return_value="Building DAG of jobs...\n" + "".join(ENV_LINE.format(f) for f in env_files) +
                              "".join(f"rule {rule}:\n    resources: tmpdir=/tmp\n" for rule in rules))


@pytest.fixture
def conda_defaults(monkeypatch):
    monkeypatch.setattr(smk, "_conda_channel_settings", lambda: (smk.CONDA_CHANNEL_ALIAS, {}))
    monkeypatch.delenv("PIP_INDEX_URL", raising=False)


def env_file(name):
    return str(Path(wf.get_package_path("smk/envs")) / name)


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_a_run_with_everything_in_place_asks_nothing(setup, snakemake, monkeypatch, tmp_path, genomes_dir, capsys):
    asked = answering(monkeypatch, nothing=True)
    set_up_dbs(tmp_path / "dbs", ["tax-db"])
    with planned():
        run_cli(monkeypatch, "classify-seqs", "-f", str(genomes_dir), "-d", str(tmp_path / "dbs"), "-o", str(tmp_path / "out"))
    assert asked == []
    setup.assert_not_called()
    snakemake.assert_called_once()


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_a_run_missing_a_database_stops_before_anything_starts(setup, snakemake, monkeypatch, tmp_path, genomes_dir, capsys):
    monkeypatch.setattr("bit.modules.ref_dbs.GTDBTK_RUNS_HERE", True)
    asked = answering(monkeypatch, nothing=True)
    set_up_dbs(tmp_path / "dbs", ["checkm2-db"])
    with planned(), pytest.raises(SystemExit):
        run_cli(monkeypatch, "characterize-genomes", "-g", str(genomes_dir), "-d", str(tmp_path / "dbs"),
                "-o", str(tmp_path / "out"))
    out = said(capsys)
    assert NO_CONNECTION in out
    # only what this run uses and doesn't have: GTDB-Tk's, not CheckM2's (set up) or the others (not used)
    assert f"This run needs this reference database, which isn't set up yet in {tmp_path / 'dbs'}/:" in out
    assert "bit data get gtdbtk-db -d <dir>" in out
    assert "checkm2-db" not in out and "tax-db" not in out
    assert "https://data.ace.uq.edu.au" in asked and "https://zenodo.org" not in asked
    setup.assert_not_called()
    snakemake.assert_not_called()


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_the_databases_listed_are_the_ones_then_set_up(setup, snakemake, monkeypatch, tmp_path, genomes_dir, capsys):
    """ the message's list and what is set up before snakemake starts come from the same place """
    monkeypatch.setattr("bit.modules.ref_dbs.GTDBTK_RUNS_HERE", True)
    given = ("characterize-genomes", "-g", str(genomes_dir), "--euk", "-d", str(tmp_path / "dbs"), "-o", str(tmp_path / "out"))

    answering(monkeypatch, "ftp.ebi.ac.uk")
    with planned(), pytest.raises(SystemExit):
        run_cli(monkeypatch, *given)
    out = said(capsys)
    listed = [name for name in ref_dbs.REF_DBS if f"bit data get {name} -d <dir>" in out]

    answering(monkeypatch)
    with planned():
        run_cli(monkeypatch, *given)
    assert sorted(listed) == sorted(call.args[0] for call in setup.call_args_list) == ["eukcc-db", "tax-db"]


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_setup_dbs_only_checks_too(setup, snakemake, monkeypatch, tmp_path, capsys):
    answering(monkeypatch, nothing=True)
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, "classify-seqs", "--setup-dbs-only", "-d", str(tmp_path / "dbs"))
    assert "bit data get tax-db -d <dir>" in said(capsys)
    setup.assert_not_called()


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_amplicon_references_headed_for_the_output_directory_are_said_to_be(setup, snakemake, monkeypatch, tmp_path, capsys):
    monkeypatch.delenv("BIT_REF_DBS_DIR", raising=False)
    answering(monkeypatch, nothing=True)
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, "amplicon", "--setup-dbs-only", "-o", str(tmp_path / "amp"))
    out = said(capsys)
    assert (f"This run needs this reference database, which isn't set up yet in {tmp_path / 'amp' / 'ref-dbs'}/ "
            "(no reference-database directory is set, so it would go in the output directory):") in out
    assert "bit data get amplicon-refs -d <dir>" in out
    assert "then point BIT_REF_DBS_DIR (or -d/--ref-dbs-dir) at that directory." in out
    assert not (tmp_path / "amp").exists()


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_a_first_run_checks_the_conda_channels(setup, snakemake, conda_defaults, monkeypatch, tmp_path, genomes_dir, capsys):
    set_up_dbs(tmp_path / "dbs", ["tax-db"])
    given = ("classify-seqs", "-f", str(genomes_dir), "-d", str(tmp_path / "dbs"), "-o", str(tmp_path / "out"))

    asked = answering(monkeypatch, "conda.anaconda.org")
    with planned(env_file("classify-seqs.yaml")), pytest.raises(SystemExit):
        run_cli(monkeypatch, *given)
    out = said(capsys)
    assert "The internet is reachable, but conda.anaconda.org isn't answering right now" in out
    assert "You can try again later, or run this from a machine or network that can access the site(s)." in out
    assert "reference database" not in out            # none is missing
    # how the environments could be created from another machine sharing this filesystem
    assert ("The conda environments this needs (classify-seqs) aren't set up yet. If this machine shares its "
            "filesystem with one that can access the site(s), running the same command there creates them.") in out
    assert "--conda-create-envs-only" not in out      # classify-seqs takes no snakemake pass-through
    assert "https://conda.anaconda.org" in asked
    snakemake.assert_not_called()

    # once snakemake has the environments, nothing is asked
    asked = answering(monkeypatch, nothing=True)
    with planned():
        run_cli(monkeypatch, *given)
    assert asked == []
    snakemake.assert_called_once()


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_run_wf_environments_note_says_how_to_create_just_them(setup, snakemake, conda_defaults, monkeypatch, tmp_path,
                                                               capsys):
    reads = tmp_path / "reads"
    reads.mkdir()
    for r in ("s1_R1.fastq.gz", "s1_R2.fastq.gz"):
        (reads / r).touch()
    set_up_dbs(tmp_path / "dbs", ["amplicon-refs"])
    answering(monkeypatch, nothing=True)
    with planned(env_file("qc.yaml")), pytest.raises(SystemExit):
        run_cli(monkeypatch, "amplicon", "-r", str(reads), "-d", str(tmp_path / "dbs"), "-o", str(tmp_path / "out"),
                "--conda-prefix", str(tmp_path / "envs"))
    out = said(capsys)
    assert NO_CONNECTION in out
    assert ("The conda environments this needs (qc) aren't set up yet. If this machine shares its filesystem with "
            "one that can access the site(s), adding `-- --conda-create-envs-only` to the same command there creates "
            "just them.") in out
    snakemake.assert_not_called()


def test_assemble_says_how_to_bring_its_environment_from_another_machine(monkeypatch, tmp_path, capsys):
    from types import SimpleNamespace
    from bit.modules import assemble
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path / "conda"))
    monkeypatch.setattr(smk, "_conda_channel_settings", lambda: (smk.CONDA_CHANNEL_ALIAS, {}))
    created = []
    monkeypatch.setattr(assemble, "run", lambda cmd, *a, **k: created.append(cmd))
    answering(monkeypatch, nothing=True)
    with pytest.raises(SystemExit):
        assemble.check_conda_env(SimpleNamespace(dry_run=False))
    out = said(capsys)
    env_path = f"{tmp_path / 'conda'}/envs/bit-assemble"
    assert NO_CONNECTION in out
    assert ("The conda environment this needs (bit-assemble) isn't set up yet. If you can get it there yourself, "
            "you could use or match its yaml (") in out
    # (the yaml is assemble.yaml, or assemble-osx-64.yaml on macOS)
    assert f".yaml) and install it here: {env_path}" in out
    assert created == []                                # conda env create never ran

    # with the environment there, nothing is asked
    (tmp_path / "conda" / "envs" / "bit-assemble").mkdir(parents=True)
    asked = answering(monkeypatch, nothing=True)
    assemble.check_conda_env(SimpleNamespace(dry_run=False))
    assert asked == [] and created == []


@patch("bit.modules.workflows.run_snakemake")
def test_dl_sra_reads_checks_ncbi(snakemake, monkeypatch, tmp_path, capsys):
    accs = tmp_path / "accs.txt"
    accs.write_text("SRR390728\n")
    given = ("dl-sra-reads", "-a", str(accs), "-o", str(tmp_path / "out"))

    answering(monkeypatch, "locate.ncbi.nlm.nih.gov")
    with planned(rules=["prefetch", "fasterq_dump_and_gzip", "all"]), pytest.raises(SystemExit):
        run_cli(monkeypatch, *given)
    out = said(capsys)
    assert "but locate.ncbi.nlm.nih.gov isn't answering right now" in out
    assert "you can get the data from another machine" in out
    snakemake.assert_not_called()

    asked = answering(monkeypatch)
    with planned(rules=["prefetch", "fasterq_dump_and_gzip", "all"]):
        run_cli(monkeypatch, *given)
    assert asked.count("https://locate.ncbi.nlm.nih.gov") == 1
    snakemake.assert_called_once()

    # nothing left to fetch (e.g. run again once done): nothing asked
    asked = answering(monkeypatch, nothing=True)
    with planned(rules=["all"]):
        run_cli(monkeypatch, *given)
    assert asked == []


def test_the_rules_said_to_download_exist_in_their_workflows():
    for name in wf.WORKFLOW_DIRS:
        workflow = wf.load_workflow(name)
        snakefile = (workflow.directory / "Snakefile").read_text()
        for rule in workflow.downloads_from:
            assert f"\nrule {rule}:" in snakefile, (name, rule)


@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.setup")
def test_a_dry_run_checks_nothing_and_says_what_would_need_the_internet(setup, snakemake, conda_defaults, monkeypatch,
                                                                         tmp_path, genomes_dir, capsys):
    asked = answering(monkeypatch, nothing=True)
    with planned(env_file("classify-seqs.yaml")):
        run_cli(monkeypatch, "classify-seqs", "-f", str(genomes_dir), "-d", str(tmp_path / "dbs"),
                "-o", str(tmp_path / "out"), "--dry-run")
    out = said(capsys)
    assert asked == []
    assert "DRY-RUN COMPLETE!" in out
    assert f"These reference databases would be set up first in {tmp_path / 'dbs'}/: tax-db" in out
    assert "These conda environments would be created first: classify-seqs" in out
    assert "The run would need the internet, to download from: github.com and conda.anaconda.org" in out
    setup.assert_not_called()


############################################################
################### the conda channels #####################
############################################################

def test_planned_conda_envs_reads_what_snakemake_says_it_will_create():
    report = ("Building DAG of jobs...\n"
              "Conda environment /a/b/envs/qc.yaml will be created.\n"
              "Incomplete Conda environment /a b/envs/mapping.yaml will be recreated.\n"
              "Conda environment /a/b/envs/qc.yaml will be created.\n"
              "Job stats:\n")
    assert smk.planned_conda_envs(report) == ["/a/b/envs/qc.yaml", "/a b/envs/mapping.yaml"]
    assert smk.planned_conda_envs("Building DAG of jobs...\nNothing to be done.\n") == []
    assert smk.planned_conda_envs(None) == []        # the dry run failed: the real run reports why


def test_conda_env_urls_follow_the_environment_files_and_the_conda_settings(monkeypatch, tmp_path):
    monkeypatch.delenv("PIP_INDEX_URL", raising=False)
    monkeypatch.setattr(smk, "_conda_channel_settings", lambda: (smk.CONDA_CHANNEL_ALIAS, {}))
    assert smk.conda_env_urls([]) == []
    assert smk.conda_env_urls([env_file("qc.yaml")]) == ["https://conda.anaconda.org"]
    # an environment with a `pip:` section also gets packages from PyPI
    assert smk.conda_env_urls([env_file("keggdecoder.yaml")]) == ["https://conda.anaconda.org", "https://pypi.org"]

    # a channel given as a URL is itself; a mirror set in conda's configuration is used for the channels it covers
    own = tmp_path / "own.yaml"
    own.write_text("channels:\n  - https://repo.example.org/pkgs\n  - conda-forge\n  - bioconda\ndependencies:\n  - x\n")
    monkeypatch.setattr(smk, "_conda_channel_settings",
                        lambda: ("https://mirror.example.edu/conda", {"bioconda": "https://bio.example.edu/channels"}))
    assert [connectivity.host(u) for u in smk.conda_env_urls([str(own)])] == \
           ["repo.example.org", "mirror.example.edu", "bio.example.edu"]


def test_every_workflow_environment_file_names_its_channels():
    """ so the first-run check knows where each is created from """
    monkey = lambda: (smk.CONDA_CHANNEL_ALIAS, {})
    with patch.object(smk, "_conda_channel_settings", monkey):
        for path in sorted(Path(wf.get_package_path("smk/envs")).glob("*.yaml")):
            assert "https://conda.anaconda.org" in smk.conda_env_urls([str(path)]), path.name


############################################################
################### the other downloads ####################
############################################################

def test_hosted_tables_check_their_site_only_when_fetched(monkeypatch, tmp_path, capsys):
    from bit.modules.ncbi.get_ncbi_assembly_data import NCBI_ASSET, get_ncbi_assembly_data
    monkeypatch.setenv("NCBI_assembly_data_dir", str(tmp_path))
    asked = answering(monkeypatch, nothing=True)

    with patch("bit.modules.hosted_parquet_asset.download_with_tqdm") as download, pytest.raises(SystemExit):
        get_ncbi_assembly_data(quiet=True)
    download.assert_not_called()
    assert NO_CONNECTION in said(capsys)
    assert connectivity.site(NCBI_ASSET.data_url) in asked

    # cached: no check
    for name in (NCBI_ASSET.parquet_filename, NCBI_ASSET.sidecar_filename):
        (tmp_path / name).write_text("x")
    asked.clear()
    get_ncbi_assembly_data(quiet=True)
    assert asked == []


def test_ncbi_tax_and_go_data_check_their_sites(monkeypatch, tmp_path, capsys):
    from bit.modules.go.get_go_dbs import download_go_data
    from bit.modules.ncbi.get_ncbi_tax_data import download_ncbi_tax_data

    answering(monkeypatch, "ftp.ncbi.nih.gov", "current.geneontology.org")
    with patch("bit.modules.ncbi.get_ncbi_tax_data.download_with_tqdm") as download, pytest.raises(SystemExit):
        download_ncbi_tax_data(str(tmp_path))
    download.assert_not_called()
    assert "but ftp.ncbi.nih.gov isn't answering right now" in said(capsys)

    with patch("bit.modules.go.get_go_dbs.download_with_tqdm") as download, pytest.raises(SystemExit):
        download_go_data(str(tmp_path))
    download.assert_not_called()
    assert "but current.geneontology.org isn't answering right now" in said(capsys)


def test_test_data_checks_the_site_it_comes_from(monkeypatch, tmp_path, capsys):
    from bit.modules.get_test_data import dl_test_data
    monkeypatch.chdir(tmp_path)
    answering(monkeypatch, "ftp.ncbi.nlm.nih.gov")
    with patch("bit.modules.get_test_data.subprocess.run") as run, pytest.raises(SystemExit):
        dl_test_data(SimpleNamespace(datatype="genome"))
    run.assert_not_called()
    out = said(capsys)
    assert "but ftp.ncbi.nlm.nih.gov isn't answering right now" in out
    assert "you can get the data from another machine" in out


def test_dl_ncbi_assemblies_checks_before_downloading_unless_all_are_here(monkeypatch, tmp_path, capsys):
    from bit.modules.ncbi import dl_ncbi_assemblies as dl
    table = tmp_path / "info.tsv"
    here, absent = tmp_path / "GCF_1.fasta.gz", tmp_path / "GCF_2.fasta.gz"
    link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/x/x_genomic.fna.gz"
    here.write_text("x")
    run_data = SimpleNamespace(quiet=True, ncbi_sub_table_path=table, num_jobs=1)

    asked = answering(monkeypatch, nothing=True)
    with patch.object(dl, "run_download_pass", return_value=([], [], 1)) as passes:
        # everything already downloaded: nothing asked
        table.write_text(f"target_link\tlocal_destination\n{link}\t{here}\n")
        try:
            dl.download_assemblies(run_data)
        except Exception:
            pass                                   # what follows the downloads isn't under test here
        assert asked == [] and passes.call_count == 1

        table.write_text(f"target_link\tlocal_destination\n{link}\t{here}\n{link}\t{absent}\n")
        with pytest.raises(SystemExit):
            dl.download_assemblies(run_data)
        assert passes.call_count == 1
    assert NO_CONNECTION in said(capsys)
    assert "https://ftp.ncbi.nlm.nih.gov" in asked


def test_assemble_dry_run_installs_nothing_and_says_it_would(monkeypatch, tmp_path, capsys):
    from types import SimpleNamespace
    from bit.modules import assemble
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path / "bit"))
    monkeypatch.setattr(smk, "_conda_channel_settings", lambda: (smk.CONDA_CHANNEL_ALIAS, {}))
    created = []
    monkeypatch.setattr(assemble, "run", lambda cmd, *a, **k: created.append(cmd))
    asked = answering(monkeypatch, nothing=True)
    config = assemble.check_conda_env(SimpleNamespace(dry_run=True))
    out = said(capsys)
    env_path = f"{tmp_path / 'bit'}/envs/bit-assemble"
    assert "The conda environment bit-assemble is not set up yet. The run would first create it at:" in out
    assert env_path in out and "downloading from" not in out
    assert created == [] and asked == []                # nothing installed, nothing asked
    assert config.conda_env == env_path
