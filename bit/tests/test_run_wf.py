import re
import sys
import shutil
import subprocess
import pytest # type: ignore
import yaml # type: ignore
from unittest.mock import patch
from bit.cli import run_wf
from bit.modules import workflows as wf
from bit.modules.snakemake import build_snakemake_cmd


@pytest.fixture
def accs(tmp_path):
    f = tmp_path / "accs.txt"
    f.write_text("SRR390728\nERR000001\n")
    return f


@pytest.fixture
def reads_dir(tmp_path):
    d = tmp_path / "reads"
    d.mkdir()
    for sample in ("s1", "s2"):
        for r in ("R1", "R2"):
            (d / f"{sample}_{r}.fastq.gz").touch()
    return d


@pytest.fixture
def genomes_dir(tmp_path):
    d = tmp_path / "genomes"
    d.mkdir()
    for genome in ("g1", "g2"):
        (d / f"{genome}.fasta").write_text(">c1\nACGT\n")
    return d


@pytest.fixture
def assemblies_dir(tmp_path):
    """ an assembly per sample of the reads_dir fixture """
    d = tmp_path / "assemblies"
    d.mkdir()
    for sample in ("s1", "s2"):
        (d / f"{sample}.fasta").write_text(">k141_1\nACGT\n")
    return d


@pytest.fixture
def ref_dbs(tmp_path):
    """A reference-database directory where every database looks set up (their marker files, empty)."""
    from bit.modules.ref_dbs import REF_DBS
    root = tmp_path / "ref-dbs"
    for db in REF_DBS.values():
        db.path(root).mkdir(parents=True)
        for name in db.ready:
            (db.path(root) / name).touch()
    return root


def run_cli(monkeypatch, *argv):
    """Runs a packaged workflow as typed after `bit`: under run-wf, or as its own command."""
    if argv[0] in wf.WORKFLOW_COMMANDS:
        monkeypatch.setattr(sys, "argv", [f"bit {argv[0]}", *argv[1:]])
        run_wf.command_main(argv[0])
    else:
        monkeypatch.setattr(sys, "argv", ["bit run-wf", *argv])
        run_wf.main()


def workflow_parser(name, show_fine=False):
    if name in wf.WORKFLOW_COMMANDS:
        return run_wf.command_build_parser(name, show_fine=show_fine)
    return run_wf.build_parser(show_fine=show_fine).subparsers.choices[name]


def dry_runs(*runs):
    """
    Makes the given snakemake dry runs at the same time (each takes a second or two, mostly snakemake starting), and
    returns their results in the order given. A run is its command, or (command, settings for subprocess.run).
    """
    from concurrent.futures import ThreadPoolExecutor
    runs = [run if isinstance(run, tuple) else (run, {}) for run in runs]
    with ThreadPoolExecutor(len(runs)) as pool:
        return list(pool.map(lambda run: subprocess.run(run[0], capture_output=True, text=True, **run[1]), runs))


@pytest.mark.parametrize("name", list(wf.WORKFLOW_DIRS))
def test_parameter_tables_are_consistent(name):
    workflow = wf.load_workflow(name)
    names = [p.name for p in workflow.params]
    assert len(names) == len(set(names))
    shorts = [p.short for p in workflow.params if p.short]
    assert len(shorts) == len(set(shorts))
    # short flags used by every workflow's run options
    assert not set(shorts) & {"-o", "-c", "-F", "-H", "-h", "-v", "-j"}
    for p in workflow.params:
        assert p.kind in ("value", "flag", "file", "dir", "path", "paths")
        assert p.help.strip()
        if p.kind == "flag":
            assert p.default is False
    assert (name in run_wf.WORKFLOW_SUMMARIES) == (name in wf.WORKFLOWS)
    for src, _ in wf.workflow_files(workflow):
        assert src.is_file(), f"{src} is listed for the {name} workflow but doesn't exist"


@pytest.mark.parametrize("name", list(wf.WORKFLOW_DIRS))
def test_concise_help_hides_what_detailed_help_shows(name, capsys):
    workflow = wf.load_workflow(name)
    concise = workflow_parser(name).format_help()
    detailed = workflow_parser(name, show_fine=True).format_help()
    for p in workflow.params:
        assert p.flag in detailed
        assert (p.flag in concise) == p.main


def test_detailed_help_flag(monkeypatch, capsys):
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, "dl-sra-reads", "-H")
    assert e.value.code == 0
    assert "--prefetch-max-size" in capsys.readouterr().err


def test_settings_precedence(tmp_path, accs):
    workflow = wf.load_workflow("dl-sra-reads")
    config = tmp_path / "my.yaml"
    config.write_text(f"accessions_file: {accs.name}\nthreads: 3\nkeep_sra_files: true\n")

    settings = wf.resolve_settings(workflow.params, {"threads": 5}, config)

    assert settings["threads"] == 5                      # the flag beats the config file
    assert settings["keep_sra_files"] is True            # the config file beats the default
    assert settings["prefetch_max_size"] == "500G"       # the default
    assert settings["accessions_file"] == str(accs)      # relative to the config file, stored absolute


def test_missing_required_setting_exits(capsys):
    workflow = wf.load_workflow("dl-sra-reads")
    with pytest.raises(SystemExit):
        wf.resolve_settings(workflow.params, {})
    assert "--accessions-file" in capsys.readouterr().out


def test_unknown_config_key_and_bad_values_exit(tmp_path, accs, capsys):
    workflow = wf.load_workflow("dl-sra-reads")
    bad = tmp_path / "bad.yaml"
    bad.write_text("not_a_setting: 1\n")
    with pytest.raises(SystemExit):
        wf.resolve_settings(workflow.params, {"accessions_file": str(accs)}, bad)
    assert "not_a_setting" in capsys.readouterr().out

    bad.write_text("threads: many\n")
    with pytest.raises(SystemExit):
        wf.resolve_settings(workflow.params, {"accessions_file": str(accs)}, bad)

    with pytest.raises(SystemExit):
        wf.resolve_settings(workflow.params, {"accessions_file": str(tmp_path / "nope.txt")})


def test_written_config_reads_back_the_same(tmp_path, accs):
    workflow = wf.load_workflow("dl-sra-reads")
    settings = wf.resolve_settings(workflow.params, {"accessions_file": str(accs), "keep_sra_files": True})
    out = tmp_path / "config.yaml"
    wf.write_config(workflow, settings, out, command="bit dl-sra-reads -a accs.txt")

    assert yaml.safe_load(out.read_text()) == settings
    assert wf.resolve_settings(workflow.params, {}, out) == settings
    # every setting carries its help text as a comment
    assert "prefetch" in out.read_text()


def test_get_workflow_writes_a_complete_copy(monkeypatch, tmp_path, reads_dir, capsys):
    dest = tmp_path / "saved"
    given = ("metagenomics", "-r", str(reads_dir), "-d", str(tmp_path / "ref-dbs"))
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, *given, "-t", "2", "--get-workflow", str(dest))
    assert e.value.code == 0

    for rel in ("Snakefile", "config.yaml", "rules/common.smk", "envs/checkm2.yaml"):
        assert (dest / rel).is_file()
    assert yaml.safe_load((dest / "config.yaml").read_text())["threads"] == 2

    # an empty directory that's already there is fine
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, *given, "--get-workflow", str(empty))
    assert e.value.code == 0
    assert (empty / "Snakefile").is_file()


def test_get_workflow_never_writes_into_a_directory_holding_files(monkeypatch, tmp_path, reads_dir, capsys):
    given = ("metagenomics", "-r", str(reads_dir), "-d", str(tmp_path / "ref-dbs"))

    # the user's own directory: refused, and what's in it is left exactly as it was
    mine = tmp_path / "my-project"
    mine.mkdir()
    (mine / "Snakefile").write_text("my own rules\n")
    (mine / "config.yaml").write_text("mine: true\n")
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, *given, "--get-workflow", str(mine))
    assert e.value.code == 1
    assert "name a different directory" in " ".join(capsys.readouterr().out.split())
    assert (mine / "Snakefile").read_text() == "my own rules\n"
    assert (mine / "config.yaml").read_text() == "mine: true\n"
    assert sorted(p.name for p in mine.iterdir()) == ["Snakefile", "config.yaml"]

    # a copy saved earlier (maybe edited since) isn't written over either, and the workflows take no -F
    saved = tmp_path / "saved"
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, *given, "--get-workflow", str(saved))
    (saved / "Snakefile").write_text("edited\n")
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, *given, "--get-workflow", str(saved))
    assert e.value.code == 1
    for flag in ("-F", "--force-overwrite"):
        with pytest.raises(SystemExit) as e:
            run_cli(monkeypatch, *given, "--get-workflow", str(saved), flag)
        assert e.value.code == 2
    assert (saved / "Snakefile").read_text() == "edited\n"


@pytest.mark.parametrize("name", sorted(wf.WORKFLOW_COMMANDS) + sorted(wf.WORKFLOWS))
def test_no_workflow_takes_force_overwrite(name):
    # the snakemake-based commands and workflows never remove or write over an output directory:
    # run again with the same one, they carry on from where they stopped
    parser = run_wf.command_build_parser(name) if name in wf.WORKFLOW_COMMANDS else workflow_parser(name, show_fine=True)
    flags = {flag for action in parser._actions for flag in action.option_strings}
    assert not {"-F", "--force-overwrite"} & flags


@pytest.mark.parametrize("extra", (["--get-workflow", "saved"], ["--workflow-dir", "saved"], ["--", "--quiet"]))
def test_workflows_that_are_commands_run_as_packaged_only(extra, monkeypatch, tmp_path, accs):
    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, "dl-sra-reads", "-a", str(accs), "-o", str(tmp_path / "out"), "--dry-run", *extra)
    assert e.value.code == 2


@patch("bit.modules.workflows.run_snakemake")
def test_run_builds_the_expected_snakemake_command(mock_run, monkeypatch, tmp_path, accs):
    out = tmp_path / "out"
    # (the dry run bit makes first is made for real in test_a_dry_run_says_what_would_be_set_up_and_plans_the_rest)
    with fake_dry_run(""):
        run_cli(monkeypatch, "dl-sra-reads", "-a", str(accs), "-t", "4", "-j", "3", "-o", str(out), "--dry-run")

    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "snakemake"
    assert cmd[cmd.index("--cores") + 1] == "12"             # jobs x threads
    assert cmd[cmd.index("--directory") + 1] == str(out)
    assert cmd[cmd.index("--software-deployment-method") + 1] == "conda"
    assert "--dry-run" in cmd
    assert cmd[cmd.index("--snakefile") + 1].endswith("dl_sra_reads/Snakefile")
    assert yaml.safe_load((out / "config.yaml").read_text())["threads"] == 4


@patch("bit.modules.workflows.run_snakemake")
def test_workflow_dir_uses_the_saved_copy_and_its_config(mock_run, monkeypatch, tmp_path, reads_dir):
    dest = tmp_path / "saved"
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, "metagenomics", "-r", str(reads_dir), "-d", str(tmp_path / "ref-dbs"), "-t", "3",
                "--get-workflow", str(dest))

    out = tmp_path / "out"
    with fake_dry_run(""):
        run_cli(monkeypatch, "metagenomics", "--workflow-dir", str(dest), "-o", str(out), "--dry-run", "--", "--quiet")

    cmd = mock_run.call_args[0][0]
    assert cmd[cmd.index("--snakefile") + 1] == str(dest / "Snakefile")
    assert "--quiet" in cmd                                   # given after `--`, passed to snakemake as is
    assert yaml.safe_load((out / "config.yaml").read_text())["threads"] == 3


def test_config_pairs_go_last_in_the_snakemake_command():
    cmd = build_snakemake_cmd("x.smk", "out", 4, config_pairs=["a=1", "b=2"], rerun_incomplete=True)
    assert cmd[-3:] == ["--config", "a=1", "b=2"]
    assert "--rerun-incomplete" in cmd
    assert "--software-deployment-method" not in cmd


def test_passed_through_targets_are_not_taken_as_option_values():
    cmd = build_snakemake_cmd("x.smk", "out", 4, resources=["mem_mb=8000", "cpus=4"], config_pairs=["a=1"],
                              extra_args=["setup_checkm2_db", "setup_KOFamScan_db", "--executor", "slurm"])
    assert cmd[:5] == ["snakemake", "setup_checkm2_db", "setup_KOFamScan_db", "--executor", "slurm"]
    resources = cmd.index("--resources")
    assert cmd[resources + 1:resources + 3] == ["mem_mb=8000", "cpus=4"]
    assert cmd[resources + 3].startswith("--")


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
@pytest.mark.parametrize("name", list(wf.WORKFLOW_DIRS))
def test_packaged_and_saved_workflows_dry_run(name, tmp_path, accs, reads_dir, genomes_dir, ref_dbs):
    """Snakemake itself parses and plans the workflow, from the package and from a saved copy."""
    workflow = wf.load_workflow(name)
    given = {
        "dl-sra-reads": {"accessions_file": str(accs)},
        "metagenomics": {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs)},
        "amplicon": {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs)},
        "characterize-genomes": {"genomes_dir": str(genomes_dir), "ref_dbs_dir": str(ref_dbs)},
        "classify-seqs": {"fasta_files": [str(genomes_dir)], "ref_dbs_dir": str(ref_dbs)},
        "magrec": {"reads_dir": str(reads_dir), "assemblies": str(genomes_dir), "ref_dbs_dir": str(ref_dbs)},
    }[name]
    if name == "magrec":
        for sample in ("s1", "s2"):
            (genomes_dir / f"{sample}.fasta").write_text(">c1\nACGT\n")
    settings = wf.resolve_settings(workflow.params, given)

    out = tmp_path / "out"
    out.mkdir()
    wf.write_config(workflow, settings, out / "config.yaml")
    cmd = build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True)
    saved = wf.get_workflow(workflow, settings, tmp_path / "saved")
    packaged_run, saved_run = dry_runs(cmd, (["snakemake", "--dry-run", "--cores", "2"], {"cwd": saved}))
    assert packaged_run.returncode == 0, packaged_run.stdout + packaged_run.stderr
    assert saved_run.returncode == 0, saved_run.stdout + saved_run.stderr


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
@pytest.mark.parametrize("euk, expected, absent", [(False, ["gtdbtk_classify"], ["run_eukcc", "tax_label_seqs"]),
                                                   (True, ["run_eukcc", "tax_search_proteins", "tax_label_seqs",
                                                           "combine_euk_tax_outputs"], ["gtdbtk_classify", "CAT"])])
def test_characterize_genomes_plans_the_steps_for_the_kind_of_genome(euk, expected, absent, tmp_path, genomes_dir, ref_dbs,
                                                                    on_linux):
    workflow = wf.load_workflow("characterize-genomes")
    given = {"genomes_dir": str(genomes_dir), "ref_dbs_dir": str(ref_dbs), "euk": euk}
    settings = wf.resolve_settings(workflow.params, given)
    out = tmp_path / "out"
    out.mkdir()
    wf.write_config(workflow, settings, out / "config.yaml")
    cmd = build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True)
    import os
    result = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, **on_linux})
    assert result.returncode == 0, result.stdout + result.stderr
    assert all(step in result.stdout for step in expected) and not any(step in result.stdout for step in absent)


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_output_prefix_gets_a_dash_unless_it_ends_in_one(tmp_path, genomes_dir, ref_dbs, on_linux):
    prefixes = {"plain": ("Swamp", "Swamp-genome-summaries.tsv"), "dash": ("Swamp-", "Swamp-genome-summaries.tsv"),
                "underscore": ("Swamp_", "Swamp_genome-summaries.tsv"), "none": ("", " genome-summaries.tsv")}
    import os
    cmds = []
    for out_name, (prefix, _) in prefixes.items():
        workflow = wf.load_workflow("characterize-genomes")
        settings = wf.resolve_settings(workflow.params, {"genomes_dir": str(genomes_dir), "ref_dbs_dir": str(ref_dbs),
                                                         "output_prefix": prefix})
        out = tmp_path / out_name
        out.mkdir()
        wf.write_config(workflow, settings, out / "config.yaml")
        cmds.append((build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True),
                     {"env": {**os.environ, **on_linux}}))
    for (out_name, (prefix, expected)), result in zip(prefixes.items(), dry_runs(*cmds)):
        assert result.returncode == 0, result.stdout + result.stderr
        assert expected in result.stdout and "Swamp--" not in result.stdout, (prefix, result.stdout)


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
@pytest.mark.parametrize("given, expected, absent", [
    ({}, ["-a ^GTGCCAGCMGCCGCGGTAA...TTAGAWACCCBDGTAGTCC -A ^GGACTACHVGGGTWTCTAA...TTACCGCGGCKGCTGGCAC", "--discard-untrimmed"], []),
    ({"primers_not_linked": True, "keep_untrimmed": True}, ["-g ^GTGCCAGCMGCCGCGGTAA -G ^GGACTACHVGGGTWTCTAA"], ["--discard-untrimmed"]),
    ({"single_end": True}, ["-a ^GTGCCAGCMGCCGCGGTAA...TTAGAWACCCBDGTAGTCC -o", "Illumina-SE-R-processing.R"], [" -A ", " -p trimmed-reads"]),
    ({"skip_primer_trimming": True}, ["Illumina-PE-R-processing.R", '"FALSE" "sample-IDs.txt" "raw-reads/"'], ["cutadapt"]),
])
def test_amplicon_plans_the_expected_commands(given, expected, absent, tmp_path, reads_dir, ref_dbs):
    workflow = wf.load_workflow("amplicon")
    if given.get("single_end"):
        for f in reads_dir.glob("*_R2*"):
            f.unlink()
    settings = wf.resolve_settings(workflow.params, {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs), **given})
    out = tmp_path / "out"
    out.mkdir()
    wf.write_config(workflow, settings, out / "config.yaml")
    cmd = build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True, extra_args=["-p"])
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    for text in expected:
        assert text in result.stdout
    for text in absent:
        assert text not in result.stdout


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_dl_sra_reads_combines_runs_when_sample_names_are_given(tmp_path):
    workflow = wf.load_workflow("dl-sra-reads")

    def dry_run_cmd(text, name):
        accs = tmp_path / f"{name}.tsv"
        accs.write_text(text)
        out = tmp_path / name
        out.mkdir()
        wf.write_config(workflow, wf.resolve_settings(workflow.params, {"accessions_file": str(accs)}), out / "config.yaml")
        return build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True)

    bad = ("S1\tSRR1\nSRR2\n", "S1\tSRR1\nS2\tSRR1\n", "SRR9\tSRR1\n")
    mapped, plain, *bad_runs = dry_runs(dry_run_cmd("S1\tSRR1\nS1\tSRR2\nS2\tERR3\n", "mapped"),
                                        dry_run_cmd("SRR1\nSRR2\n", "plain"),
                                        *(dry_run_cmd(text, f"bad{n}") for n, text in enumerate(bad)))

    assert mapped.returncode == 0, mapped.stdout + mapped.stderr
    counts = {line.split()[0]: line.split()[1] for line in mapped.stdout.splitlines() if len(line.split()) == 2}
    assert counts["combine_runs"] == "2" and counts["prefetch"] == "3"

    assert plain.returncode == 0 and "combine_runs" not in plain.stdout

    assert len(bad_runs) == len(bad)
    for result in bad_runs:
        assert result.returncode != 0


@pytest.mark.parametrize("name, extra, expected", [
    ("metagenomics", [], ["tax-db", "kofamscan-db", "checkm2-db", "gtdbtk-db"]),
    ("metagenomics", ["--skip-binning"], ["tax-db", "kofamscan-db"]),
    ("characterize-genomes", [], ["checkm2-db", "gtdbtk-db"]),
    ("characterize-genomes", ["--euk"], ["eukcc-db", "tax-db"]),
    ("amplicon", [], ["amplicon-refs"]),
    ("classify-seqs", [], ["tax-db"]),
    ("magrec", [], ["checkm2-db", "gtdbtk-db"]),
    ("metagenomics", ["--binning-only"], ["checkm2-db", "gtdbtk-db"]),
])
@patch("bit.modules.workflows.run_snakemake")
@patch("bit.modules.ref_dbs.set_up_missing")
def test_setup_dbs_only_sets_up_what_a_run_would_use_and_stops(set_up, snakemake, name, extra, expected, monkeypatch, tmp_path):
    """No inputs are given; the databases are set up by bit itself, and snakemake isn't run."""
    monkeypatch.setattr("bit.modules.ref_dbs.GTDBTK_RUNS_HERE", True)
    run_cli(monkeypatch, name, "-d", str(tmp_path / "dbs"), "--setup-dbs-only", "-o", str(tmp_path / "out"), *extra)
    root, names = set_up.call_args[0]
    assert root == tmp_path / "dbs" and names == expected
    snakemake.assert_not_called()
    assert not (tmp_path / "out").exists()


def test_setup_dbs_only_is_not_offered_without_reference_databases():
    assert "--setup-dbs-only" not in workflow_parser("dl-sra-reads", show_fine=True).format_help()
    assert "--setup-dbs-only" in workflow_parser("metagenomics").format_help()


@patch("bit.modules.ref_dbs.set_up_missing")
def test_the_reference_directory_comes_from_the_variable_unless_given(set_up, monkeypatch, tmp_path):
    monkeypatch.setenv("BIT_REF_DBS_DIR", str(tmp_path / "from-variable"))
    run_cli(monkeypatch, "classify-seqs", "--setup-dbs-only")
    assert set_up.call_args[0][0] == tmp_path / "from-variable"
    run_cli(monkeypatch, "classify-seqs", "--setup-dbs-only", "-d", str(tmp_path / "given"))
    assert set_up.call_args[0][0] == tmp_path / "given"


@patch("bit.modules.ref_dbs.set_up_missing")
def test_without_a_reference_directory_the_run_stops_except_amplicons(set_up, monkeypatch, tmp_path, capsys):
    monkeypatch.delenv("BIT_REF_DBS_DIR", raising=False)
    with pytest.raises(SystemExit):
        run_cli(monkeypatch, "classify-seqs", "--setup-dbs-only")
    assert "bit data locations set" in " ".join(capsys.readouterr().out.split())
    set_up.assert_not_called()

    # the amplicon references are small: they go in the output directory instead
    run_cli(monkeypatch, "amplicon", "--setup-dbs-only", "-o", str(tmp_path / "amp"))
    assert set_up.call_args[0][0] == tmp_path / "amp" / "ref-dbs"


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
@patch("bit.modules.ref_dbs.set_up_missing")
def test_a_dry_run_says_what_would_be_set_up_and_plans_the_rest(set_up, monkeypatch, tmp_path, genomes_dir, capsys):
    run_cli(monkeypatch, "classify-seqs", "-f", str(genomes_dir), "-d", str(tmp_path / "nothing-here"),
            "-o", str(tmp_path / "out"), "--dry-run")
    set_up.assert_not_called()
    assert "would be set up first" in capsys.readouterr().out


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_a_saved_copy_run_directly_names_the_missing_databases(tmp_path, genomes_dir, ref_dbs):
    workflow = wf.load_workflow("classify-seqs")
    (ref_dbs / "bit-tax-db-2026-11" / "proteins.dmnd").unlink()
    settings = wf.resolve_settings(workflow.params, {"fasta_files": [str(genomes_dir)], "ref_dbs_dir": str(ref_dbs)})
    saved = wf.get_workflow(workflow, settings, tmp_path / "saved")
    result = subprocess.run(["snakemake", "--dry-run", "--cores", "2"], cwd=saved, capture_output=True, text=True)
    assert result.returncode != 0
    assert f"bit data get tax-db -d {ref_dbs}/" in result.stdout + result.stderr


def planned(workflow_name, given, tmp_path, out_name="out", env=None):
    """ the steps a dry run of the packaged workflow plans: {rule: number of jobs} (or the failed run) """
    return planned_each([(workflow_name, given, out_name)], tmp_path, env=env)[0]


def planned_each(runs, tmp_path, env=None):
    """ planned() for each (workflow name, settings given, output directory name), the dry runs made at the same time """
    import os
    cmds = []
    for workflow_name, given, out_name in runs:
        workflow = wf.load_workflow(workflow_name)
        settings = wf.resolve_settings(workflow.params, given)
        out = tmp_path / out_name
        out.mkdir()
        wf.write_config(workflow, settings, out / "config.yaml")
        cmds.append(build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True))
    results = dry_runs(*((cmd, {"env": {**os.environ, **(env or {})}}) for cmd in cmds))
    return [result if result.returncode else
            {line.split()[0]: int(line.split()[1]) for line in result.stdout.splitlines()
             if len(line.split()) == 2 and line.split()[1].isdigit()}
            for result in results]


def _platform_env(tmp_path, system):
    """ an environment for a subprocess in which platform.system() says `system` (what GTDB-Tk's use depends on) """
    import os
    d = tmp_path / f"platform-{system}"
    d.mkdir()
    (d / "sitecustomize.py").write_text(f"import platform\nplatform.system = lambda: {system!r}\n")
    return {"PYTHONPATH": os.pathsep.join([str(d)] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p])}


@pytest.fixture
def on_macos(tmp_path):
    return _platform_env(tmp_path, "Darwin")


@pytest.fixture
def on_linux(tmp_path):
    """ for the tests of the steps planned where GTDB-Tk runs, whatever system the tests themselves run on """
    return _platform_env(tmp_path, "Linux")


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_magrec_bins_each_samples_assembly_with_its_own_reads(tmp_path, reads_dir, assemblies_dir, ref_dbs, on_linux):
    jobs = planned("magrec", {"reads_dir": str(reads_dir), "assemblies": str(assemblies_dir), "ref_dbs_dir": str(ref_dbs)},
                   tmp_path, env=on_linux)
    assert isinstance(jobs, dict), jobs.stdout + jobs.stderr
    assert {rule: jobs[rule] for rule in ("copy_assembly", "index_assembly", "mapping", "metabat_binning")} == \
           {"copy_assembly": 2, "index_assembly": 2, "mapping": 2, "metabat_binning": 2}
    assert jobs["gtdbtk_on_MAGs"] == 1 and "assemble" not in jobs and "bbduk" not in jobs


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_magrec_bins_a_co_assembly_with_every_samples_reads(tmp_path, reads_dir, ref_dbs):
    co = tmp_path / "co-assembly.fa"
    co.write_text(">k141_1\nACGT\n")
    jobs = planned("magrec", {"reads_dir": str(reads_dir), "assemblies": str(co), "ref_dbs_dir": str(ref_dbs)}, tmp_path)
    assert isinstance(jobs, dict), jobs.stdout + jobs.stderr
    assert {rule: jobs[rule] for rule in ("copy_assembly", "index_assembly", "mapping", "metabat_binning")} == \
           {"copy_assembly": 1, "index_assembly": 1, "mapping": 2, "metabat_binning": 1}


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_magrec_says_which_samples_have_no_assembly(tmp_path, reads_dir, assemblies_dir, ref_dbs):
    (assemblies_dir / "s2.fasta").unlink()
    result = planned("magrec", {"reads_dir": str(reads_dir), "assemblies": str(assemblies_dir), "ref_dbs_dir": str(ref_dbs)}, tmp_path)
    assert not isinstance(result, dict)
    assert "for these samples: s2" in " ".join((result.stdout + result.stderr).split())


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_metagenomics_binning_only_skips_annotation_and_taxonomy(tmp_path, reads_dir, ref_dbs, on_linux):
    given = {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs)}
    full, binning, both = planned_each([("metagenomics", given, "full"),
                                        ("metagenomics", {**given, "binning_only": True}, "binning"),
                                        ("metagenomics", {**given, "binning_only": True, "skip_binning": True}, "both")],
                                       tmp_path, env=on_linux)
    assert isinstance(binning, dict), binning.stdout + binning.stderr
    skipped = {"call_genes", "KO_annotation", "tax_search_proteins", "tax_label_seqs", "make_combined_gene_level_tables",
               "summarize_MAG_KO_annots_with_KEGG_Decoder"}
    assert skipped <= set(full) and not skipped & set(binning)
    assert {"assemble", "bbduk", "mapping", "metabat_binning", "gtdbtk_on_MAGs", "dereplicate_MAGs"} <= set(binning)

    assert not isinstance(both, dict) and "can't both be given" in both.stdout + both.stderr


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_on_macos_genome_taxonomy_comes_from_the_tax_db(tmp_path, reads_dir, genomes_dir, ref_dbs, on_macos):
    mg, genomes = planned_each(
        [("metagenomics", {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs), "binning_only": True}, "mg"),
         ("characterize-genomes", {"genomes_dir": str(genomes_dir), "ref_dbs_dir": str(ref_dbs)}, "characterize")],
        tmp_path, env=on_macos)
    assert isinstance(mg, dict), mg.stdout + mg.stderr
    assert {"tax_db_on_MAGs_search", "tax_db_on_MAGs_label"} <= set(mg) and "gtdbtk_on_MAGs" not in mg

    assert isinstance(genomes, dict), genomes.stdout + genomes.stderr
    assert {"tax_db_on_genomes_search", "tax_db_on_genomes_label", "run_checkm2"} <= set(genomes)
    assert "gtdbtk_classify" not in genomes


@patch("bit.modules.ref_dbs.GTDBTK_RUNS_HERE", False)
def test_on_macos_the_tax_db_is_set_up_instead_of_gtdbtks():
    for name, settings, expected in (
            ("metagenomics", {"skip_binning": False, "binning_only": True}, ["checkm2-db", "tax-db"]),
            ("metagenomics", {"skip_binning": False, "binning_only": False}, ["tax-db", "kofamscan-db", "checkm2-db"]),
            ("magrec", {}, ["checkm2-db", "tax-db"]),
            ("characterize-genomes", {"euk": False}, ["checkm2-db", "tax-db"])):
        assert wf.load_workflow(name).ref_dbs_needed(settings) == expected


@pytest.mark.parametrize("name, extra, said", [
    ("magrec", [], True), ("metagenomics", [], True), ("metagenomics", ["--skip-binning"], False),
    ("characterize-genomes", [], True), ("characterize-genomes", ["--euk"], False), ("classify-seqs", [], False),
])
@patch("bit.modules.ref_dbs.set_up_missing")
def test_on_macos_a_run_that_would_use_gtdbtk_says_why_it_doesnt(set_up, name, extra, said, monkeypatch, tmp_path, capsys):
    monkeypatch.setattr("bit.modules.ref_dbs.GTDBTK_RUNS_HERE", False)
    run_cli(monkeypatch, name, "-d", str(tmp_path / "dbs"), "--setup-dbs-only", *extra)
    out = " ".join(capsys.readouterr().out.split())
    assert ("GTDB-Tk can't be run on macOS" in out and "best way" in out) == said


def test_an_option_value_starting_with_a_dash_is_taken_as_its_value(monkeypatch, tmp_path, reads_dir, assemblies_dir):
    parser = workflow_parser("magrec")
    assert run_wf.attach_dash_values(parser, ["-x", "-assembly.fasta", "-r", "r", "-t", "-h", "--", "-x", "-y"]) == \
           ["-x=-assembly.fasta", "-r", "r", "-t", "-h", "--", "-x", "-y"]     # (-h is an option: left alone)

    with pytest.raises(SystemExit) as e:
        run_cli(monkeypatch, "magrec", "-r", str(reads_dir), "-a", str(assemblies_dir), "-x", "-assembly.fasta",
                "--get-workflow", str(tmp_path / "saved"))
    assert e.value.code == 0
    assert yaml.safe_load((tmp_path / "saved" / "config.yaml").read_text())["assembly_extension"] == "-assembly.fasta"


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_classify_seqs_genomes_mode_gives_one_table(tmp_path, genomes_dir, ref_dbs):
    given = {"fasta_files": [str(genomes_dir)], "ref_dbs_dir": str(ref_dbs)}
    per_seq, genomes = planned_each([("classify-seqs", given, "per-seq"),
                                     ("classify-seqs", {**given, "genomes": True}, "whole")], tmp_path)
    assert isinstance(genomes, dict), genomes.stdout + genomes.stderr
    assert genomes["combine_genome_tables"] == 1 and "combine_genome_tables" not in per_seq
    assert {"tax_kmer_step", "tax_search_proteins", "tax_label_seqs"} <= set(genomes)
    # genes are only called after the k-mers, for the genomes they didn't label
    assert "tax_genes_of_undecided" in genomes and "call_genes" not in genomes
    assert "call_genes" in per_seq and "tax_genes_of_undecided" not in per_seq


DRY_RUN_OUTPUT = """
rule assemble:
    input: reads/s1_R1.fastq.gz
    resources: tmpdir=/tmp, cpus=8, mem_mb=100000, mem=100 GB, mem_mib=95368
rule assemble:
    resources: tmpdir=/tmp, cpus=8, mem_mb=100000, mem=100 GB, mem_mib=95368
localrule setup_gtdbtk_db:
    resources: tmpdir=/tmp
rule gtdbtk_on_MAGs:
    resources: tmpdir=/tmp, cpus=8, mem_mb=175000, mem=175 GB, mem_mib=166893
rule all:
    resources: tmpdir=/tmp
"""


def fake_dry_run(output=DRY_RUN_OUTPUT, returncode=0):
    return patch("bit.modules.workflows.subprocess.run",
                 return_value=subprocess.CompletedProcess([], returncode, stdout=output, stderr=""))


def test_planned_jobs_memory_reads_each_steps_request():
    with fake_dry_run():
        steps = wf.planned_jobs_memory(wf.dry_run_report(["snakemake", "-n"]))
    assert steps == {"assemble": 100000, "setup_gtdbtk_db": None, "gtdbtk_on_MAGs": 175000, "all": None}
    with fake_dry_run(returncode=1):
        assert wf.planned_jobs_memory(wf.dry_run_report(["snakemake", "-n"])) is None


def test_a_step_needing_more_memory_than_there_is_stops_the_run(capsys):
    with fake_dry_run(), pytest.raises(SystemExit):
        wf.check_steps_can_run(wf.dry_run_report(["snakemake", "-n"]), 64000, memory_given=False)
    err = capsys.readouterr()
    message = " ".join((err.out + err.err).split())
    assert "'assemble' (100 GB)" in message and "'gtdbtk_on_MAGs' (175 GB)" in message
    assert "this machine has (64 GB)" in message

    with fake_dry_run():
        wf.check_steps_can_run(wf.dry_run_report(["snakemake", "-n"]), 245000, memory_given=False)     # everything fits: no exit


@patch("bit.modules.workflows.run_snakemake")
def test_memory_total_goes_to_snakemake_unless_jobs_go_to_a_cluster(mock_run, monkeypatch, tmp_path, reads_dir, ref_dbs):
    with fake_dry_run(""):
        run_cli(monkeypatch, "metagenomics", "-r", str(reads_dir), "-d", str(ref_dbs), "-o", str(tmp_path / "a"),
                "--memory", "48", "--dry-run")
    cmd = mock_run.call_args[0][0]
    assert cmd[cmd.index("--resources") + 1] == "mem_mb=48000"
    assert "memory_total_mb=48000" in cmd[cmd.index("--config"):]

    # (the dry run made first is still made, for the conda environments to create, but its steps' memory isn't
    # held against this machine's: here one asks for more than any machine has)
    with fake_dry_run("rule assemble:\n    resources: tmpdir=/tmp, mem_mb=99000000\n"):
        run_cli(monkeypatch, "metagenomics", "-r", str(reads_dir), "-d", str(ref_dbs), "-o", str(tmp_path / "b"),
                "--dry-run", "--", "--executor", "slurm")
    cmd = mock_run.call_args[0][0]
    assert "--resources" not in cmd
    assert "memory_total_mb" not in " ".join(cmd)


@pytest.mark.full
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not available")
def test_megahit_memory_is_kept_within_the_runs_total(tmp_path, reads_dir, ref_dbs):
    """megahit's 32 GB default is a ceiling, so on a 16 GB machine it's given 16 GB rather than stopping the run."""
    workflow = wf.load_workflow("metagenomics")
    settings = wf.resolve_settings(workflow.params, {"reads_dir": str(reads_dir), "ref_dbs_dir": str(ref_dbs)})

    def assemble_mem_mb(config_pairs):
        out = tmp_path / str(len(config_pairs))
        out.mkdir()
        wf.write_config(workflow, settings, out / "config.yaml")
        cmd = build_snakemake_cmd(workflow.directory / "Snakefile", out, 2, dry_run=True, config_pairs=config_pairs)
        return wf.planned_jobs_memory(wf.dry_run_report(cmd))["assemble"]

    # (the two dry runs are made at the same time)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(2) as pool:
        as_set, within_16_gb = pool.map(assemble_mem_mb, ([], ["memory_total_mb=16000"]))
    assert as_set == 32000
    assert within_16_gb == 16000
