import pytest
import pandas as pd
from argparse import Namespace
from unittest.mock import patch, MagicMock

from bit.modules.ncbi.get_accessions_from_ncbi import (
    FIELD_MAP,
    build_summary_command,
    run_ncbi_summary,
    get_accessions,
    get_taxon_count,
    get_accessions_from_ncbi,
)


# ─── helpers / fixtures ───────────────────────────────────────────────────────

def make_args(**kwargs):
    defaults = {
        'target_taxon': None,
        'source': 'both',
        'reference_genomes_only': False,
        'assembly_level': None,
        'annotated_only': False,
        'get_taxon_counts': False,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


# the human-readable header row that `dataformat tsv genome` emits, in the same order
# as FIELD_MAP (we re-key by position in the module, so only the order matters here)
DATAFORMAT_HEADER = "\t".join([
    "Assembly Accession", "Organism Name", "Organism Taxonomic ID",
    "Assembly Name", "Assembly Level", "Assembly RefSeq Category",
    "Assembly Stats Total Sequence Length", "CheckM Completeness",
    "CheckM Contamination", "Source Database",
])

# rows chosen to exercise: GenBank + RefSeq source values, a blank refseq_category,
# blank/whitespace checkm values, and an "na" placeholder (all should normalize to "NA").
# column order matches FIELD_MAP: accession, organism_name, tax_id, assembly_name,
# assembly_level, refseq_category, total_sequence_length, checkm_completeness,
# checkm_contamination, source_database
DATAFORMAT_ROWS = [
    "GCA_000196795.1\tAlteromonas macleodii\t28108\tASM19679v1\tComplete Genome\tna\t4467846\t99.5\t0.3\tSOURCE_DATABASE_GENBANK",
    "GCF_000196795.1\tAlteromonas macleodii\t28108\tASM19679v1\tComplete Genome\treference genome\t4467846\t99.5\t0.3\tSOURCE_DATABASE_REFSEQ",
    "GCA_111111111.1\tAlteromonas sp.\t28108\tASM1v1\tScaffold\t\t3500000\t  \t\tSOURCE_DATABASE_GENBANK",
]


def make_dataformat_tsv(rows=None):
    """Build a full dataformat-style TSV string (header + given rows)."""
    rows = DATAFORMAT_ROWS if rows is None else rows
    return DATAFORMAT_HEADER + "\n" + "\n".join(rows) + "\n"


def make_popen_mock(stdout=b"", stderr=b"", returncode=0):
    """A MagicMock standing in for a subprocess.Popen instance."""
    proc = MagicMock()
    proc.communicate.return_value = (stdout, stderr)
    proc.stderr.read.return_value = stderr
    proc.returncode = returncode
    return proc


def patch_pipeline(summary_proc, dataformat_proc):
    """Patch subprocess.Popen so the two calls return the given mocks in order.

    run_ncbi_summary calls Popen twice: first for `datasets summary`, then for
    `dataformat`. side_effect returns them in that call order.
    """
    return patch(
        "bit.modules.ncbi.get_accessions_from_ncbi.subprocess.Popen",
        side_effect=[summary_proc, dataformat_proc],
    )


# ─── build_summary_command ────────────────────────────────────────────────────

def test_build_summary_command_minimal():
    cmd = build_summary_command(make_args(target_taxon="Alteromonas"))
    assert cmd == ["datasets", "summary", "genome", "taxon", "Alteromonas", "--as-json-lines"]


def test_build_summary_command_taxon_with_space_is_single_arg():
    cmd = build_summary_command(make_args(target_taxon="Escherichia coli"))
    # the taxon must be a single argv element, not split on the space
    assert "Escherichia coli" in cmd


def test_build_summary_command_assembly_source():
    cmd = build_summary_command(make_args(target_taxon="Alteromonas", source="refseq"))
    assert "--assembly-source" in cmd
    assert cmd[cmd.index("--assembly-source") + 1] == "RefSeq"


def test_build_summary_command_assembly_source_all_omitted():
    cmd = build_summary_command(make_args(target_taxon="Alteromonas", source="both"))
    assert "--assembly-source" not in cmd


def test_build_summary_command_reference_genomes_only():
    cmd = build_summary_command(make_args(target_taxon="Alteromonas", reference_genomes_only=True))
    assert "--reference" in cmd


def test_build_summary_command_assembly_level_joined_with_commas():
    # CLI accepts space-delimited (a list); datasets wants comma-delimited
    cmd = build_summary_command(make_args(
        target_taxon="Alteromonas", assembly_level=["complete", "chromosome"]))
    assert "--assembly-level" in cmd
    assert cmd[cmd.index("--assembly-level") + 1] == "complete,chromosome"


def test_build_summary_command_annotated_only():
    cmd = build_summary_command(make_args(target_taxon="Alteromonas", annotated_only=True))
    assert "--annotated" in cmd


def test_build_summary_command_all_options_combined():
    cmd = build_summary_command(make_args(
        target_taxon="Alteromonas", source="genbank", reference_genomes_only=True,
        assembly_level=["complete"], annotated_only=True))
    for expected in ("--assembly-source", "--reference", "--assembly-level",
                     "--annotated"):
        assert expected in cmd


# ─── run_ncbi_summary (parsing / normalization) ───────────────────────────────

def test_run_ncbi_summary_parses_and_renames_columns():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv().encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert list(tab.columns) == list(FIELD_MAP.values())
    assert len(tab) == 3


def test_run_ncbi_summary_accessions_intact():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv().encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert tab["accession"].tolist() == [
        "GCA_000196795.1", "GCF_000196795.1", "GCA_111111111.1"]


def test_run_ncbi_summary_blank_and_na_become_NA():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv().encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    # row 0 had refseq_category "na"; row 2 had blank/whitespace cells
    assert tab.loc[0, "refseq_category"] == "NA"
    assert tab.loc[2, "refseq_category"] == "NA"
    assert tab.loc[2, "checkm_completeness"] == "NA"
    assert tab.loc[2, "checkm_contamination"] == "NA"
    # no empty strings should remain anywhere
    assert not (tab == "").any().any()


def test_run_ncbi_summary_real_values_preserved():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv().encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert tab.loc[1, "refseq_category"] == "reference genome"
    assert tab.loc[0, "organism_name"] == "Alteromonas macleodii"


def test_run_ncbi_summary_source_database_cleanup():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv().encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert tab["source_database"].tolist() == ["genbank", "refseq", "genbank"]


def test_run_ncbi_summary_unknown_source_database_passes_through():
    rows = ["GCA_999.1\tx sp.\t1\tASM9\tContig\tna\t100\t90\t1\tSOMETHING_ELSE"]
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=make_dataformat_tsv(rows).encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="x"))
    assert tab.loc[0, "source_database"] == "SOMETHING_ELSE"


def test_run_ncbi_summary_empty_stdout_returns_empty_df():
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=b"", returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        tab = run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert tab.empty
    assert list(tab.columns) == list(FIELD_MAP.values())


def test_run_ncbi_summary_datasets_real_error_exits():
    # datasets exits non-zero with an unrecognized error and no stdout
    summary_proc = make_popen_mock(stderr=b"Error: connection refused", returncode=1)
    dataformat_proc = make_popen_mock(stdout=b"", returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        with pytest.raises(SystemExit):
            run_ncbi_summary(make_args(target_taxon="Alteromonas"))


def test_run_ncbi_summary_dataformat_error_exits(capsys):
    # dataformat exits non-zero (e.g. bad field name)
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(
        stderr=b"Error: field(s) [source-database] not recognized.", returncode=1)
    with patch_pipeline(summary_proc, dataformat_proc):
        with pytest.raises(SystemExit):
            run_ncbi_summary(make_args(target_taxon="Alteromonas"))
    assert "dataformat" in capsys.readouterr().out


def test_run_ncbi_summary_column_count_drift_exits():
    # header with the wrong number of columns triggers the drift guard
    bad = "ColA\tColB\nval1\tval2\n"
    summary_proc = make_popen_mock(returncode=0)
    dataformat_proc = make_popen_mock(stdout=bad.encode(), returncode=0)
    with patch_pipeline(summary_proc, dataformat_proc):
        with pytest.raises(SystemExit):
            run_ncbi_summary(make_args(target_taxon="Alteromonas"))


# ─── get_accessions ───────────────────────────────────────────────────────────

def test_get_accessions_writes_both_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({
        "accession": ["GCA_1.1", "GCF_2.1"],
        "organism_name": ["x", "y"],
    })
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_accessions(make_args(target_taxon="Alteromonas"))
    assert (tmp_path / "NCBI-Alteromonas-accessions.txt").exists()
    assert (tmp_path / "NCBI-Alteromonas-metadata.tsv").exists()


def test_get_accessions_accs_file_contents(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({"accession": ["GCA_1.1", "GCF_2.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_accessions(make_args(target_taxon="Alteromonas"))
    accs = (tmp_path / "NCBI-Alteromonas-accessions.txt").read_text().splitlines()
    assert accs == ["GCA_1.1", "GCF_2.1"]


def test_get_accessions_taxon_with_space_in_filename(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({"accession": ["GCA_1.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_accessions(make_args(target_taxon="Escherichia coli"))
    assert (tmp_path / "NCBI-Escherichia-coli-accessions.txt").exists()


def test_get_accessions_reference_genomes_only_suffix(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({"accession": ["GCA_1.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_accessions(make_args(target_taxon="Alteromonas", reference_genomes_only=True))
    assert (tmp_path / "NCBI-Alteromonas-refseq-ref-accessions.txt").exists()


def test_get_accessions_assembly_source_suffix(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({"accession": ["GCF_1.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_accessions(make_args(target_taxon="Alteromonas", source="refseq"))
    assert (tmp_path / "NCBI-Alteromonas-refseq-accessions.txt").exists()


def test_get_accessions_empty_result_exits_without_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    empty = pd.DataFrame(columns=list(FIELD_MAP.values()))
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=empty):
        with pytest.raises(SystemExit):
            get_accessions(make_args(target_taxon="Alteromonas"))
    assert not (tmp_path / "NCBI-Alteromonas-accessions.txt").exists()


# ─── get_taxon_count ──────────────────────────────────────────────────────────

def test_get_taxon_count_reports_number(capsys):
    tab = pd.DataFrame({"accession": ["GCA_1.1", "GCF_2.1", "GCA_3.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        get_taxon_count(make_args(target_taxon="Alteromonas"))
    assert "3" in capsys.readouterr().out


def test_get_taxon_count_zero(capsys):
    empty = pd.DataFrame(columns=list(FIELD_MAP.values()))
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=empty):
        get_taxon_count(make_args(target_taxon="Alteromonas"))
    assert "0" in capsys.readouterr().out


# ─── get_accessions_from_ncbi (orchestrator) ──────────────────────────────────

def test_orchestrator_taxon_counts_path(capsys):
    tab = pd.DataFrame({"accession": ["GCA_1.1", "GCF_2.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        with pytest.raises(SystemExit):
            get_accessions_from_ncbi(make_args(
                target_taxon="Alteromonas", get_taxon_counts=True))
    assert "2" in capsys.readouterr().out


def test_orchestrator_get_accessions_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tab = pd.DataFrame({"accession": ["GCA_1.1"]})
    with patch("bit.modules.ncbi.get_accessions_from_ncbi.run_ncbi_summary",
               return_value=tab):
        with pytest.raises(SystemExit):
            get_accessions_from_ncbi(make_args(target_taxon="Alteromonas"))
    assert (tmp_path / "NCBI-Alteromonas-accessions.txt").exists()
