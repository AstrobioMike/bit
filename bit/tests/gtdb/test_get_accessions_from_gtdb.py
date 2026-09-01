import pytest # type: ignore
import pandas as pd # type: ignore
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
from argparse import Namespace
from unittest.mock import patch
from bit.modules.gtdb.get_accessions_from_gtdb import (
    report_gtdb_version_info,
    copy_gtdb_table,
    report_taxon_counts,
    report_rank_counts_for_taxon,
    report_unique_taxa_counts_of_all_ranks,
    get_accessions_from_gtdb,
    _resolve_or_exit
)
from bit.modules.gtdb.build_gtdb_data_parquet import PARQUET_FILENAME, VERSION_FILENAME


# ─── helpers / fixtures ───────────────────────────────────────────────────────

def make_args(**kwargs):
    defaults = {
        'get_table': False,
        'get_taxon_counts': False,
        'target_taxon': None,
        'gtdb_representatives_only': False,
        'refseq_reference_genomes_only': False,
        'get_rank_counts': False,
        'target_rank': None,
        'derep_rank': 'off',
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


ROWS = [
    {
        "accession": "RS_GCF_000002125.1",
        "domain": "Archaea", "phylum": "Halobacteriota", "class": "Halobacteria",
        "order": "Halobacteriales", "family": "Haloarculaceae", "genus": "Haloarcula",
        "species": "Haloarcula hispanica",
        "gtdb_representative": "t", "ncbi_refseq_category": "reference genome",
        "ncbi_genbank_assembly_accession": "GCA_000002125.1",
        "checkm2_completeness": "99.0", "checkm2_contamination": "0.5",
        "genome_size": "4000000", "contig_count": "1",
    },
    {
        "accession": "RS_GCF_000009905.1",
        "domain": "Archaea", "phylum": "Halobacteriota", "class": "Halobacteria",
        "order": "Halobacteriales", "family": "Haloarculaceae", "genus": "Haloarcula",
        "species": "Haloarcula marismortui",
        "gtdb_representative": "f", "ncbi_refseq_category": "na",
        "ncbi_genbank_assembly_accession": "GCA_000009905.1",
        "checkm2_completeness": "85.0", "checkm2_contamination": "1.0",
        "genome_size": "4000000", "contig_count": "5",
    },
    {
        "accession": "GB_GCA_000001405.1",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Escherichia",
        "species": "Escherichia coli",
        "gtdb_representative": "t", "ncbi_refseq_category": "reference genome",
        "ncbi_genbank_assembly_accession": "GCA_000001405.1",
        "checkm2_completeness": "99.5", "checkm2_contamination": "0.2",
        "genome_size": "5000000", "contig_count": "1",
    },
    {
        "accession": "GB_GCA_000005845.2",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Escherichia",
        "species": "Escherichia coli",
        "gtdb_representative": "f", "ncbi_refseq_category": "na",
        "ncbi_genbank_assembly_accession": "GCA_000005845.2",
        "checkm2_completeness": "90.0", "checkm2_contamination": "0.8",
        "genome_size": "5000000", "contig_count": "3",
    },
    {
        "accession": "GB_GCA_000006925.2",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Salmonella",
        "species": "Salmonella enterica",
        "gtdb_representative": "t", "ncbi_refseq_category": "reference genome",
        "ncbi_genbank_assembly_accession": "GCA_000006925.2",
        "checkm2_completeness": "98.0", "checkm2_contamination": "0.4",
        "genome_size": "4800000", "contig_count": "2",
    },
]


# Every GENBANK accession in ROWS. GTDB is a snapshot, so its entries are screened
# against the NCBI asset for liveness before anything is written; these tests stub that
# asset out so the screen is a no-op unless a test deliberately suppresses something.
ALL_GENBANK_ACCS = [r["ncbi_genbank_assembly_accession"] for r in ROWS]


def _write_ncbi_screen_asset(path, accessions):
    """A minimal NCBI asset carrying just the column the liveness screen reads."""
    pq.write_table(
        pa.table({"assembly_accession": pa.array(list(accessions), type=pa.string())}),
        str(path))
    return str(path)


@pytest.fixture(autouse=True)
def live_ncbi_screen(tmp_path_factory):
    """
    Stub the NCBI liveness screen so it finds every GTDB accession alive.

    Without this, get-accs-from-gtdb would try to FETCH the NCBI asset (the screen
    always-fetches by design), which needs a set env var and network. Autouse so the
    existing expectations -- written before the screen existed -- keep holding.
    """
    path = tmp_path_factory.mktemp("ncbi-screen") / "ncbi-data.parquet"
    _write_ncbi_screen_asset(path, ALL_GENBANK_ACCS)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.ncbi_screen_path",
               return_value=str(path)):
        yield str(path)


@pytest.fixture
def gtdb_tab():
    return pd.DataFrame(ROWS)


@pytest.fixture
def gtdb_dir(tmp_path):
    """A tmp directory holding just the GTDB version-info file (for the version reader)."""
    (tmp_path / VERSION_FILENAME).write_text("R220\n2024-04-24\n")
    return tmp_path

@pytest.fixture
def gtdb_parquet(tmp_path):
    """
    The real R4 contract: get_gtdb_data returns the .parquet PATH, with the
    version-info file sitting beside it. Orchestrator tests use this.
    """
    (tmp_path / VERSION_FILENAME).write_text("R220\n2024-04-24\n")
    df = pd.DataFrame(ROWS)
    path = tmp_path / PARQUET_FILENAME
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), str(path))
    return str(path)


# ─── report_gtdb_version_info ─────────────────────────────────────────────────

def test_report_gtdb_version_info_prints_version(gtdb_dir, capsys):
    report_gtdb_version_info(str(gtdb_dir))
    out = capsys.readouterr().out
    assert "R220" in out


# ─── copy_gtdb_table ──────────────────────────────────────────────────────────

def test_copy_gtdb_table_copies_file(gtdb_parquet, tmp_path, monkeypatch):
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    monkeypatch.chdir(out_dir)
    copy_gtdb_table(gtdb_parquet)
    assert (out_dir / "gtdb-arc-and-bac-metadata.tsv").exists()


# ─── selection + writing (end to end, through the orchestrator) ───────────────

def _run(gtdb_parquet, **kw):
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(**kw))


def test_taxon_writes_accs_and_metadata(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Escherichia")
    accs = (tmp_path / "gtdb-escherichia-genus-accs.txt").read_text().splitlines()
    assert len(accs) == 2
    assert (tmp_path / "gtdb-escherichia-genus-metadata.tsv").exists()


def test_metadata_tsv_keeps_every_column_of_the_asset(gtdb_parquet, tmp_path, monkeypatch):
    """The selection path carries the full asset schema, not a chosen subset."""
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Escherichia")
    header = (tmp_path / "gtdb-escherichia-genus-metadata.tsv").read_text().splitlines()[0]
    assert header.split("\t") == list(pq.ParquetFile(gtdb_parquet).schema_arrow.names)


def test_metadata_tsv_is_plain_unquoted_tsv(gtdb_parquet, tmp_path, monkeypatch):
    """
    Written through Arrow now rather than pandas. Arrow quotes column names, and its
    "needed" quoting style quotes every string value, so a naive port would put quotes
    around every field -- this pins the plain output.
    """
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Escherichia")
    text = (tmp_path / "gtdb-escherichia-genus-metadata.tsv").read_text()
    assert '"' not in text
    assert text.splitlines()[1].count("\t") == text.splitlines()[0].count("\t")


def test_species_with_a_space_is_dashed_in_the_filename(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Escherichia coli", target_rank="species")
    assert (tmp_path / "gtdb-escherichia-coli-species-accs.txt").exists()


def test_representatives_only_filters_and_names_the_file(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Escherichia", gtdb_representatives_only=True)
    accs = (tmp_path / "gtdb-escherichia-genus-gtdb-rep-accs.txt").read_text().splitlines()
    assert len(accs) == 1


def test_derep_keeps_one_genome_per_group(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="Enterobacteriaceae", derep_rank="genus")
    accs = (tmp_path / "gtdb-enterobacteriaceae-family-accs.txt").read_text().splitlines()
    assert len(accs) == 2      # Escherichia + Salmonella


def test_all_writes_accessions_and_no_metadata(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all")
    accs = (tmp_path / "gtdb-arc-and-bac-accessions.txt").read_text().splitlines()
    assert len(accs) == 5
    assert not list(tmp_path.glob("*-rep-metadata.tsv"))


def test_all_with_gtdb_reps_names_the_source_actually_used(gtdb_parquet, tmp_path, monkeypatch):
    """
    The 'all' branch hardcoded "refseq-rep" in its filenames, so
    `-t all --gtdb-representatives-only` wrote GTDB representatives into a file named
    for RefSeq reference genomes.
    """
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all", gtdb_representatives_only=True)
    accs = (tmp_path / "gtdb-arc-and-bac-gtdb-rep-accessions.txt").read_text().splitlines()
    assert len(accs) == 3
    assert (tmp_path / "gtdb-arc-and-bac-gtdb-rep-metadata.tsv").exists()
    assert not (tmp_path / "gtdb-arc-and-bac-refseq-rep-accessions.txt").exists()


def test_all_with_refseq_reps_names_refseq(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all", refseq_reference_genomes_only=True)
    accs = (tmp_path / "gtdb-arc-and-bac-refseq-rep-accessions.txt").read_text().splitlines()
    assert len(accs) == 3


def test_taxon_not_found_exits_without_files(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="NonExistent")
    assert not list(tmp_path.glob("gtdb-*accs.txt"))


# ─── report_taxon_counts ──────────────────────────────────────────────────────

def test_report_taxon_counts_all_total(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "all")
    assert "5" in capsys.readouterr().out


def test_report_taxon_counts_all_with_reps(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "all", representatives_source="gtdb")
    out = capsys.readouterr().out
    assert "5" in out   # total
    assert "3" in out   # gtdb_representative == "t"


def test_report_taxon_counts_specific_taxon(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "Escherichia")
    assert "2" in capsys.readouterr().out


def test_report_taxon_counts_not_found_exits(gtdb_parquet):
    with pytest.raises(SystemExit):
        report_taxon_counts(gtdb_parquet, "NonExistent")


def test_report_taxon_counts_reports_the_dereplicated_size(gtdb_parquet, capsys):
    """--derep-rank used to be ignored here, so the number quoted was not the number
    a pull would return."""
    report_taxon_counts(gtdb_parquet, "Enterobacteriaceae",
                        args=make_args(derep_rank="genus"))
    out = capsys.readouterr().out
    assert "The rank 'family' has 3 Enterobacteriaceae entries." in out
    assert "Dereplicated at 'genus', that would be 2 genome(s)." in out


def test_report_taxon_counts_without_derep_says_nothing_about_derep(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "Escherichia", args=make_args())
    assert "Dereplicated at" not in capsys.readouterr().out


def test_report_taxon_counts_coarser_derep_rank_is_explained(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "Escherichia",
                        args=make_args(derep_rank="phylum"))
    out = capsys.readouterr().out
    assert "broader rank than" in out
    assert "--derep-rank rank must be the same or finer" in out


# ─── report_rank_counts_for_taxon ─────────────────────────────────────────────

def test_rank_counts_for_taxon_starts_at_its_own_rank(gtdb_parquet, capsys):
    report_rank_counts_for_taxon(gtdb_parquet, "Enterobacteriaceae")
    out = capsys.readouterr().out
    assert "under 'Enterobacteriaceae'" in out
    assert "domain" not in out          # coarser than the taxon's own rank
    assert "genus      2" in out        # Escherichia + Salmonella


def test_rank_counts_for_taxon_not_found_exits(gtdb_parquet):
    with pytest.raises(SystemExit):
        report_rank_counts_for_taxon(gtdb_parquet, "NonExistent")


# ─── report_unique_taxa_counts_of_all_ranks ───────────────────────────────────

def test_unique_taxa_counts_of_all_ranks_prints_all_ranks(gtdb_parquet, capsys):
    report_unique_taxa_counts_of_all_ranks(gtdb_parquet)
    out = capsys.readouterr().out
    for rank in ("domain", "phylum", "class", "order", "family", "genus", "species"):
        assert rank in out


def test_unique_taxa_counts_of_all_ranks_with_refseq_rep(gtdb_parquet, capsys):
    """RefSeq reference genomes are a sparse subset, so they get their own (genuinely
    different) counts table."""
    report_unique_taxa_counts_of_all_ranks(gtdb_parquet, representatives_source="refseq")
    out = capsys.readouterr().out
    assert "RefSeq reference genomes" in out
    assert "Num. Unique Ref. Taxa" in out


def test_unique_taxa_counts_of_all_ranks_with_gtdb_rep_shows_note_not_table(gtdb_parquet, capsys):
    """
    GTDB representatives can't change a *unique-taxa* count at any rank (every taxon
    has a representative), so the command prints an explanatory note instead of a
    duplicate second table.
    """
    report_unique_taxa_counts_of_all_ranks(gtdb_parquet, representatives_source="gtdb")
    out = capsys.readouterr().out
    assert "doesn't change these counts" in out
    assert "Num. Unique Ref. Taxa" not in out
    assert out.count("Num. Unique Taxa") == 1


# ─── get_accessions_from_gtdb (orchestrator) ──────────────────────────────────

def test_orchestrator_get_table_copies_file_and_exits(gtdb_parquet, tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.chdir(out_dir)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_table=True))
    assert (out_dir / "gtdb-arc-and-bac-metadata.tsv").exists()


def test_orchestrator_get_taxon_counts_without_taxon_exits(tmp_path):
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=str(tmp_path)):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_taxon_counts=True, target_taxon=None))


def test_orchestrator_conflicting_rep_flags_exits(tmp_path):
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=str(tmp_path)):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(
                gtdb_representatives_only=True,
                refseq_reference_genomes_only=True,
            ))


def test_orchestrator_get_rank_counts_exits(gtdb_parquet, capsys):
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_rank_counts=True))
    assert "domain" in capsys.readouterr().out


def test_orchestrator_rank_counts_gtdb_reps_shows_note(gtdb_parquet, capsys):
    """Full command path: --get-rank-counts --gtdb-representatives-only prints one
    table plus a note (GTDB reps can't change unique-taxa counts), not a duplicate."""
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_rank_counts=True,
                                               gtdb_representatives_only=True))
    out = capsys.readouterr().out
    assert "doesn't change these counts" in out
    assert "Num. Unique Ref. Taxa" not in out


def test_orchestrator_rank_counts_refseq_reps_shows_subtable(gtdb_parquet, capsys):
    """RefSeq reference genomes DO differ, so --get-rank-counts -R keeps the second
    table (contrast with -G, which shows a note)."""
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_rank_counts=True,
                                               refseq_reference_genomes_only=True))
    out = capsys.readouterr().out
    assert "RefSeq reference genomes" in out
    assert "Num. Unique Ref. Taxa" in out


def test_orchestrator_get_taxon_counts_exits(gtdb_parquet, capsys):
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(get_taxon_counts=True, target_taxon="Escherichia"))
    assert "2" in capsys.readouterr().out


def test_orchestrator_get_accessions_exits(gtdb_parquet, tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.chdir(out_dir)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(target_taxon="Escherichia"))
    assert (out_dir / "gtdb-escherichia-genus-accs.txt").exists()


def test_orchestrator_gtdb_rep_only_filters_to_rep_genomes(gtdb_parquet, tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.chdir(out_dir)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(
                target_taxon="Escherichia",
                gtdb_representatives_only=True,
            ))
    accs = (out_dir / "gtdb-escherichia-genus-gtdb-rep-accs.txt").read_text().splitlines()
    assert len(accs) == 1


def test_orchestrator_refseq_ref_only_filters_to_ref_genomes(gtdb_parquet, tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    monkeypatch.chdir(out_dir)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(
                target_taxon="Escherichia",
                refseq_reference_genomes_only=True,
            ))
    accs = (out_dir / "gtdb-escherichia-genus-refseq-rep-accs.txt").read_text().splitlines()
    assert len(accs) == 1


# ─── deliberate-read orchestrator path (new: reads only what each path needs) ──

def test_orchestrator_taxon_slice_matches_full_filter(gtdb_parquet, tmp_path, monkeypatch):
    """
    The deliberate taxon-slice read must yield the SAME accessions as loading the
    whole table and filtering in pandas (guards the resolve->select seam that
    silently mis-sliced when resolve_taxon's rank was treated as a list).
    """
    monkeypatch.chdir(tmp_path)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(target_taxon="Escherichia"))

    new_accs = sorted((tmp_path / "gtdb-escherichia-genus-accs.txt").read_text().split())
    full = pd.DataFrame(ROWS)
    old_accs = sorted(full[full["genus"] == "Escherichia"]["ncbi_genbank_assembly_accession"].tolist())
    assert new_accs == old_accs


def test_orchestrator_metadata_tsv_carries_all_asset_columns(gtdb_parquet, tmp_path, monkeypatch):
    """The metadata TSV should carry every column in the parquet asset."""
    monkeypatch.chdir(tmp_path)
    with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
               return_value=gtdb_parquet):
        with pytest.raises(SystemExit):
            get_accessions_from_gtdb(make_args(target_taxon="Escherichia"))

    meta = pd.read_csv(tmp_path / "gtdb-escherichia-genus-metadata.tsv", sep="\t")
    asset_cols = set(pq.ParquetFile(gtdb_parquet).schema_arrow.names)
    assert set(meta.columns) == asset_cols


def test_resolve_or_exit_returns_canonical_and_rank(gtdb_parquet):
    """
    _resolve_or_exit returns (canonical, rank, domain) -- rank is a single string,
    not a list, and domain is the taxon's sole assigned domain.
    """
    canonical, rank, domain = _resolve_or_exit(gtdb_parquet, "escherichia")
    assert canonical == "Escherichia"
    assert rank == "genus"
    assert domain == "Bacteria"


def test_resolve_or_exit_all_returns_none_rank(gtdb_parquet):
    assert _resolve_or_exit(gtdb_parquet, "all") == ("all", None, None)


def test_resolve_or_exit_ambiguous_exits(gtdb_parquet, tmp_path):
    """A name at multiple ranks with no -r must print-and-exit, not raise."""
    rows = [{
        "accession": "x", "domain": "Bacteria", "phylum": "P", "class": "Dup",
        "order": "O", "family": "Dup", "genus": "G", "species": "G s",
        "gtdb_representative": "t", "ncbi_refseq_category": "na",
        "ncbi_genbank_assembly_accession": "GCA_1.1",
    }]
    p = tmp_path / "amb.parquet"
    pq.write_table(pa.Table.from_pandas(pd.DataFrame(rows), preserve_index=False), str(p))
    with pytest.raises(SystemExit):
        _resolve_or_exit(str(p), "Dup")


# ─── `-t all` with --derep-rank ───────────────────────────────────────────────

def test_all_with_derep_rank_dereplicates_within_each_domain(gtdb_parquet, tmp_path,
                                                             monkeypatch, capsys):
    """
    'all' has no rank of its own to group beneath, so it used to reject --derep-rank
    (and before that, silently ignore it). It now runs the per-domain selection and
    merges -- the fixture holds Archaea and Bacteria.
    """
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all", derep_rank="domain")
    accs = (tmp_path / "gtdb-arc-and-bac-accessions.txt").read_text().splitlines()
    assert len(accs) == 2      # one genome per domain
    out = capsys.readouterr().out
    assert "Dereplicating within each domain (Archaea, Bacteria)" in out
    assert (tmp_path / "gtdb-arc-and-bac-metadata.tsv").exists()


def test_all_with_derep_rank_at_a_finer_rank(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all", derep_rank="genus")
    accs = (tmp_path / "gtdb-arc-and-bac-accessions.txt").read_text().splitlines()
    assert len(accs) == 3      # Haloarcula + Escherichia + Salmonella


def test_all_without_derep_rank_is_still_a_bulk_dump(gtdb_parquet, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(gtdb_parquet, target_taxon="all")
    accs = (tmp_path / "gtdb-arc-and-bac-accessions.txt").read_text().splitlines()
    assert len(accs) == 5


def test_all_taxon_counts_report_the_dereplicated_size(gtdb_parquet, capsys):
    report_taxon_counts(gtdb_parquet, "all", args=make_args(derep_rank="genus"))
    assert "Dereplicated within each domain, that would be 3 genome(s)." in capsys.readouterr().out


# ─── NCBI liveness screening of the GTDB snapshot ─────────────────────────────
#
# GTDB is a fixed snapshot and still lists assemblies NCBI has since suppressed or
# removed. NCBI drops those from its assembly summary, so absence there == suppressed.
# Screening here means a user doesn't get handed accessions that can't be downloaded.

class TestSuppressedGtdbEntriesAreScreenedOut:

    @pytest.fixture
    def screen_missing_one(self, tmp_path):
        """An NCBI asset where GCA_000005845.2 (a non-rep E. coli) has been suppressed."""
        path = tmp_path / "ncbi-screen-partial.parquet"
        live = [a for a in ALL_GENBANK_ACCS if a != "GCA_000005845.2"]
        _write_ncbi_screen_asset(path, live)
        return str(path)

    def _run(self, gtdb_parquet, screen_path, tmp_path, monkeypatch, **args):
        monkeypatch.chdir(tmp_path)
        with patch("bit.modules.gtdb.get_accessions_from_gtdb.get_gtdb_data",
                   return_value=gtdb_parquet), \
             patch("bit.modules.gtdb.get_accessions_from_gtdb.ncbi_screen_path",
                   return_value=screen_path):
            with pytest.raises(SystemExit):
                get_accessions_from_gtdb(make_args(**args))

    def _accessions_written(self, tmp_path):
        # taxon pulls write "*-accs.txt"; the bulk 'all' dump writes "*-accessions.txt"
        hits = sorted(tmp_path.glob("*accs.txt")) + sorted(tmp_path.glob("*accessions.txt"))
        assert len(hits) == 1, f"expected one accessions file, got {hits}"
        return hits[0].read_text().split()

    def test_a_suppressed_genome_is_dropped_from_a_taxon_pull(
            self, gtdb_parquet, screen_missing_one, tmp_path, monkeypatch):
        self._run(gtdb_parquet, screen_missing_one, tmp_path, monkeypatch,
                  target_taxon="Escherichia")
        accs = self._accessions_written(tmp_path)
        assert "GCA_000005845.2" not in accs
        assert "GCA_000001405.1" in accs

    def test_nothing_is_dropped_when_everything_is_live(
            self, gtdb_parquet, live_ncbi_screen, tmp_path, monkeypatch):
        self._run(gtdb_parquet, live_ncbi_screen, tmp_path, monkeypatch,
                  target_taxon="Escherichia")
        accs = self._accessions_written(tmp_path)
        assert set(accs) == {"GCA_000001405.1", "GCA_000005845.2"}

    def test_the_bulk_all_dump_is_screened_too(
            self, gtdb_parquet, screen_missing_one, tmp_path, monkeypatch):
        self._run(gtdb_parquet, screen_missing_one, tmp_path, monkeypatch,
                  target_taxon="all")
        accs = self._accessions_written(tmp_path)
        assert "GCA_000005845.2" not in accs
        assert len(accs) == len(ALL_GENBANK_ACCS) - 1

    def test_the_derep_path_is_screened_too(
            self, gtdb_parquet, tmp_path, monkeypatch):
        """
        With the better E. coli suppressed, derep must fall back to the surviving one
        rather than picking the dead genome and then dropping it (which would leave the
        species contributing nothing).
        """
        path = tmp_path / "ncbi-screen-no-ref.parquet"
        live = [a for a in ALL_GENBANK_ACCS if a != "GCA_000001405.1"]
        _write_ncbi_screen_asset(path, live)
        self._run(gtdb_parquet, str(path), tmp_path, monkeypatch,
                  target_taxon="Escherichia", derep_rank="species")
        accs = self._accessions_written(tmp_path)
        assert accs == ["GCA_000005845.2"]

    def test_a_suppressed_genome_is_reported_to_the_user(
            self, gtdb_parquet, screen_missing_one, tmp_path, monkeypatch, capsys):
        self._run(gtdb_parquet, screen_missing_one, tmp_path, monkeypatch,
                  target_taxon="Escherichia")
        out = capsys.readouterr().out
        assert "no longer available at NCBI" in out
