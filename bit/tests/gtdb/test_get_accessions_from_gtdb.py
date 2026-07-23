import pytest # type: ignore
import pandas as pd # type: ignore
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
from argparse import Namespace
from unittest.mock import patch
from bit.modules.gtdb.get_accessions_from_gtdb import (
    report_gtdb_version_info,
    copy_gtdb_table,
    find_ranks_for_taxon,
    get_accessions,
    get_unique_taxon_counts,
    get_unique_taxa_counts_of_all_ranks,
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
    },
    {
        "accession": "RS_GCF_000009905.1",
        "domain": "Archaea", "phylum": "Halobacteriota", "class": "Halobacteria",
        "order": "Halobacteriales", "family": "Haloarculaceae", "genus": "Haloarcula",
        "species": "Haloarcula marismortui",
        "gtdb_representative": "f", "ncbi_refseq_category": "na",
        "ncbi_genbank_assembly_accession": "GCA_000009905.1",
    },
    {
        "accession": "GB_GCA_000001405.1",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Escherichia",
        "species": "Escherichia coli",
        "gtdb_representative": "t", "ncbi_refseq_category": "reference genome",
        "ncbi_genbank_assembly_accession": "GCA_000001405.1",
    },
    {
        "accession": "GB_GCA_000005845.2",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Escherichia",
        "species": "Escherichia coli",
        "gtdb_representative": "f", "ncbi_refseq_category": "na",
        "ncbi_genbank_assembly_accession": "GCA_000005845.2",
    },
    {
        "accession": "GB_GCA_000006925.2",
        "domain": "Bacteria", "phylum": "Proteobacteria", "class": "Gammaproteobacteria",
        "order": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Salmonella",
        "species": "Salmonella enterica",
        "gtdb_representative": "t", "ncbi_refseq_category": "reference genome",
        "ncbi_genbank_assembly_accession": "GCA_000006925.2",
    },
]


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


# ─── find_ranks_for_taxon ─────────────────────────────────────────────────────

def test_find_ranks_for_taxon_single_rank(gtdb_tab):
    assert find_ranks_for_taxon("Escherichia", gtdb_tab) == ["genus"]


def test_find_ranks_for_taxon_domain(gtdb_tab):
    assert find_ranks_for_taxon("Archaea", gtdb_tab) == ["domain"]


def test_find_ranks_for_taxon_not_found(gtdb_tab):
    assert find_ranks_for_taxon("NonExistent", gtdb_tab) == []


def test_find_ranks_for_taxon_multiple_ranks():
    multi_tab = pd.DataFrame([{
        "domain": "Bacteria", "phylum": "A", "class": "Duplicate",
        "order": "C", "family": "Duplicate", "genus": "G", "species": "G sp",
        "accession": "acc1", "gtdb_representative": "t",
        "ncbi_refseq_category": "na", "ncbi_genbank_assembly_accession": "GCA_001.1",
    }])
    result = find_ranks_for_taxon("Duplicate", multi_tab)
    assert "class" in result
    assert "family" in result


# ─── get_accessions ───────────────────────────────────────────────────────────

def test_get_accessions_all_writes_accs_file(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("all", gtdb_tab)
    accs = (tmp_path / "gtdb-arc-and-bac-accessions.txt").read_text().splitlines()
    assert len(accs) == 5


def test_get_accessions_all_no_rep_source_no_tab_file(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("all", gtdb_tab)
    assert not (tmp_path / "gtdb-arc-and-bac-refseq-rep-metadata.tsv").exists()


def test_get_accessions_all_with_rep_source_writes_both_files(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rep_tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]
    get_accessions("all", gtdb_tab, gtdb_rep_tab=rep_tab, representatives_source="gtdb")
    assert (tmp_path / "gtdb-arc-and-bac-refseq-rep-metadata.tsv").exists()
    accs = (tmp_path / "gtdb-arc-and-bac-refseq-rep-accessions.txt").read_text().splitlines()
    assert len(accs) == 3


def test_get_accessions_specific_taxon_writes_accs_file(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("Escherichia", gtdb_tab)
    accs = (tmp_path / "gtdb-escherichia-genus-accs.txt").read_text().splitlines()
    assert len(accs) == 2


def test_get_accessions_specific_taxon_writes_metadata_file(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("Escherichia", gtdb_tab)
    assert (tmp_path / "gtdb-escherichia-genus-metadata.tsv").exists()


def test_get_accessions_rank_specified_explicitly(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("Escherichia", gtdb_tab, rank="genus")
    accs = (tmp_path / "gtdb-escherichia-genus-accs.txt").read_text().splitlines()
    assert len(accs) == 2


def test_get_accessions_species_with_space_replaced_by_dash(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions("Escherichia coli", gtdb_tab, rank="species")
    assert (tmp_path / "gtdb-escherichia-coli-species-accs.txt").exists()


def test_get_accessions_with_rep_source_in_filename(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rep_tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]
    get_accessions("Escherichia", gtdb_tab, gtdb_rep_tab=rep_tab, representatives_source="gtdb")
    assert (tmp_path / "gtdb-escherichia-genus-gtdb-rep-accs.txt").exists()
    accs = (tmp_path / "gtdb-escherichia-genus-gtdb-rep-accs.txt").read_text().splitlines()
    assert len(accs) == 1


def test_get_accessions_multi_rank_no_rank_exits():
    multi_tab = pd.DataFrame([{
        "domain": "Bacteria", "phylum": "A", "class": "Duplicate",
        "order": "C", "family": "Duplicate", "genus": "G", "species": "G sp",
        "accession": "acc1", "gtdb_representative": "t",
        "ncbi_refseq_category": "na", "ncbi_genbank_assembly_accession": "GCA_001.1",
    }])
    with pytest.raises(SystemExit):
        get_accessions("Duplicate", multi_tab)


def test_get_accessions_taxon_not_found_exits(gtdb_tab, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions("NonExistent", gtdb_tab)


# ─── get_unique_taxon_counts ──────────────────────────────────────────────────

def test_get_unique_taxon_counts_all_total(gtdb_tab, capsys):
    get_unique_taxon_counts("all", gtdb_tab)
    assert "5" in capsys.readouterr().out


def test_get_unique_taxon_counts_all_with_rep_tab(gtdb_tab, capsys):
    rep_tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]
    get_unique_taxon_counts("all", gtdb_tab, gtdb_rep_tab=rep_tab, representatives_source="RefSeq")
    out = capsys.readouterr().out
    assert "5" in out   # total
    assert "3" in out   # rep count


def test_get_unique_taxon_counts_specific_taxon(gtdb_tab, capsys):
    get_unique_taxon_counts("Escherichia", gtdb_tab)
    assert "2" in capsys.readouterr().out


def test_get_unique_taxon_counts_not_found_exits(gtdb_tab):
    with pytest.raises(SystemExit):
        get_unique_taxon_counts("NonExistent", gtdb_tab)


# ─── get_unique_taxa_counts_of_all_ranks ──────────────────────────────────────

def test_get_unique_taxa_counts_of_all_ranks_prints_all_ranks(gtdb_tab, capsys):
    get_unique_taxa_counts_of_all_ranks(gtdb_tab)
    out = capsys.readouterr().out
    for rank in ("domain", "phylum", "class", "order", "family", "genus", "species"):
        assert rank in out


def test_get_unique_taxa_counts_of_all_ranks_with_refseq_rep(gtdb_tab, capsys):
    """RefSeq reference genomes are a sparse subset, so they get their own (genuinely
    different) counts table."""
    rep_tab = gtdb_tab[gtdb_tab["ncbi_refseq_category"] == "reference genome"]
    get_unique_taxa_counts_of_all_ranks(gtdb_tab, gtdb_rep_tab=rep_tab, representatives_source="RefSeq")
    out = capsys.readouterr().out
    assert "RefSeq reference genomes" in out
    assert "Num. Unique Ref. Taxa" in out


def test_get_unique_taxa_counts_of_all_ranks_with_gtdb_rep_shows_note_not_table(gtdb_tab, capsys):
    """
    GTDB representatives can't change a *unique-taxa* count at any rank (every taxon
    has a representative), so the command prints an explanatory note instead of a
    duplicate second table.
    """
    rep_tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]
    get_unique_taxa_counts_of_all_ranks(gtdb_tab, gtdb_rep_tab=rep_tab, representatives_source="gtdb")
    out = capsys.readouterr().out
    assert "doesn't change these counts" in out
    # and NOT a second table
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
    accs = (out_dir / "gtdb-escherichia-genus-RefSeq-rep-accs.txt").read_text().splitlines()
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
    """_resolve_or_exit returns (canonical, rank) -- rank is a single string, not a list."""
    canonical, rank = _resolve_or_exit(gtdb_parquet, "escherichia")
    assert canonical == "Escherichia"
    assert rank == "genus"


def test_resolve_or_exit_all_returns_none_rank(gtdb_parquet):
    assert _resolve_or_exit(gtdb_parquet, "all") == ("all", None)


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
