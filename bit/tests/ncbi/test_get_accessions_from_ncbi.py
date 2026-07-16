import os

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore
from types import SimpleNamespace

import bit.modules.ncbi.get_accessions_from_ncbi as M
from bit.modules.ncbi.get_accessions_from_ncbi import (
    parse_assembly_levels,
    get_accessions_from_ncbi,
    _apply_filters,
    _select_by_taxid,
)


# The module reads the hosted NCBI Parquet; fixtures are Parquet built on the same
# columns the reader requests (_COLUMNS) plus the per-rank *_taxid columns the taxid
# path needs.

_COLS = [
    "assembly_accession", "organism_name", "taxid", "asm_name", "assembly_level",
    "refseq_category", "checkm_completeness", "checkm_contamination", "genome_size",
    "domain", "phylum", "class", "order", "family", "genus", "species",
    "domain_taxid", "phylum_taxid", "class_taxid", "order_taxid", "family_taxid",
    "genus_taxid", "species_taxid",
]


def _row(acc, genus="Alteromonas", level="Complete Genome", refseq="",
         genus_taxid="28108", **over):
    r = {c: "NA" for c in _COLS}
    r.update({
        "assembly_accession": acc, "organism_name": f"{genus} sp.", "taxid": "28108",
        "asm_name": "ASM1", "assembly_level": level, "refseq_category": refseq,
        "checkm_completeness": "99.0", "checkm_contamination": "0.5",
        "genome_size": "4000000",
        "domain": "Bacteria", "phylum": "Pseudomonadota", "class": "Gammaproteobacteria",
        "order": "Alteromonadales", "family": "Alteromonadaceae", "genus": genus,
        "species": f"{genus} macleodii",
        "domain_taxid": "2", "phylum_taxid": "1224", "class_taxid": "1236",
        "order_taxid": "135622", "family_taxid": "72275", "genus_taxid": genus_taxid,
        "species_taxid": "28108",
    })
    r.update(over)
    return r


def _make_table(tmp_path, rows):
    data = {c: pa.array([str(r[c]) for r in rows]) for c in _COLS}
    p = tmp_path / "ncbi-data.parquet"
    pq.write_table(pa.table(data), str(p))
    # the date stamp the module now reports at startup sits beside the table
    (tmp_path / "date-retrieved.txt").write_text("2026,01,05\n")
    return str(p)


def _args(**kw):
    base = dict(target_taxon=None, target_rank=None, source="refseq",
                refseq_reference_genomes_only=False, assembly_level=None,
                get_taxon_counts=False, get_rank_counts=False)
    base.update(kw)
    return SimpleNamespace(**base)


@pytest.fixture
def table(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1", refseq="reference genome"),
        _row("GCF_000000002.1", refseq=""),
        _row("GCA_000000003.1", refseq="", level="Scaffold"),
    ])
    # skip the download; point the module at our fixture
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    return p


# --- parse_assembly_levels ------------------------------------------------

def test_parse_assembly_levels_maps_friendly_names():
    assert parse_assembly_levels(["complete", "contig"]) == ["Complete Genome", "Contig"]


def test_parse_assembly_levels_comma_string():
    assert parse_assembly_levels("complete,scaffold") == ["Complete Genome", "Scaffold"]


def test_parse_assembly_levels_empty():
    assert parse_assembly_levels(None) == []
    assert parse_assembly_levels([]) == []


def test_parse_assembly_levels_unknown_raises():
    with pytest.raises(ValueError):
        parse_assembly_levels(["banana"])


# --- taxon path -----------------------------------------------------------

def test_taxon_writes_both_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas"))
    assert (tmp_path / "NCBI-Alteromonas-refseq-accessions.txt").exists()
    assert (tmp_path / "NCBI-Alteromonas-refseq-metadata.tsv").exists()


def test_taxon_accs_file_contents(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # source=both so GCA rows aren't prefix-filtered out
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", source="both"))
    accs = (tmp_path / "NCBI-Alteromonas-accessions.txt").read_text().split()
    assert "GCF_000000001.1" in accs
    assert "GCA_000000003.1" in accs


def test_reference_only_filters_and_names_file(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas",
                                   refseq_reference_genomes_only=True))
    accs = (tmp_path / "NCBI-Alteromonas-refseq-ref-accessions.txt").read_text().split()
    assert accs == ["GCF_000000001.1"]     # only the reference genome


def test_case_insensitive_match(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="alteromonas"))
    assert (tmp_path / "NCBI-alteromonas-refseq-accessions.txt").exists()


def test_taxon_not_found_exits_without_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(target_taxon="Nonexistentia"))
    assert not list(tmp_path.glob("NCBI-*"))


def test_counts_mode_reports_and_writes_nothing(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(target_taxon="Alteromonas", get_taxon_counts=True))
    out = capsys.readouterr().out
    assert "genome(s) under" in out
    assert not list(tmp_path.glob("NCBI-*"))


# --- taxid path -----------------------------------------------------------

def test_taxid_path_matches_by_rank_taxid(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # 28108 is the genus_taxid in the fixture
    get_accessions_from_ncbi(_args(target_taxon="28108", source="both"))
    assert (tmp_path / "NCBI-28108-accessions.txt").exists()


def test_select_by_taxid_matches_specific_rank(tmp_path):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", genus_taxid="561"),
        _row("GCF_2.1", genus_taxid="590"),
    ])
    tab = _select_by_taxid(p, "561")
    accs = tab.column("assembly_accession").to_pylist()
    assert accs == ["GCF_1.1"]


# --- filters --------------------------------------------------------------

def test_source_refseq_keeps_only_gcf(table, tmp_path):
    tab = pq.read_table(table, columns=_COLS)
    out = _apply_filters(tab, _args(source="refseq"))
    accs = out.column("assembly_accession").to_pylist()
    assert all(a.startswith("GCF_") for a in accs)


def test_source_genbank_keeps_only_gca(table, tmp_path):
    tab = pq.read_table(table, columns=_COLS)
    out = _apply_filters(tab, _args(source="genbank"))
    accs = out.column("assembly_accession").to_pylist()
    assert all(a.startswith("GCA_") for a in accs)


def test_assembly_level_filter(table, tmp_path):
    tab = pq.read_table(table, columns=_COLS)
    out = _apply_filters(tab, _args(source="both", assembly_level=["scaffold"]))
    levels = set(out.column("assembly_level").to_pylist())
    assert levels == {"Scaffold"}


# --- --get-rank-counts (whole-database, no taxon) -------------------------

def test_get_rank_counts_prints_all_ranks(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    # a line per rank, with the header
    assert "Num. Unique Taxa" in out
    for rank in ["domain", "phylum", "genus", "species"]:
        assert rank in out
    # 3 fixture rows: 1 genus (Alteromonas) -> genus count 1
    assert "genus      1" in out


def test_get_rank_counts_writes_no_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True))
    assert not list(tmp_path.glob("NCBI-*"))


def test_get_rank_counts_reps_only_adds_second_table(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True,
                                       refseq_reference_genomes_only=True))
    out = capsys.readouterr().out
    assert "RefSeq reference genomes" in out
    assert "Num. Unique Ref. Taxa" in out


def test_get_rank_counts_needs_no_taxon(table, tmp_path, monkeypatch):
    """The rank-counts path must run with target_taxon=None (whole-DB)."""
    monkeypatch.chdir(tmp_path)
    args = _args(get_rank_counts=True)
    assert args.target_taxon is None
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(args)


def test_get_rank_counts_honors_source_refseq(table, tmp_path, monkeypatch, capsys):
    """--source refseq counts only GCF_ rows. Fixture: 2 GCF + 1 GCA, all Alteromonas
    -> genus count is 1 either way, but the header labels the source."""
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True, source="refseq"))
    out = capsys.readouterr().out
    assert "(RefSeq)" in out


def test_get_rank_counts_source_genbank_excludes_gcf(tmp_path, monkeypatch, capsys):
    """A genus present ONLY in a GCF row must not be counted under --source genbank."""
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1", genus="RefseqOnly"),   # refseq
        _row("GCA_000000002.1", genus="GenbankOne"),   # genbank
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True, source="genbank"))
    out = capsys.readouterr().out
    assert "(GenBank)" in out
    # genbank scope: exactly 1 genus (GenbankOne); the GCF-only genus is excluded
    assert "genus      1" in out


def test_get_rank_counts_source_both_counts_everything(tmp_path, monkeypatch, capsys):
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1", genus="RefseqOnly"),
        _row("GCA_000000002.1", genus="GenbankOne"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True, source="both"))
    out = capsys.readouterr().out
    assert "(all)" in out
    assert "genus      2" in out          # both genera counted


# --- -t all (every genome, skips resolution) ------------------------------

def test_all_writes_every_row(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", genus="Aaa"),
        _row("GCF_2.1", genus="Bbb"),
        _row("GCA_3.1", genus="Ccc"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    # source=both so all three rows survive
    get_accessions_from_ncbi(_args(target_taxon="all", source="both"))
    accs = (tmp_path / "NCBI-all-accessions.txt").read_text().split()
    assert sorted(accs) == ["GCA_3.1", "GCF_1.1", "GCF_2.1"]


def test_all_default_refseq_scopes_to_gcf(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", genus="Aaa"),
        _row("GCA_2.1", genus="Bbb"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all"))    # default source refseq
    accs = (tmp_path / "NCBI-all-refseq-accessions.txt").read_text().split()
    assert accs == ["GCF_1.1"]                             # GCA excluded


def test_all_is_case_insensitive(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [_row("GCF_1.1")])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="ALL"))
    # filename normalizes to lowercase 'all' regardless of input casing
    assert (tmp_path / "NCBI-all-refseq-accessions.txt").exists()


def test_all_reference_only(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", refseq="reference genome"),
        _row("GCF_2.1", refseq=""),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all", refseq_reference_genomes_only=True))
    accs = (tmp_path / "NCBI-all-refseq-ref-accessions.txt").read_text().split()
    assert accs == ["GCF_1.1"]


def test_all_counts_mode(tmp_path, monkeypatch, capsys):
    p = _make_table(tmp_path, [_row("GCF_1.1"), _row("GCF_2.1")])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(target_taxon="all", get_taxon_counts=True))
    assert "2 genome(s) under all genomes" in capsys.readouterr().out
    assert not list(tmp_path.glob("NCBI-*"))


# --- date reporting -------------------------------------------------------

def test_read_date_retrieved_formats_nicely(tmp_path):
    from bit.modules.ncbi.get_ncbi_assembly_data import read_date_retrieved
    (tmp_path / "date-retrieved.txt").write_text("2026,01,05\n")
    assert read_date_retrieved(str(tmp_path)) == "Jan 05, 2026"


def test_read_date_retrieved_falls_back_to_raw(tmp_path):
    from bit.modules.ncbi.get_ncbi_assembly_data import read_date_retrieved
    (tmp_path / "date-retrieved.txt").write_text("not-a-date\n")
    assert read_date_retrieved(str(tmp_path)) == "not-a-date"


def test_startup_prints_date_and_update_hint(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        get_accessions_from_ncbi(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    assert "Date NCBI data retrieved: Jan 05, 2026" in out
