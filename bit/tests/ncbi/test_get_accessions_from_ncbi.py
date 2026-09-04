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
    "genome_size_ungapped", "contig_count",
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
        "genome_size": "4000000", "genome_size_ungapped": "4000000",
        "contig_count": "1",
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
    base = dict(target_taxon=None, target_rank=None, ncbi_section="refseq",
                refseq_reference_genomes_only=False, assembly_level=None,
                get_taxon_counts=False, get_rank_counts=False, get_table=False,
                derep_rank="off")
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


def test_cli_assembly_level_is_repeatable():
    # repeated --assembly-level must accumulate (argparse append), not overwrite down
    # to the last value; parse then turns the tokens into the NCBI display strings
    from bit.cli.get_accessions_from_ncbi import build_parser
    args = build_parser().parse_args(
        ["-t", "Alteromonas", "--assembly-level", "complete",
         "--assembly-level", "contig"])
    assert args.assembly_level == ["complete", "contig"]
    assert parse_assembly_levels(args.assembly_level) == ["Complete Genome", "Contig"]


# --- taxon path -----------------------------------------------------------

def test_taxon_writes_both_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas"))
    assert (tmp_path / "ncbi-alteromonas-genus-refseq-accs.txt").exists()
    assert (tmp_path / "ncbi-alteromonas-genus-refseq-metadata.tsv").exists()


def test_taxon_accs_file_contents(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # source=both so GCA rows aren't prefix-filtered out
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", ncbi_section="both"))
    accs = (tmp_path / "ncbi-alteromonas-genus-accs.txt").read_text().split()
    assert "GCF_000000001.1" in accs
    assert "GCA_000000003.1" in accs


def test_reference_only_filters_and_names_file(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas",
                                   refseq_reference_genomes_only=True))
    accs = (tmp_path / "ncbi-alteromonas-genus-refseq-ref-accs.txt").read_text().split()
    assert accs == ["GCF_000000001.1"]     # only the reference genome


def test_case_insensitive_match(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="alteromonas"))
    assert (tmp_path / "ncbi-alteromonas-genus-refseq-accs.txt").exists()


def test_taxon_not_found_exits_without_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Nonexistentia"))
    assert not list(tmp_path.glob("ncbi-*-accs.txt"))
    assert not list(tmp_path.glob("ncbi-*-metadata.tsv"))

def test_counts_mode_reports_and_writes_nothing(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Alteromonas", get_taxon_counts=True))
    out = capsys.readouterr().out
    # per-rank format, matching get-accs-from-gtdb
    assert "The rank 'genus' has" in out
    assert "Alteromonas entries" in out
    assert not list(tmp_path.glob("ncbi-*-accs.txt"))
    assert not list(tmp_path.glob("ncbi-*-metadata.tsv"))


def _match_lines(out):
    return [l.strip() for l in out.splitlines() if "entries" in l]


def test_counts_mode_match_counts_are_not_collapsed_by_derep(table, tmp_path,
                                                             monkeypatch, capsys):
    """
    --get-taxon-counts reports how many genomes MATCH the filters, so --derep-rank
    must not change that number -- but the dereplicated size is now reported
    alongside it (it used to be omitted entirely).
    """
    monkeypatch.chdir(tmp_path)
    seen = []
    for derep in ("off", "species", "auto"):
        with pytest.raises(SystemExit):
            get_accessions_from_ncbi(_args(target_taxon="Alteromonas",
                                           get_taxon_counts=True, derep_rank=derep))
        seen.append(capsys.readouterr().out)

    assert _match_lines(seen[0]) == _match_lines(seen[1]) == _match_lines(seen[2])
    assert "Dereplicated at" not in seen[0]          # derep off
    assert "Dereplicated at 'species'" in seen[1]


def test_counts_mode_reports_the_dereplicated_size(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Alteromonas", ncbi_section="both",
                                       get_taxon_counts=True, derep_rank="species"))
    out = capsys.readouterr().out
    # all 3 fixture rows share one species -> one genome survives
    assert "Dereplicated at 'species', that would be 1 genome(s)." in out


def test_counts_mode_derep_size_matches_what_a_pull_returns(table, tmp_path,
                                                            monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Alteromonas", ncbi_section="both",
                                       get_taxon_counts=True, derep_rank="species"))
    reported = capsys.readouterr().out

    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", ncbi_section="both",
                                   derep_rank="species"))
    accs = (tmp_path / "ncbi-alteromonas-genus-accs.txt").read_text().split()
    assert f"that would be {len(accs)} genome(s)." in reported


def test_counts_mode_reflects_source_filter(table, tmp_path, monkeypatch, capsys):
    # refseq scoping drops the one GCA row, so the genus count is 2 (not 3), and every
    # count line carries the generic "(after any specified filters)" tag rather than an
    # enumerated scope note
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Alteromonas",
                                       get_taxon_counts=True, ncbi_section="refseq"))
    out = capsys.readouterr().out
    assert "(after any specified filters)" in out
    assert "The rank 'genus' has 2 Alteromonas entries" in out


def test_counts_mode_reports_every_rank_a_name_occurs_at(tmp_path, monkeypatch, capsys):
    """An ambiguous name is informative here rather than an error: one line per rank."""
    p = _make_table(tmp_path, [
        _row("GCF_1.1", genus="Dup", family="Dup"),
        _row("GCF_2.1", genus="Other", family="Dup"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="Dup", get_taxon_counts=True,
                                       ncbi_section="both"))
    out = capsys.readouterr().out
    assert "The rank 'family' has 2 Dup entries" in out
    assert "The rank 'genus' has 1 Dup entries" in out


# --- taxid path -----------------------------------------------------------

def test_taxid_path_matches_by_rank_taxid(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # 28108 is the genus_taxid in the fixture
    get_accessions_from_ncbi(_args(target_taxon="28108", ncbi_section="both"))
    assert (tmp_path / "ncbi-taxid-28108-accs.txt").exists()


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
    out = _apply_filters(tab, _args(ncbi_section="refseq"))
    accs = out.column("assembly_accession").to_pylist()
    assert all(a.startswith("GCF_") for a in accs)


def test_source_genbank_keeps_only_gca(table, tmp_path):
    tab = pq.read_table(table, columns=_COLS)
    out = _apply_filters(tab, _args(ncbi_section="genbank"))
    accs = out.column("assembly_accession").to_pylist()
    assert all(a.startswith("GCA_") for a in accs)


def test_apply_filters_no_longer_handles_assembly_level(table, tmp_path):
    """
    --assembly-level is a PRE-filter now (see the derep-ordering tests below), so
    _apply_filters is source scoping only and must not silently drop rows by level.
    """
    tab = pq.read_table(table, columns=_COLS)
    out = _apply_filters(tab, _args(ncbi_section="both", assembly_level=["scaffold"]))
    assert out.num_rows == tab.num_rows


def test_assembly_level_still_filters_end_to_end(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", ncbi_section="both",
                                   assembly_level=["scaffold"]))
    accs = (tmp_path / "ncbi-alteromonas-genus-accs.txt").read_text().split()
    assert accs == ["GCA_000000003.1"]


def test_assembly_level_still_filters_for_all(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all", ncbi_section="both",
                                   assembly_level=["scaffold"]))
    accs = (tmp_path / "ncbi-all-accs.txt").read_text().split()
    assert accs == ["GCA_000000003.1"]


def test_assembly_level_still_filters_for_a_taxid(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="28108", ncbi_section="both",
                                   assembly_level=["scaffold"]))
    accs = (tmp_path / "ncbi-taxid-28108-accs.txt").read_text().split()
    assert accs == ["GCA_000000003.1"]


def test_assembly_level_is_applied_before_dereplication(tmp_path, monkeypatch):
    """
    FamB holds a Complete Genome and a higher-quality Scaffold. Filtering on level
    AFTER dereplication picked the Scaffold as FamB's best and then deleted it, so
    FamB contributed nothing -- instead of contributing the best genome that actually
    meets the requested level.
    """
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1", genus="GenA", family="FamA",
             checkm_completeness="90.0"),
        _row("GCF_000000003.1", genus="GenB", family="FamB",
             checkm_completeness="90.0"),
        _row("GCA_000000004.1", genus="GenB", family="FamB", level="Scaffold",
             checkm_completeness="99.9"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonadales", ncbi_section="both",
                                   derep_rank="family", assembly_level=["complete"]))
    accs = sorted((tmp_path / "ncbi-alteromonadales-order-accs.txt").read_text().split())
    assert accs == ["GCF_000000001.1", "GCF_000000003.1"]


# --- --get-rank-counts (whole-database, no taxon) -------------------------

def test_get_rank_counts_prints_all_ranks(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _run(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    # a line per rank, with the header
    assert "Num. Unique Taxa" in out
    for rank in ["domain", "phylum", "genus", "species"]:
        assert rank in out
    # 3 fixture rows: 1 genus (Alteromonas) -> genus count 1
    assert "genus      1" in out


def test_get_rank_counts_writes_no_files(table, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _run(_args(get_rank_counts=True))
    assert not list(tmp_path.glob("ncbi-*-accs.txt"))
    assert not list(tmp_path.glob("ncbi-*-metadata.tsv"))


def test_get_rank_counts_reps_only_adds_second_table(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _run(_args(get_rank_counts=True,
                                       refseq_reference_genomes_only=True))
    out = capsys.readouterr().out
    assert "RefSeq reference genomes" in out
    assert "Num. Unique Ref. Taxa" in out


def test_get_rank_counts_needs_no_taxon(table, tmp_path, monkeypatch):
    """The rank-counts path must run with target_taxon=None (whole-DB)."""
    monkeypatch.chdir(tmp_path)
    args = _args(get_rank_counts=True)
    assert args.target_taxon is None
    _run(args)


def test_get_rank_counts_honors_source_refseq(table, tmp_path, monkeypatch, capsys):
    """--source refseq counts only GCF_ rows. Fixture: 2 GCF + 1 GCA, all Alteromonas
    -> genus count is 1 either way, but the header labels the source."""
    monkeypatch.chdir(tmp_path)
    _run(_args(get_rank_counts=True, ncbi_section="refseq"))
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
    _run(_args(get_rank_counts=True, ncbi_section="genbank"))
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
    _run(_args(get_rank_counts=True, ncbi_section="both"))
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
    get_accessions_from_ncbi(_args(target_taxon="all", ncbi_section="both"))
    accs = (tmp_path / "ncbi-all-accs.txt").read_text().split()
    assert sorted(accs) == ["GCA_3.1", "GCF_1.1", "GCF_2.1"]


def test_all_default_refseq_scopes_to_gcf(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", genus="Aaa"),
        _row("GCA_2.1", genus="Bbb"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all"))    # default source refseq
    accs = (tmp_path / "ncbi-all-refseq-accs.txt").read_text().split()
    assert accs == ["GCF_1.1"]                             # GCA excluded


def test_all_is_case_insensitive(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [_row("GCF_1.1")])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="ALL"))
    # filename normalizes to lowercase 'all' regardless of input casing
    assert (tmp_path / "ncbi-all-refseq-accs.txt").exists()


def test_all_reference_only(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_1.1", refseq="reference genome"),
        _row("GCF_2.1", refseq=""),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all", refseq_reference_genomes_only=True))
    accs = (tmp_path / "ncbi-all-refseq-ref-accs.txt").read_text().split()
    assert accs == ["GCF_1.1"]


def test_all_counts_mode(tmp_path, monkeypatch, capsys):
    p = _make_table(tmp_path, [_row("GCF_1.1"), _row("GCF_2.1")])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="all", get_taxon_counts=True))
    assert "2 genome(s) under all genomes" in capsys.readouterr().out
    assert not list(tmp_path.glob("ncbi-*-accs.txt"))
    assert not list(tmp_path.glob("ncbi-*-metadata.tsv"))


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
    _run(_args(get_rank_counts=True))
    out = capsys.readouterr().out
    assert "Using NCBI assembly-data retrieved: Jan 05, 2026" in out


# --- preflight: --derep-rank applicability --------------------------------

def test_all_with_derep_rank_is_allowed_now(capsys):
    from bit.cli.get_accessions_from_ncbi import check_derep_rank_is_applicable
    check_derep_rank_is_applicable(_args(target_taxon="all", derep_rank="family"))


def test_all_with_derep_rank_includes_eukaryotes(tmp_path, monkeypatch, capsys):
    """
    The domain list is read from the table, not hardcoded to the prokaryotic pair --
    NCBI carries eukaryotes and they must not be silently dropped.
    """
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1"),
        _row("GCF_000000002.1"),
        _row("GCF_000000020.1", genus="Saccharomyces", domain="Eukaryota",
             phylum="Ascomycota"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="all", ncbi_section="both",
                                   derep_rank="domain"))
    accs = (tmp_path / "ncbi-all-accs.txt").read_text().split()
    assert len(accs) == 2                      # one bacterium, one eukaryote
    assert "GCF_000000020.1" in accs
    assert "Bacteria, Eukaryota" in capsys.readouterr().out


def test_derep_rank_with_a_taxid_is_refused(capsys):
    from bit.cli.get_accessions_from_ncbi import check_derep_rank_is_applicable
    with pytest.raises(SystemExit):
        check_derep_rank_is_applicable(_args(target_taxon="28108", derep_rank="family"))
    assert "`--derep-rank` can't be applied with a taxid" in capsys.readouterr().out


def test_derep_rank_off_with_a_taxid_is_allowed():
    from bit.cli.get_accessions_from_ncbi import check_derep_rank_is_applicable
    check_derep_rank_is_applicable(_args(target_taxon="28108", derep_rank="off"))


def test_derep_rank_with_a_named_taxon_is_allowed():
    from bit.cli.get_accessions_from_ncbi import check_derep_rank_is_applicable
    check_derep_rank_is_applicable(_args(target_taxon="Alteromonas", derep_rank="genus"))


################################################################################
# 'all' and the genomes it can't reach
#
# The NCBI table carries rows with NO assigned domain (viral and
# metagenome/uncultured entries). 'all' is expanded to the table's DOMAINS, so those
# rows are unreachable by it -- which has to be both CONSISTENT (the same pool with
# and without --derep-rank) and REPORTED (an unscoped --get-rank-counts otherwise
# quotes a number no pull will produce).
################################################################################

def _viral_row(acc):
    """A domain-less row whose class exists nowhere else in the fixture."""
    return _row(acc, domain="NA", phylum="Uroviricota", **{"class": "Caudoviricetes"},
                order="NA", family="NA", genus="NA", species="Phage sp.")


@pytest.fixture
def table_with_viruses(tmp_path, monkeypatch):
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1", refseq="reference genome"),
        _row("GCF_000000002.1", genus="Pseudoalteromonas", genus_taxid="53246"),
        _viral_row("GCF_000000030.1"),
        _viral_row("GCF_000000031.1"),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    return p


def _accs(tmp_path, name):
    return [l.strip() for l in open(tmp_path / name) if l.strip()]


def _run(args):
    """The counts/table paths exit; the pull paths return. Tolerate either."""
    try:
        get_accessions_from_ncbi(args)
    except SystemExit:
        pass


def test_all_without_derep_excludes_domainless_genomes(table_with_viruses, tmp_path):
    """
    'all' has to mean the same pool whether or not --derep-rank is set. It used to be
    a whole-table dump without derep (viruses in) and a domain walk with it (viruses
    out), so an unrelated flag decided whether viral genomes appeared.
    """
    _run(_args(target_taxon="all"))
    accs = _accs(tmp_path, "ncbi-all-refseq-accs.txt")
    assert sorted(accs) == ["GCF_000000001.1", "GCF_000000002.1"]


def test_all_reports_the_genomes_it_left_behind(table_with_viruses, capsys):
    _run(_args(target_taxon="all", derep_rank="class"))
    out = capsys.readouterr().out
    assert "2 genome(s)" in out
    assert "no assigned domain" in out


def test_rank_counts_with_all_are_scoped_to_what_all_pulls(table_with_viruses,
                                                           tmp_path, capsys):
    """
    The reconciliation: the genus count reported for `-t all` must equal the number of
    accessions `-t all --derep-rank genus` writes. Unscoped it counted the viral rows
    too, which is the shape of "--get-rank-counts says 326 but we get 280".
    """
    _run(_args(target_taxon="all", get_rank_counts=True))
    assert "genus      2" in capsys.readouterr().out

    _run(_args(target_taxon="all", derep_rank="genus"))
    assert len(_accs(tmp_path, "ncbi-all-refseq-accs.txt")) == 2


def test_rank_counts_without_a_taxon_still_counts_the_whole_table(table_with_viruses,
                                                                  capsys):
    """Unscoped `--get-rank-counts` is documented as the whole database -- unchanged."""
    _run(_args(get_rank_counts=True))
    assert "class      2" in capsys.readouterr().out   # includes Caudoviricetes


def test_domainless_genomes_are_still_reachable_by_name(table_with_viruses, tmp_path):
    """They're excluded from 'all', not from the tool."""
    _run(_args(target_taxon="Uroviricota"))
    accs = _accs(tmp_path, "ncbi-uroviricota-phylum-refseq-accs.txt")
    assert sorted(accs) == ["GCF_000000030.1", "GCF_000000031.1"]


def test_source_both_no_longer_warns_about_the_overlap(table, tmp_path, monkeypatch, capsys):
    # the `--ncbi-section both` overlap note was dropped now that --derep-rank
    # defaults to auto, so `both` no longer returns near-duplicate assemblies
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", ncbi_section="both"))
    out = capsys.readouterr().out
    assert "overlap between RefSeq and GenBank" not in out


def test_single_source_says_nothing_about_overlap(table, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    get_accessions_from_ncbi(_args(target_taxon="Alteromonas", ncbi_section="refseq"))
    assert "overlap between RefSeq and GenBank" not in capsys.readouterr().out


def test_eukaryote_alias_resolves(tmp_path, monkeypatch):
    """`-t eukarya` (and friends) resolve to Eukaryota wherever a taxon is taken."""
    p = _make_table(tmp_path, [
        _row("GCF_000000001.1"),
        _row("GCF_000000020.1", domain="Eukaryota", phylum="Ascomycota",
             genus="Saccharomyces", **{"class": "Saccharomycetes"}),
    ])
    monkeypatch.setattr(M, "ncbi_table_path", lambda **k: p)
    monkeypatch.chdir(tmp_path)
    _run(_args(target_taxon="eukarya"))
    assert _accs(tmp_path, "ncbi-eukaryota-domain-refseq-accs.txt") == ["GCF_000000020.1"]
