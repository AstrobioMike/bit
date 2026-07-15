import pytest # type: ignore
from bit.modules.taxonomy.tax_ranks import NA
from bit.modules.taxonomy.tax_select import (
    SOURCES,
    AmbiguousTaxon,
    TaxonNotFound,
    resolve_taxon,
    select_accessions,
    select_by_taxid,
    live_accession_cores,
)
from bit.modules.gtdb.build_gtdb_data_parquet import (
    build_table as gtdb_build,
    write_parquet as gtdb_write,
)


GHDR = ("accession\tgtdb_taxonomy\tncbi_genbank_assembly_accession\tncbi_taxid\t"
        "gtdb_representative\tncbi_refseq_category\tcheckm2_completeness\t"
        "checkm2_contamination\tcheckm_completeness\tcheckm_contamination\t"
        "genome_size\tcontig_count\tgc_count\tgc_percentage\tambiguous_bases\t"
        "coding_bases\tcoding_density\n")


def _grow(acc, tax, gca, rep="t", comp="99.0", cont="0.5"):
    return (f"{acc}\t{tax}\t{gca}\t1\t{rep}\tna\t{comp}\t{cont}\t90\t1\t"
            f"1000\t1\t500\t50.0\t0\t900\t0.9\n")


def _tax(d="Bacteria", p="P1", c="C1", o="O1", f="F1", g="G1", s="G1 s1"):
    return f"d__{d};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"


@pytest.fixture
def gtdb(tmp_path):
    """
    Three classes under Bacteria:
      C1 -- two representatives of differing quality (tests the quality pick)
      C2 -- one representative
      C3 -- present only as a NON-representative (tests the reps pool)
    Plus one genome with an UNNAMED class (tests the NA-group guard).
    """
    src = tmp_path / "bac120.tsv"
    src.write_text(
        GHDR
        + _grow("RS_GCF_1.1", _tax(c="C1", g="G1", s="G1 a"), "GCA_1.1", comp="80.0", cont="4.0")
        + _grow("RS_GCF_2.1", _tax(c="C1", g="G2", s="G2 b"), "GCA_2.1", comp="99.5", cont="0.1")
        + _grow("RS_GCF_3.1", _tax(c="C2", g="G3", s="G3 c"), "GCA_3.1", comp="95.0", cont="1.0")
        + _grow("RS_GCF_4.1", _tax(c="C3", g="G4", s="G4 d"), "GCA_4.1", rep="f")
        + _grow("RS_GCF_5.1", _tax(c="", g="G5", s="G5 e"), "GCA_5.1")
    )
    arc = tmp_path / "ar53.tsv"
    arc.write_text(GHDR)
    out = tmp_path / "gtdb.parquet"
    gtdb_write(gtdb_build(str(arc), str(src)), str(out), row_group_size=2)
    return str(out)


def test_resolve_taxon_is_case_insensitive(gtdb):
    canonical, rank = resolve_taxon(gtdb, "bacteria")
    assert canonical == "Bacteria"
    assert rank == "domain"


def test_unknown_taxon_raises(gtdb):
    with pytest.raises(TaxonNotFound):
        resolve_taxon(gtdb, "Nonexistus")


def test_explicit_rank_that_doesnt_match_raises(gtdb):
    with pytest.raises(TaxonNotFound):
        resolve_taxon(gtdb, "Bacteria", rank="phylum")


def test_homonym_raises_ambiguous(tmp_path):
    """A name living at two ranks must force the caller to disambiguate."""
    src = tmp_path / "b.tsv"
    # 'X' is used as BOTH an order and a family
    src.write_text(GHDR + _grow("RS_GCF_9.1", _tax(o="X", f="X"), "GCA_9.1"))
    arc = tmp_path / "a.tsv"
    arc.write_text(GHDR)
    out = tmp_path / "g.parquet"
    gtdb_write(gtdb_build(str(arc), str(src)), str(out))

    with pytest.raises(AmbiguousTaxon) as e:
        resolve_taxon(str(out), "X")
    assert set(e.value.ranks_found) == {"order", "family"}

    # ...but an explicit rank resolves it
    assert resolve_taxon(str(out), "X", rank="family") == ("X", "family")


def test_select_accessions_returns_ncbi_accs(gtdb):
    accs = select_accessions(gtdb, "gtdb", "domain", "Bacteria")
    assert sorted(accs) == ["GCA_1.1", "GCA_2.1", "GCA_3.1", "GCA_4.1", "GCA_5.1"]


def test_reps_only_excludes_non_representatives(gtdb):
    accs = select_accessions(gtdb, "gtdb", "domain", "Bacteria", reps_only=True)
    assert "GCA_4.1" not in accs      # the only rep='f' genome


def test_source_specs_cover_both_and_differ_only_in_columns():
    assert set(SOURCES) == {"gtdb", "ncbi"}
    assert SOURCES["gtdb"].rep_filter == ("gtdb_representative", "t")
    assert SOURCES["ncbi"].rep_filter == ("refseq_category", "reference genome")
    # BOTH sources now have checkm and a refseq_category, so the pick key is shared
    assert SOURCES["gtdb"].quality_cols == ("checkm2_completeness", "checkm2_contamination")
    assert SOURCES["ncbi"].quality_cols == ("checkm_completeness", "checkm_contamination")
    assert SOURCES["gtdb"].ref_col == "ncbi_refseq_category"
    assert SOURCES["ncbi"].ref_col == "refseq_category"
    # only NCBI needs the euk fallback (GTDB is prokaryotes only, always has checkm2)
    assert SOURCES["ncbi"].level_col == "assembly_level"
    assert SOURCES["ncbi"].contig_col == "contig_count"
    assert SOURCES["ncbi"].size_col == "genome_size_ungapped"
    assert SOURCES["ncbi"].size_fallback_col == "genome_size"
    assert SOURCES["gtdb"].level_col is None


def test_live_accession_cores_is_prefix_and_version_agnostic(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq_
    t = pa.table({"assembly_accession": pa.array(["GCF_000005845.2", "GCA_000012285.1"])})
    p = tmp_path / "n.parquet"
    pq_.write_table(t, str(p))
    cores = live_accession_cores(str(p))
    # GTDB-style prefixed/versioned accessions must resolve to the same cores
    from bit.modules.taxonomy.tax_ranks import accession_core
    assert accession_core("RS_GCF_000005845.2") in cores
    assert accession_core("GB_GCA_000005845.3") in cores   # version drift still matches
    assert accession_core("GCA_999999999.1") not in cores


def test_gtdb_spec_now_carries_size_and_contig_columns():
    assert SOURCES["gtdb"].size_col == "genome_size"
    assert SOURCES["gtdb"].contig_col == "contig_count"


# --- select_by_taxid (NCBI lineage-taxid selection) -----------------------

import pyarrow as _pa # type: ignore
import pyarrow.parquet as _pq # type: ignore
from bit.modules.taxonomy.tax_select import select_by_taxid

# minimal NCBI-shaped columns select_by_taxid needs: an accession + the per-rank
# {rank}_taxid columns it filters on. (The real table has ~28 columns; the reader
# only ever touches these, so the fixture stays small on purpose.)
_NCBI_TAXID_COLS = ["assembly_accession", "genus_taxid", "family_taxid"]


def _ncbi_taxid_table(tmp_path, rows):
    """rows: list of (accession, genus_taxid, family_taxid) as strings."""
    data = {c: [] for c in _NCBI_TAXID_COLS}
    for acc, g, f in rows:
        data["assembly_accession"].append(acc)
        data["genus_taxid"].append(g)
        data["family_taxid"].append(f)
    p = tmp_path / "ncbi.parquet"
    _pq.write_table(_pa.table({c: _pa.array(data[c]) for c in _NCBI_TAXID_COLS}), str(p))
    return str(p)


def test_select_by_taxid_returns_matching_accessions(tmp_path):
    p = _ncbi_taxid_table(tmp_path, [
        ("GCF_1.1", "561", "543"),      # genus Escherichia (561)
        ("GCF_2.1", "561", "543"),      # same genus
        ("GCF_3.1", "590", "543"),      # different genus (Salmonella 590)
    ])
    accs = select_by_taxid(p, "genus", 561)
    assert sorted(accs) == ["GCF_1.1", "GCF_2.1"]


def test_select_by_taxid_no_match_returns_empty(tmp_path):
    p = _ncbi_taxid_table(tmp_path, [("GCF_1.1", "561", "543")])
    assert select_by_taxid(p, "genus", 99999) == []


def test_select_by_taxid_accepts_int_or_str(tmp_path):
    """the filter stringifies the taxid, so an int and its str form match the same rows."""
    p = _ncbi_taxid_table(tmp_path, [("GCF_1.1", "561", "543")])
    assert select_by_taxid(p, "genus", 561) == select_by_taxid(p, "genus", "561")


def test_select_by_taxid_is_rank_specific(tmp_path):
    """
    A taxid value must only match within its own rank column. Here 543 is a FAMILY
    taxid; querying it at genus rank (genus_taxid) must NOT return the row, even
    though 543 appears in the table (under family_taxid).
    """
    p = _ncbi_taxid_table(tmp_path, [("GCF_1.1", "561", "543")])
    assert select_by_taxid(p, "family", 543) == ["GCF_1.1"]   # right column
    assert select_by_taxid(p, "genus", 543) == []             # wrong column, no cross-match
