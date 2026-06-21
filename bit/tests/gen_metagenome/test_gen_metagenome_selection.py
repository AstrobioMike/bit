import pandas as pd
import pytest

from bit.modules.gen_mg.selection import (
    _select_gtdb_one_per_rank,
    _normalize_gtdb_rows,
    _normalize_ncbi_rows,
    _split_ncbi_lineage,
    merge_sources,
    user_filled_groups,
    SCHEMA,
    RANKS,
)


# ─── fixtures ──────────────────────────────────────────────────────────────

def gtdb_row(acc, dom, gen, sp, rep, refseq="", size=4_000_000):
    return {
        "ncbi_genbank_assembly_accession": acc,
        "domain": dom, "phylum": "P", "class": "C", "order": "O",
        "family": "F", "genus": gen, "species": sp,
        "gtdb_representative": rep, "ncbi_refseq_category": refseq,
        "genome_size": size,
    }


@pytest.fixture
def gtdb_tab():
    return pd.DataFrame([
        gtdb_row("GCF_001.1", "Bacteria", "Genus_A", "sp A1", "t", "reference genome", 5_000_000),
        gtdb_row("GCF_002.1", "Bacteria", "Genus_A", "sp A2", "t", "", 4_500_000),
        gtdb_row("GCF_003.1", "Bacteria", "Genus_B", "sp B1", "t", "", 3_000_000),
        gtdb_row("GCF_004.1", "Bacteria", "Pseudomonas", "sp P1", "t", "", 6_300_000),
        gtdb_row("GCF_005.1", "Archaea", "Genus_C", "sp C1", "t", "", 1_800_000),
        gtdb_row("GCF_006.1", "Bacteria", "Genus_B", "sp B2", "f", "", 3_100_000),  # non-rep
    ])


def ncbi_row(acc, name, taxid, size, src="genbank"):
    return {
        "accession": acc, "organism_name": name, "tax_id": taxid,
        "assembly_name": "x", "assembly_level": "Complete",
        "refseq_category": "NA", "total_sequence_length": str(size),
        "checkm_completeness": "NA", "checkm_contamination": "NA",
        "source_database": src,
    }


# ─── GTDB selection ────────────────────────────────────────────────────────

def test_one_per_genus_picks_one_per_group(gtdb_tab):
    sel, warnings = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus", seed=1)
    assert set(sel["genus"]) == {"Genus_A", "Genus_B", "Pseudomonas", "Genus_C"}
    assert len(sel) == 4
    assert warnings == []


def test_within_group_prefers_refseq_reference(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus", seed=1)
    a_pick = sel[sel["genus"] == "Genus_A"]["ncbi_genbank_assembly_accession"].iloc[0]
    assert a_pick == "GCF_001.1"


def test_species_derep_excludes_non_representatives(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="species", seed=1)
    accs = set(sel["ncbi_genbank_assembly_accession"])
    assert "GCF_006.1" not in accs       # the gtdb_representative == "f" row
    assert len(sel) == 5


def test_domain_slice(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus",
                                       domains=["Bacteria"], seed=1)
    assert set(sel["genus"]) == {"Genus_A", "Genus_B", "Pseudomonas"}


def test_num_exceeds_groups_warns_and_returns_fewer(gtdb_tab):
    sel, warnings = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus",
                                              num_genomes=10, seed=1)
    assert len(sel) == 4
    assert any("only 4 unique" in w for w in warnings)


def test_num_caps_below_groups(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus",
                                       num_genomes=2, seed=1)
    assert len(sel) == 2
    assert sel["genus"].nunique() == 2


def test_absent_domain_warns_and_empty(gtdb_tab):
    sel, warnings = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus",
                                              domains=["Eukaryota"], seed=1)
    assert len(sel) == 0
    assert any("Eukaryota" in w for w in warnings)


def test_derep_off_returns_all_representatives(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="off", seed=1)
    assert len(sel) == 5      # all reps, no derep
    assert "GCF_006.1" not in set(sel["ncbi_genbank_assembly_accession"])


def test_exclude_groups_skips_user_filled(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus",
                                       exclude_groups={"Pseudomonas"}, seed=1)
    assert "Pseudomonas" not in set(sel["genus"])


def test_invalid_derep_rank_raises(gtdb_tab):
    with pytest.raises(ValueError):
        _select_gtdb_one_per_rank(gtdb_tab, derep_rank="kingdom")


# ─── normalization ─────────────────────────────────────────────────────────

def test_normalize_gtdb_schema_and_sizes(gtdb_tab):
    sel, _ = _select_gtdb_one_per_rank(gtdb_tab, derep_rank="genus", seed=1)
    norm = _normalize_gtdb_rows(sel)
    assert list(norm.columns) == SCHEMA
    assert (~norm["user_supplied"]).all()
    # gtdb_* ranks populated, ncbi_* left NA (resolved later by taxonomy layer)
    assert (norm["gtdb_domain"] != "NA").all()
    assert set(norm["gtdb_domain"]) <= {"Bacteria", "Archaea"}
    assert (norm["ncbi_domain"] == "NA").all()
    # dropped columns are absent
    assert "source_db" not in norm.columns
    assert "taxonomy_source" not in norm.columns
    assert "metadata_genome_size" not in norm.columns


def test_normalize_ncbi_without_lineage_gives_na_ranks():
    tab = pd.DataFrame([
        ncbi_row("GCA_900.1", "Saccharomyces cerevisiae", "4932", 12_100_000),
        ncbi_row("GCA_901.1", "Aspergillus niger", "5061", 34_000_000),
    ])
    norm = _normalize_ncbi_rows(tab, lineage_lookup=None, user_supplied=True)
    # both rank sets NA pre-resolution; gtdb_* always NA for NCBI-normalized rows
    assert (norm["ncbi_domain"] == "NA").all()
    assert (norm["gtdb_domain"] == "NA").all()
    assert norm["user_supplied"].all()
    assert "source_db" not in norm.columns
    assert "metadata_genome_size" not in norm.columns


def test_normalize_ncbi_with_lineage_populates_ranks():
    tab = pd.DataFrame([ncbi_row("GCA_900.1", "Saccharomyces cerevisiae", "4932", 12_100_000)])
    lineage = {"GCA_900.1": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
                            "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"}
    norm = _normalize_ncbi_rows(tab, lineage_lookup=lineage, user_supplied=True)
    row = norm.iloc[0]
    assert row["ncbi_domain"] == "Eukaryota"
    assert row["ncbi_genus"] == "Saccharomyces"
    assert row["ncbi_species"] == "Saccharomyces cerevisiae"
    assert row["gtdb_domain"] == "NA"      # GTDB filled later if in GTDB


def test_split_lineage_handles_missing_ranks():
    vals = _split_ncbi_lineage("Eukaryota;Ascomycota")
    assert vals["domain"] == "Eukaryota"
    assert vals["phylum"] == "Ascomycota"
    assert vals["species"] == "NA"
    # empty / None -> all NA
    assert all(v == "NA" for v in _split_ncbi_lineage(None).values())


# ─── merge ─────────────────────────────────────────────────────────────────

def test_merge_pools_sources():
    g = _normalize_gtdb_rows(pd.DataFrame([gtdb_row("GCF_001.1", "Bacteria", "Genus_A", "sp", "t")]))
    u = _normalize_ncbi_rows(pd.DataFrame([ncbi_row("GCA_900.1", "Fungus", "1", 12_000_000)]),
                             user_supplied=True)
    merged, warnings = merge_sources(g, u)
    assert set(merged["accession"]) == {"GCF_001.1", "GCA_900.1"}


def test_merge_dedup_prefers_user_row():
    g = _normalize_gtdb_rows(pd.DataFrame([gtdb_row("GCF_001.1", "Bacteria", "Genus_A", "sp", "t")]))
    dup = g.copy()
    dup["user_supplied"] = True
    merged, warnings = merge_sources(g, dup)
    assert len(merged) == 1
    assert merged["user_supplied"].iloc[0]
    assert any("duplicate" in w for w in warnings)


def test_merge_empty_returns_warning():
    merged, warnings = merge_sources(None, None)
    assert len(merged) == 0
    assert warnings


def test_user_filled_groups():
    # user_filled_groups reads gtdb_<rank> (generative derep is GTDB-based); in
    # the real flow fill_gtdb_taxonomy populates this before exclusion runs.
    u = _normalize_ncbi_rows(
        pd.DataFrame([ncbi_row("GCA_999.1", "Pseudomonas aeruginosa", "287", 6_300_000)]),
        lineage_lookup=None, user_supplied=True)
    # before GTDB resolution, gtdb ranks are NA -> no groups filled
    assert user_filled_groups(u, "genus") == set()
    # simulate fill_gtdb_taxonomy having resolved the GTDB genus
    u = u.copy()
    u["gtdb_genus"] = "Pseudomonas"
    assert user_filled_groups(u, "genus") == {"Pseudomonas"}
    assert user_filled_groups(u, "off") == set()
    assert user_filled_groups(None, "genus") == set()
