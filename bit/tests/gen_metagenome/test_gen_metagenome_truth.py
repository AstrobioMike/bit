import os
import pandas as pd
import pytest

from bit.modules.gen_metagenome.truth import (
    build_gen_reads_args,
    build_per_genome_table,
    build_per_rank_tables,
    build_read_truth,
    RANKS,
)


@pytest.fixture
def merged():
    """ a normalized + abundance-assigned 3-genome table (incl. an unresolved euk). """
    return pd.DataFrame([
        dict(accession="GCF_A", source_db="GTDB", domain="Bacteria", phylum="P",
             **{"class": "C"}, order="O", family="F", genus="Genus_A", species="sp A",
             genome_size=1_000_000, assigned_rel_abundance=0.5, assigned_coverage=10.0,
             assigned_reads=5000, taxonomy_source="GTDB", user_supplied=False),
        dict(accession="GCF_B", source_db="GTDB", domain="Bacteria", phylum="P",
             **{"class": "C"}, order="O", family="F", genus="Genus_B", species="sp B",
             genome_size=500_000, assigned_rel_abundance=0.3, assigned_coverage=12.0,
             assigned_reads=3000, taxonomy_source="GTDB", user_supplied=False),
        dict(accession="GCA_E", source_db="genbank", domain="NA", phylum="NA",
             **{"class": "NA"}, order="NA", family="NA", genus="NA", species="NA",
             genome_size=300_000, assigned_rel_abundance=0.2, assigned_coverage=8.0,
             assigned_reads=2000, taxonomy_source="NCBI", user_supplied=True),
    ])


# ─── gen-reads args ────────────────────────────────────────────────────────

def test_build_gen_reads_args_defaults():
    args = build_gen_reads_args(["a.fasta", "b.fasta"], "cov.tsv", "pref")
    assert args.coverage == "cov.tsv"          # coverage-specified mode
    assert args.source_tsv is True             # provenance always on
    assert args.read_length == 150
    assert args.type == "paired-end"


def test_build_gen_reads_args_long_default_length():
    args = build_gen_reads_args(["a.fasta"], "cov.tsv", "pref", read_type="long")
    assert args.read_length == 5000


# ─── per-genome truth ──────────────────────────────────────────────────────

def test_per_genome_table_merges_mutation(merged):
    mut = {
        "GCF_A": {"rate": 0.02, "num_substitutions": 100, "num_indels": 5, "num_total_changes": 105},
        "GCF_B": {"rate": 0.0},
        "GCA_E": {"rate": 0.0},
    }
    pg = build_per_genome_table(merged, mut)
    a = pg[pg["accession"] == "GCF_A"].iloc[0]
    assert a["mutation_rate"] == 0.02
    assert a["num_substitutions"] == 100
    assert a["num_total_changes"] == 105
    b = pg[pg["accession"] == "GCF_B"].iloc[0]
    assert b["mutation_rate"] == 0.0
    assert b["num_substitutions"] == 0


def test_per_genome_table_column_order(merged):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]})
    assert pg.columns[0] == "accession"
    for r in RANKS:
        assert r in pg.columns
    assert "assigned_rel_abundance" in pg.columns
    assert "taxonomy_source" in pg.columns


# ─── per-rank collapse ─────────────────────────────────────────────────────

def test_per_rank_sums_to_one(merged):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]})
    tables = build_per_rank_tables(pg)
    for rank, tab in tables.items():
        assert abs(tab["rel_abundance"].sum() - 1.0) < 1e-9


def test_per_rank_groups_shared_taxa(merged):
    # add a 4th genome sharing Genus_A so genus-collapse merges them
    extra = merged.iloc[[0]].copy()
    extra["accession"] = "GCF_A2"
    extra["species"] = "sp A2"
    extra["assigned_rel_abundance"] = 0.1
    merged.loc[0, "assigned_rel_abundance"] = 0.4
    pg = build_per_genome_table(pd.concat([merged, extra], ignore_index=True),
                                {a: {"rate": 0.0} for a in list(merged["accession"]) + ["GCF_A2"]})
    genus = build_per_rank_tables(pg)["genus"]
    a_row = genus[genus["genus"] == "Genus_A"]
    assert len(a_row) == 1
    assert abs(a_row["rel_abundance"].iloc[0] - 0.5) < 1e-9   # 0.4 + 0.1


def test_per_rank_keeps_na_bucket(merged):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]})
    genus = build_per_rank_tables(pg)["genus"]
    assert "NA" in set(genus["genus"])        # unresolved euk preserved


# ─── read-level truth ──────────────────────────────────────────────────────

def test_read_truth_joins_taxonomy(merged, tmp_path):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]})

    # a minimal gen-reads source TSV referencing two genomes
    src = tmp_path / "reads-read-sources.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\t/x/GCF_A.fasta\tGCF_A_ctg0\t10\t160\t+\tfalse\n"
        "r1/2\t/x/GCF_A.fasta\tGCF_A_ctg0\t300\t450\t-\tfalse\n"
        "r2/1\t/x/GCA_E.fasta\tGCA_E_ctg0\t5\t155\t+\tfalse\n"
    )
    fasta_to_acc = {
        "/x/GCF_A.fasta": "GCF_A", "GCF_A.fasta": "GCF_A",
        "/x/GCA_E.fasta": "GCA_E", "GCA_E.fasta": "GCA_E",
    }
    out = tmp_path / "read-truth.tsv"
    rt = build_read_truth(str(src), fasta_to_acc, pg, str(out))

    assert (rt["accession"] != "NA").all()
    # GCF_A reads carry Genus_A; euk read carries NA
    assert (rt[rt["accession"] == "GCF_A"]["genus"] == "Genus_A").all()
    assert (rt[rt["accession"] == "GCA_E"]["genus"] == "NA").all()
    # coordinates carried through
    assert set(["contig", "start", "end", "strand", "wrapped"]).issubset(rt.columns)
    assert out.exists()


def test_read_truth_basename_fallback(merged, tmp_path):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]})
    src = tmp_path / "src.tsv"
    # source_fasta given as bare basename; lookup has only basename key
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\tGCF_B.fasta\tGCF_B_ctg0\t1\t151\t+\tfalse\n"
    )
    rt = build_read_truth(str(src), {"GCF_B.fasta": "GCF_B"}, pg, str(tmp_path / "o.tsv"))
    assert rt["accession"].iloc[0] == "GCF_B"
    assert rt["genus"].iloc[0] == "Genus_B"
