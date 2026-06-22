import os
import pandas as pd # type: ignore
import pytest # type: ignore

from bit.modules.gen_mg.truth import (
    build_gen_reads_args,
    build_per_genome_table,
    build_per_rank_tables,
    build_read_truth,
    write_truth_outputs,
    RANKS,
    TAXONOMIES,
)


# helpers ────────────────────────────────────────────────────────────────────

def _ranks(prefix, vals):
    """ build {prefix_domain: .., ..} from a 7-value list. """
    return {f"{prefix}_{r}": v for r, v in zip(RANKS, vals)}


@pytest.fixture
def merged():
    """
    a normalized + abundance-assigned 3-genome table carrying BOTH taxonomies.
    GCF_A / GCF_B are GTDB-resolved (gtdb_* and ncbi_* both populated); GCA_E is
    a euk-style user genome with gtdb_* = NA but ncbi_* populated.
    """
    return pd.DataFrame([
        dict(accession="GCF_A", user_supplied=False,
             **_ranks("gtdb", ["Bacteria", "P", "C", "O", "F", "Genus_A", "sp A"]),
             **_ranks("ncbi", ["Bacteria", "Pn", "Cn", "On", "Fn", "Genus_A", "sp A"]),
             downloaded_genome_size=1_000_000, used_genome_size=1_000_000,
             assigned_rel_abundance=0.5, assigned_coverage=10.0, assigned_reads=5000),
        dict(accession="GCF_B", user_supplied=False,
             **_ranks("gtdb", ["Bacteria", "P", "C", "O", "F", "Genus_B", "sp B"]),
             **_ranks("ncbi", ["Bacteria", "Pn", "Cn", "On", "Fn", "Genus_B", "sp B"]),
             downloaded_genome_size=500_000, used_genome_size=500_000,
             assigned_rel_abundance=0.3, assigned_coverage=12.0, assigned_reads=3000),
        dict(accession="GCA_E", user_supplied=True,
             **_ranks("gtdb", ["NA", "NA", "NA", "NA", "NA", "NA", "NA"]),
             **_ranks("ncbi", ["Eukaryota", "Ascomycota", "Saccharomycetes",
                               "Saccharomycetales", "Saccharomycetaceae",
                               "Saccharomyces", "Saccharomyces cerevisiae"]),
             downloaded_genome_size=300_000, used_genome_size=300_000,
             assigned_rel_abundance=0.2, assigned_coverage=8.0, assigned_reads=2000),
    ])


# ─── gen-reads args ────────────────────────────────────────────────────────

def test_build_gen_reads_args_defaults():
    args = build_gen_reads_args(["a.fasta", "b.fasta"], "cov.tsv", "pref")
    assert args.coverage == "cov.tsv"          # coverage-specified mode
    assert args.source_tsv is False            # intermediate; off unless per-read wanted
    assert args.read_length == 150
    assert args.type == "paired-end"


def test_build_gen_reads_args_source_tsv_opt_in():
    args = build_gen_reads_args(["a.fasta"], "cov.tsv", "pref", source_tsv=True)
    assert args.source_tsv is True


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
    pg = build_per_genome_table(merged, mut, "gtdb")
    a = pg[pg["accession"] == "GCF_A"].iloc[0]
    assert a["mutation_rate"] == 0.02
    assert a["num_substitutions"] == 100
    assert a["num_total_changes"] == 105
    b = pg[pg["accession"] == "GCF_B"].iloc[0]
    assert b["mutation_rate"] == 0.0
    assert b["num_substitutions"] == 0


def test_per_genome_table_emits_plain_rank_cols(merged):
    """ the chosen taxonomy's ranks are emitted under plain names, and the
    OTHER taxonomy's columns are not present in the output. """
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    assert pg.columns[0] == "accession"
    for r in RANKS:
        assert r in pg.columns               # plain rank names
    assert "gtdb_genus" not in pg.columns    # source prefix stripped
    assert "ncbi_genus" not in pg.columns    # other taxonomy absent
    assert "source_db" not in pg.columns     # dropped from schema
    assert "metadata_genome_size" not in pg.columns
    assert "downloaded_genome_size" in pg.columns
    assert "used_genome_size" in pg.columns


def test_per_genome_table_taxonomy_selects_source(merged):
    """ gtdb vs ncbi pick different rank values for the same genome. """
    pg_g = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    pg_n = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "ncbi")
    # GCF_A: gtdb phylum 'P' vs ncbi phylum 'Pn'
    assert pg_g.set_index("accession").loc["GCF_A", "phylum"] == "P"
    assert pg_n.set_index("accession").loc["GCF_A", "phylum"] == "Pn"
    # GCA_E euk: gtdb is NA, ncbi is real
    assert pg_g.set_index("accession").loc["GCA_E", "genus"] == "NA"
    assert pg_n.set_index("accession").loc["GCA_E", "genus"] == "Saccharomyces"


def test_per_genome_table_rejects_bad_taxonomy(merged):
    with pytest.raises(ValueError):
        build_per_genome_table(merged, {}, "silva")


# ─── per-rank collapse ─────────────────────────────────────────────────────

def test_per_rank_sums_to_one(merged):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    tables = build_per_rank_tables(pg)
    for rank, tab in tables.items():
        assert abs(tab["rel_abundance"].sum() - 1.0) < 1e-9


def test_per_rank_groups_shared_taxa(merged):
    # add a 4th genome sharing Genus_A so genus-collapse merges them
    extra = merged.iloc[[0]].copy()
    extra["accession"] = "GCF_A2"
    extra["gtdb_species"] = "sp A2"
    extra["assigned_rel_abundance"] = 0.1
    merged.loc[0, "assigned_rel_abundance"] = 0.4
    pg = build_per_genome_table(pd.concat([merged, extra], ignore_index=True),
                                {a: {"rate": 0.0} for a in list(merged["accession"]) + ["GCF_A2"]},
                                "gtdb")
    genus = build_per_rank_tables(pg)["genus"]
    a_row = genus[genus["genus"] == "Genus_A"]
    assert len(a_row) == 1
    assert abs(a_row["rel_abundance"].iloc[0] - 0.5) < 1e-9   # 0.4 + 0.1


def test_per_rank_keeps_na_bucket(merged):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    genus = build_per_rank_tables(pg)["genus"]
    assert "NA" in set(genus["genus"])        # unresolved euk preserved under gtdb


# ─── read-level truth ──────────────────────────────────────────────────────

def test_read_truth_joins_taxonomy(merged, tmp_path):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")

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
    # GCF_A reads carry Genus_A (gtdb); euk read carries NA under gtdb
    assert (rt[rt["accession"] == "GCF_A"]["genus"] == "Genus_A").all()
    assert (rt[rt["accession"] == "GCA_E"]["genus"] == "NA").all()
    # coordinates carried through, and appear before the rank columns
    assert set(["contig", "start", "end", "strand", "wrapped"]).issubset(rt.columns)
    cols = list(rt.columns)
    assert cols.index("wrapped") < cols.index("genus")
    assert out.exists()


def test_read_truth_ncbi_taxonomy(merged, tmp_path):
    """ the euk read gets real NCBI ranks when the ncbi per-genome table is used. """
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "ncbi")
    src = tmp_path / "src.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\t/x/GCA_E.fasta\tGCA_E_ctg0\t5\t155\t+\tfalse\n"
    )
    rt = build_read_truth(str(src), {"/x/GCA_E.fasta": "GCA_E", "GCA_E.fasta": "GCA_E"},
                          pg, str(tmp_path / "o.tsv"))
    assert rt["genus"].iloc[0] == "Saccharomyces"


def test_read_truth_basename_fallback(merged, tmp_path):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    src = tmp_path / "src.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\tGCF_B.fasta\tGCF_B_ctg0\t1\t151\t+\tfalse\n"
    )
    rt = build_read_truth(str(src), {"GCF_B.fasta": "GCF_B"}, pg, str(tmp_path / "o.tsv"))
    assert rt["accession"].iloc[0] == "GCF_B"
    assert rt["genus"].iloc[0] == "Genus_B"


# ─── full split-tree writer ────────────────────────────────────────────────

def test_write_truth_outputs_split_tree(merged, tmp_path):
    src = tmp_path / "reads-read-sources.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\t/x/GCF_A.fasta\tGCF_A_ctg0\t10\t160\t+\tfalse\n"
        "r2/1\t/x/GCA_E.fasta\tGCA_E_ctg0\t5\t155\t+\tfalse\n"
    )
    fasta_to_acc = {"/x/GCF_A.fasta": "GCF_A", "GCF_A.fasta": "GCF_A",
                    "/x/GCA_E.fasta": "GCA_E", "GCA_E.fasta": "GCA_E"}
    written = write_truth_outputs(
        merged, {a: {"rate": 0.0} for a in merged["accession"]}, str(tmp_path),
        read_sources_tsv=str(src), fasta_to_accession=fasta_to_acc, per_read=True)

    # both taxonomy subdirs produced
    assert set(written.keys()) == set(TAXONOMIES)
    for tax in TAXONOMIES:
        assert os.path.exists(written[tax]["per_genome"])
        assert os.path.isdir(written[tax]["per_rank_dir"])
        assert written[tax]["per_read"].endswith("truth-per-read.tsv.gz")
        assert os.path.exists(written[tax]["per_read"])

    # gtdb tree's euk genus is NA; ncbi tree's euk genus is real
    g = pd.read_csv(written["gtdb"]["per_genome"], sep="\t", keep_default_na=False)
    n = pd.read_csv(written["ncbi"]["per_genome"], sep="\t", keep_default_na=False)
    assert g.set_index("accession").loc["GCA_E", "genus"] == "NA"
    assert n.set_index("accession").loc["GCA_E", "genus"] == "Saccharomyces"

    # per-read is gzipped and reads back
    rt = pd.read_csv(written["ncbi"]["per_read"], sep="\t")
    assert len(rt) == 2
