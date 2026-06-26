import os
import pandas as pd # type: ignore
import pytest # type: ignore

from bit.modules.gen_mg.truth import (
    add_realized_columns,
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
    a 3-genome table as it looks after add_realized_columns: carrying BOTH
    taxonomies, a taxid column, and the realized columns (rel_abundance,
    mean_coverage, reads_generated, detection). GCF_A / GCF_B are GTDB-resolved
    (gtdb_* and ncbi_* both populated); GCA_E is a euk-style user genome with
    gtdb_* = NA but ncbi_* populated.
    """
    return pd.DataFrame([
        dict(accession="GCF_A", taxid="111", user_supplied=False,
             **_ranks("gtdb", ["Bacteria", "P", "C", "O", "F", "Genus_A", "sp A"]),
             **_ranks("ncbi", ["Bacteria", "Pn", "Cn", "On", "Fn", "Genus_A", "sp A"]),
             downloaded_genome_size=1_000_000, used_genome_size=1_000_000,
             rel_abundance=0.5, mean_coverage=10.0, reads_generated=5000,
             detection=0.999),
        dict(accession="GCF_B", taxid="222", user_supplied=False,
             **_ranks("gtdb", ["Bacteria", "P", "C", "O", "F", "Genus_B", "sp B"]),
             **_ranks("ncbi", ["Bacteria", "Pn", "Cn", "On", "Fn", "Genus_B", "sp B"]),
             downloaded_genome_size=500_000, used_genome_size=500_000,
             rel_abundance=0.3, mean_coverage=12.0, reads_generated=3000,
             detection=0.998),
        dict(accession="GCA_E", taxid="NA", user_supplied=True,
             **_ranks("gtdb", ["NA", "NA", "NA", "NA", "NA", "NA", "NA"]),
             **_ranks("ncbi", ["Eukaryota", "Ascomycota", "Saccharomycetes",
                               "Saccharomycetales", "Saccharomycetaceae",
                               "Saccharomyces", "Saccharomyces cerevisiae"]),
             downloaded_genome_size=300_000, used_genome_size=300_000,
             rel_abundance=0.2, mean_coverage=8.0, reads_generated=2000,
             detection=0.997),
    ])


# ─── gen-reads args ────────────────────────────────────────────────────────

def test_build_gen_reads_args_defaults():
    args = build_gen_reads_args(["a.fasta", "b.fasta"], "cov.tsv", "pref")
    assert args.coverage == "cov.tsv"          # coverage-specified mode
    assert args.per_read_tsv is False          # intermediate; off unless per-read wanted
    assert args.read_length == 150
    assert args.type == "paired-end"


def test_build_gen_reads_args_per_read_tsv_opt_in():
    args = build_gen_reads_args(["a.fasta"], "cov.tsv", "pref", per_read_tsv=True)
    assert args.per_read_tsv is True


def test_build_gen_reads_args_long_default_length():
    args = build_gen_reads_args(["a.fasta"], "cov.tsv", "pref", read_type="long")
    assert args.read_length == 5000


# ─── realized columns from gen-reads stats ─────────────────────────────────

def test_add_realized_columns_basic_math():
    merged = pd.DataFrame({"accession": ["A", "B"],
                           "used_genome_size": [1000, 2000]})
    stats = {
        "A": {"detection": 0.9999988, "genome_size": 1000,
              "reads_generated": 300, "realized_bases": 45000},
        "B": {"detection": 0.5004, "genome_size": 2000,
              "reads_generated": 100, "realized_bases": 15000},
    }
    out = add_realized_columns(merged, stats)
    a = out[out["accession"] == "A"].iloc[0]
    assert a["reads_generated"] == 300
    assert a["mean_coverage"] == 45.0          # 45000 / 1000
    assert a["rel_abundance"] == 0.75          # 300 / 400
    assert a["detection"] == 1.0               # 0.9999988 -> 1.0 at 3 places
    b = out[out["accession"] == "B"].iloc[0]
    assert b["mean_coverage"] == 7.5           # 15000 / 2000
    assert b["rel_abundance"] == 0.25
    assert b["detection"] == 0.5               # 0.5004 -> 0.500 at 3 places
    # realized rel_abundance sums to 1 across the community
    assert abs(out["rel_abundance"].sum() - 1.0) < 1e-9


def test_add_realized_columns_missing_genome_gets_zeros():
    # a genome absent from stats (e.g. NA-coverage, produced no reads) -> zeros
    merged = pd.DataFrame({"accession": ["A", "C"],
                           "used_genome_size": [1000, 500]})
    stats = {"A": {"detection": 1.0, "genome_size": 1000,
                   "reads_generated": 200, "realized_bases": 30000}}
    out = add_realized_columns(merged, stats)
    c = out[out["accession"] == "C"].iloc[0]
    assert c["reads_generated"] == 0
    assert c["mean_coverage"] == 0.0
    assert c["rel_abundance"] == 0.0
    assert c["detection"] == 0.0


def test_add_realized_columns_zero_total_reads_no_div_by_zero():
    merged = pd.DataFrame({"accession": ["A"], "used_genome_size": [1000]})
    out = add_realized_columns(merged, {})        # no stats at all
    assert out.iloc[0]["rel_abundance"] == 0.0
    assert out.iloc[0]["mean_coverage"] == 0.0


def test_add_realized_columns_zero_genome_size_no_div_by_zero():
    merged = pd.DataFrame({"accession": ["A"], "used_genome_size": [0]})
    stats = {"A": {"detection": 0.0, "genome_size": 0,
                   "reads_generated": 0, "realized_bases": 0}}
    out = add_realized_columns(merged, stats)
    assert out.iloc[0]["mean_coverage"] == 0.0


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
    assert pg.columns[1] == "taxid"          # taxid sits right after accession
    for r in RANKS:
        assert r in pg.columns               # plain rank names
    assert "gtdb_genus" not in pg.columns    # source prefix stripped
    assert "ncbi_genus" not in pg.columns    # other taxonomy absent
    assert "source_db" not in pg.columns     # dropped from schema
    assert "metadata_genome_size" not in pg.columns
    assert "downloaded_genome_size" in pg.columns
    assert "used_genome_size" in pg.columns
    # realized columns present under their plain names; assigned_* targets absent
    for c in ("rel_abundance", "mean_coverage", "reads_generated", "detection"):
        assert c in pg.columns
    assert "assigned_rel_abundance" not in pg.columns
    assert "assigned_coverage" not in pg.columns
    assert "assigned_reads" not in pg.columns


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
    extra["rel_abundance"] = 0.1
    merged.loc[0, "rel_abundance"] = 0.4
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
    build_read_truth(str(src), fasta_to_acc, pg, str(out))
    rt = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)

    assert (rt["accession"] != "NA").all()
    # GCF_A reads carry Genus_A (gtdb); euk read carries NA under gtdb
    assert (rt[rt["accession"] == "GCF_A"]["genus"] == "Genus_A").all()
    assert (rt[rt["accession"] == "GCA_E"]["genus"] == "NA").all()
    # coordinates carried through, and appear before the rank columns
    assert set(["contig", "start", "end", "strand"]).issubset(rt.columns)
    # 'wrapped' is intentionally dropped from the truth table (gen-metagenome
    # never circularizes, so it would be uniformly 'false')
    assert "wrapped" not in rt.columns
    cols = list(rt.columns)
    assert cols.index("strand") < cols.index("genus")
    assert out.exists()


def test_read_truth_ncbi_taxonomy(merged, tmp_path):
    """ the euk read gets real NCBI ranks when the ncbi per-genome table is used. """
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "ncbi")
    src = tmp_path / "src.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\t/x/GCA_E.fasta\tGCA_E_ctg0\t5\t155\t+\tfalse\n"
    )
    out = tmp_path / "o.tsv"
    build_read_truth(str(src), {"/x/GCA_E.fasta": "GCA_E", "GCA_E.fasta": "GCA_E"},
                     pg, str(out))
    rt = pd.read_csv(out, sep="\t", dtype=str)
    assert rt["genus"].iloc[0] == "Saccharomyces"


def test_read_truth_basename_fallback(merged, tmp_path):
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    src = tmp_path / "src.tsv"
    src.write_text(
        "read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n"
        "r1/1\tGCF_B.fasta\tGCF_B_ctg0\t1\t151\t+\tfalse\n"
    )
    out = tmp_path / "o.tsv"
    build_read_truth(str(src), {"GCF_B.fasta": "GCF_B"}, pg, str(out))
    rt = pd.read_csv(out, sep="\t", dtype=str)
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


def test_build_read_truth_streams_large_input(merged, tmp_path):
    """ a larger input streams through polars to a single gzipped output with the
    header once and no rows lost/duplicated (memory is bounded by the engine). """
    import gzip
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    n = 5000
    lines = ["read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped"]
    for i in range(n):
        lines.append(f"r{i}/1\t/x/GCF_A.fasta\tGCF_A_ctg0\t{i}\t{i+150}\t+\tfalse")
    src = tmp_path / "src.tsv"
    src.write_text("\n".join(lines) + "\n")
    f2a = {"/x/GCF_A.fasta": "GCF_A", "GCF_A.fasta": "GCF_A"}

    out = tmp_path / "out.tsv.gz"
    build_read_truth(str(src), f2a, pg, str(out))

    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert len(df) == n                                  # no rows lost
    assert (df["accession"] == "GCF_A").all()            # all mapped
    assert (df["genus"] == "Genus_A").all()              # taxonomy joined
    with gzip.open(out, "rt") as fh:                      # header exactly once
        assert sum(1 for ln in fh if ln.startswith("read_id\t")) == 1


def test_build_read_truth_progress_callback_per_chunk(merged, tmp_path):
    """ the progress callback fires once per chunk, and the chunked write stays
    complete and correct (memory bounded by chunksize). """
    pg = build_per_genome_table(merged, {a: {"rate": 0.0} for a in merged["accession"]}, "gtdb")
    n = 500
    lines = ["read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped"]
    for i in range(n):
        lines.append(f"r{i}/1\t/x/GCF_A.fasta\tGCF_A_ctg0\t{i}\t{i+150}\t+\tfalse")
    src = tmp_path / "src.tsv"
    src.write_text("\n".join(lines) + "\n")

    calls = []
    out = tmp_path / "o.tsv.gz"
    build_read_truth(str(src), {"/x/GCF_A.fasta": "GCF_A", "GCF_A.fasta": "GCF_A"},
                     pg, str(out), chunksize=200, progress=lambda: calls.append(1))
    assert len(calls) == 3                              # ceil(500/200)
    df = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert len(df) == n
    assert (df["genus"] == "Genus_A").all()
