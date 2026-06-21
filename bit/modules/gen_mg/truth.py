"""
gen-metagenome reads + truth-table assembly.

Calls bit's gen-reads (in-process) across all genome fastas using a coverage TSV,
with per-read provenance on, then assembles the ground-truth outputs.

Truth outputs are written per taxonomy source into separate directories so a user
benchmarking against GTDB (or NCBI) taxonomy can point at one tree:

    ground-truth/
      gtdb/  truth-per-genome.tsv, truth-per-read.tsv.gz, truth-per-rank/<rank>.tsv
      ncbi/  truth-per-genome.tsv, truth-per-read.tsv.gz, truth-per-rank/<rank>.tsv

The merged table carries both gtdb_* and ncbi_* rank columns; each builder takes a
`taxonomy` argument ('gtdb' or 'ncbi'), selects that source's rank columns, and
emits them under plain rank names so each tree is self-consistent.
"""
import os
import pandas as pd
from types import SimpleNamespace

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
TAXONOMIES = ["gtdb", "ncbi"]


def _src_rank_cols(taxonomy):
    """ prefixed source columns (e.g. gtdb_domain..) for a taxonomy. """
    return [f"{taxonomy}_{r}" for r in RANKS]


# ---------- gen-reads invocation ----------

def build_gen_reads_args(fasta_paths, coverage_tsv, output_prefix, read_type="paired-end",
                         read_length=None, fragment_size=500, fragment_size_range=10,
                         long_read_length_range=50, seed=None, circularize=False,
                         include_Ns=False):
    """
    Construct an args-like object matching what bit.modules.gen_reads.generate_reads
    expects. Coverage mode is used (the coverage TSV drives per-genome read counts),
    and --source-tsv provenance is always on for gen-metagenome.
    """
    if read_length is None:
        read_length = 5000 if read_type == "long" else 150
    return SimpleNamespace(
        input_fastas=list(fasta_paths),
        output_prefix=output_prefix,
        type=read_type,
        num_reads=1_000_000,          # overridden internally by coverage mode
        read_length=read_length,
        coverage=coverage_tsv,        # path -> coverage-specified mode
        proportions_file=None,
        circularize=circularize,
        include_Ns=include_Ns,
        seed=seed,
        fragment_size=fragment_size,
        fragment_size_range=fragment_size_range,
        long_read_length_range=long_read_length_range,
        source_tsv=True,
        quiet=True,
    )


# ---------- per-genome truth ----------

def build_per_genome_table(merged_df, mutation_results, taxonomy):
    """
    Per-genome truth for one taxonomy source ('gtdb' or 'ncbi'). The chosen
    source's rank columns are emitted under plain rank names (domain..species).

    merged_df: normalized + abundance-assigned genome table (one row per genome),
               carrying gtdb_* and ncbi_* rank columns.
    mutation_results: dict accession -> {rate, num_substitutions, ...}.
    """
    if taxonomy not in TAXONOMIES:
        raise ValueError(f"taxonomy must be one of {TAXONOMIES}, got '{taxonomy}'")

    out = merged_df.copy()
    rates, subs, indels, totch = [], [], [], []
    for acc in out["accession"]:
        mr = mutation_results.get(acc, {})
        rates.append(mr.get("rate", 0.0))
        subs.append(mr.get("num_substitutions", 0))
        indels.append(mr.get("num_indels", 0))
        totch.append(mr.get("num_total_changes", 0))
    out["mutation_rate"] = rates
    out["num_substitutions"] = subs
    out["num_indels"] = indels
    out["num_total_changes"] = totch

    # rename this taxonomy's source columns to plain rank names
    src_cols = _src_rank_cols(taxonomy)
    for src, plain in zip(src_cols, RANKS):
        if src in out.columns:
            out[plain] = out[src]

    cols = (["accession"] + RANKS +
            ["downloaded_genome_size", "used_genome_size",
             "assigned_rel_abundance", "assigned_coverage",
             "assigned_reads", "mutation_rate", "num_substitutions", "num_indels",
             "num_total_changes", "user_supplied"])
    cols = [c for c in cols if c in out.columns]
    return out[cols]


# ---------- per-rank collapsed abundance (profiler truth) ----------

def build_per_rank_tables(per_genome_df):
    """
    Collapse assigned_rel_abundance to each rank. Expects plain rank columns
    (as produced by build_per_genome_table for a given taxonomy). Returns
    dict rank -> dataframe [<rank>, rel_abundance]. 'NA' genomes group under 'NA'.
    """
    tables = {}
    for rank in RANKS:
        if rank not in per_genome_df.columns:
            continue
        grp = (per_genome_df.groupby(rank, dropna=False)["assigned_rel_abundance"]
               .sum().reset_index()
               .rename(columns={"assigned_rel_abundance": "rel_abundance"})
               .sort_values("rel_abundance", ascending=False)
               .reset_index(drop=True))
        tables[rank] = grp
    return tables


# ---------- read-level truth ----------

def build_read_truth(read_sources_tsv, fasta_to_accession, per_genome_df, out_path):
    """
    Join gen-reads' per-read source TSV to accession + taxonomy, writing the
    read-level truth table. per_genome_df is expected to carry plain rank columns
    for a single taxonomy (so call once per source with that source's table).
    """
    src = pd.read_csv(read_sources_tsv, sep="\t", dtype=str)

    def resolve(sf):
        if sf in fasta_to_accession:
            return fasta_to_accession[sf]
        return fasta_to_accession.get(os.path.basename(sf), "NA")
    src["accession"] = src["source_fasta"].map(resolve)

    tax = per_genome_df.set_index("accession")
    rank_cols = [r for r in RANKS if r in tax.columns]
    joined = src.merge(tax[rank_cols], left_on="accession", right_index=True, how="left")

    col_order = (["read_id", "accession", "contig", "start", "end", "strand", "wrapped"] + rank_cols)
    col_order = [c for c in col_order if c in joined.columns]
    joined = joined[col_order]
    compression = "gzip" if str(out_path).endswith(".gz") else None
    joined.to_csv(out_path, sep="\t", index=False, compression=compression)
    return joined


# ---------- top-level: write the split ground-truth tree ----------

def write_truth_outputs(merged_df, mutation_results, out_dir,
                        read_sources_tsv=None, fasta_to_accession=None,
                        per_read=False, attempt_to_make_dir=os.makedirs):
    """
    Write the full ground-truth/ tree, one subdirectory per taxonomy source.

    Returns a dict taxonomy -> {per_genome, per_rank_dir, per_read?} paths.
    """
    written = {}
    gt_root = os.path.join(out_dir, "ground-truth")
    attempt_to_make_dir(gt_root, exist_ok=True) if attempt_to_make_dir is os.makedirs else attempt_to_make_dir(gt_root)

    for taxonomy in TAXONOMIES:
        tax_dir = os.path.join(gt_root, taxonomy)
        _mkdir(attempt_to_make_dir, tax_dir)

        per_genome = build_per_genome_table(merged_df, mutation_results, taxonomy)
        pg_path = os.path.join(tax_dir, "truth-per-genome.tsv")
        per_genome.to_csv(pg_path, sep="\t", index=False)
        entry = {"per_genome": pg_path}

        rank_dir = os.path.join(tax_dir, "truth-per-rank")
        _mkdir(attempt_to_make_dir, rank_dir)
        for rank, tab in build_per_rank_tables(per_genome).items():
            tab.to_csv(os.path.join(rank_dir, f"truth-{rank}-abundance.tsv"),
                       sep="\t", index=False)
        entry["per_rank_dir"] = rank_dir

        if per_read and read_sources_tsv and fasta_to_accession is not None:
            rt_path = os.path.join(tax_dir, "truth-per-read.tsv.gz")
            build_read_truth(read_sources_tsv, fasta_to_accession, per_genome, rt_path)
            entry["per_read"] = rt_path

        written[taxonomy] = entry

    return written


def _mkdir(maker, path):
    """ make a directory tolerantly whether maker is os.makedirs or bit's helper. """
    try:
        maker(path, exist_ok=True)
    except TypeError:
        if not os.path.exists(path):
            maker(path)
