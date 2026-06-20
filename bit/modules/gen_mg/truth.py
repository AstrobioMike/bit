"""
gen-metagenome reads + truth-table assembly.

Calls bit's gen-reads (in-process) across all genome fastas using a coverage TSV,
with per-read provenance on, then assembles the ground-truth outputs:
  - per-genome truth table
  - per-rank collapsed relative-abundance tables (profiler truth)
  - optional read-level truth (read_id -> accession + taxonomy + coords)
"""
import os
import pandas as pd
from types import SimpleNamespace

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


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
    )


# ---------- per-genome truth ----------

def build_per_genome_table(merged_df, mutation_results):
    """
    merged_df: normalized + abundance-assigned genome table (one row per genome).
    mutation_results: dict accession -> {rate, num_substitutions, ...} (or {rate:0.0}).
    Returns a per-genome truth dataframe.
    """
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

    cols = (["accession", "source_db"] + RANKS +
            ["metadata_genome_size", "used_genome_size",
             "assigned_rel_abundance", "assigned_coverage",
             "assigned_reads", "mutation_rate", "num_substitutions", "num_indels",
             "num_total_changes", "taxonomy_source", "user_supplied"])
    cols = [c for c in cols if c in out.columns]
    return out[cols]


# ---------- per-rank collapsed abundance (profiler truth) ----------

def build_per_rank_tables(per_genome_df):
    """
    Collapse assigned_rel_abundance to each rank. Returns dict rank -> dataframe
    with columns [<rank>, rel_abundance], summing genomes that share that taxon.
    Genomes with 'NA' at a rank are grouped under 'NA' (e.g. unresolved euks).
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
    Join gen-reads' per-read source TSV (read_id, source_fasta, contig, start, end,
    strand, wrapped) to accession + taxonomy, writing the read-level truth table.

    fasta_to_accession: dict mapping the fasta path/basename used in gen-reads to
    its accession (so source_fasta -> accession).
    """
    src = pd.read_csv(read_sources_tsv, sep="\t", dtype=str)

    # map source_fasta -> accession (try full path then basename)
    def resolve(sf):
        if sf in fasta_to_accession:
            return fasta_to_accession[sf]
        return fasta_to_accession.get(os.path.basename(sf), "NA")
    src["accession"] = src["source_fasta"].map(resolve)

    tax = per_genome_df.set_index("accession")
    rank_cols = [r for r in RANKS if r in tax.columns]
    joined = src.merge(tax[rank_cols], left_on="accession", right_index=True, how="left")

    col_order = (["read_id", "accession"] + rank_cols +
                 ["contig", "start", "end", "strand", "wrapped"])
    col_order = [c for c in col_order if c in joined.columns]
    joined = joined[col_order]
    joined.to_csv(out_path, sep="\t", index=False)
    return joined
