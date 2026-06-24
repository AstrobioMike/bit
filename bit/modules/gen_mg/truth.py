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
import gzip
import pandas as pd # type: ignore
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
                         include_Ns=False, genome_sizes=None, jobs=10, per_read_tsv=False):
    """
    Construct an args-like object matching what bit.modules.gen_reads.generate_reads
    expects. Coverage mode is used (the coverage TSV drives per-genome read counts).

    per_read_tsv: whether gen-reads should emit its per-read source TSV. This is an
    intermediate consumed by build_read_truth to make the per-read truth tables,
    so gen-metagenome only needs it when --per-read-tsv is set.

    genome_sizes: optional {fasta_path: measured_size} so gen-reads can skip
    re-measuring genomes gen-metagenome already sized
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
        per_read_tsv=per_read_tsv,
        genome_sizes=genome_sizes or {},
        quiet=True,
        jobs=jobs
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

def build_read_truth(read_sources_tsv, fasta_to_accession, per_genome_df, out_path,
                     chunksize=500_000, progress=None):
    """
    Join gen-reads' per-read source TSV to accession + taxonomy, writing the
    read-level truth table. per_genome_df is expected to carry plain rank columns
    for a single taxonomy (so call once per source with that source's table).

    The source TSV has one row per read and can be enormous (hundreds of millions
    of rows for high read counts), so it is processed in fixed-size chunks: each
    chunk is read, mapped to accession + taxonomy, and appended to the gzipped
    output, so peak memory is bounded by chunksize rather than the total read
    count. The gzip is streamed directly (no uncompressed intermediate on disk).
    Returns the output path.

    progress: optional callable invoked once per chunk written, so a caller can
    drive a progress bar without this module doing console I/O.

    fasta_to_accession maps source_fasta -> accession (the caller populates both
    the full-path and basename forms); anything unmatched becomes "NA". Unmatched
    ranks are written as the literal string "NA" to match the per-genome table's
    own NA convention.
    """
    def resolve(sf):
        if sf in fasta_to_accession:
            return fasta_to_accession[sf]
        return fasta_to_accession.get(os.path.basename(sf), "NA")

    tax = per_genome_df.set_index("accession")
    rank_cols = [r for r in RANKS if r in tax.columns]
    # 'wrapped' is intentionally omitted from the truth table: gen-metagenome
    # never circularizes (gen-reads' circularize stays False), so it is uniformly
    # 'false' and carries no information here. It still exists in gen-reads' own
    # output for standalone circularized runs. We skip it at parse time (usecols)
    # rather than reading then dropping, so it is never parsed on any chunk.
    base_cols = ["read_id", "accession", "contig", "start", "end", "strand"]

    # columns to actually read from the source: source_fasta (to derive accession)
    # plus the non-derived base columns. Intersect with the header so a future
    # gen-reads format change degrades gracefully instead of erroring on usecols.
    want_cols = ["read_id", "source_fasta", "contig", "start", "end", "strand"]
    present = pd.read_csv(read_sources_tsv, sep="\t", nrows=0).columns
    read_cols = [c for c in want_cols if c in present]

    out_path = str(out_path)
    gzipped = out_path.endswith(".gz")
    opener = gzip.open if gzipped else open

    wrote_header = False
    with opener(out_path, "wt", newline="") as out:
        for chunk in pd.read_csv(read_sources_tsv, sep="\t", dtype=str,
                                 chunksize=chunksize, usecols=read_cols):
            chunk["accession"] = chunk["source_fasta"].map(resolve)
            joined = chunk.merge(tax[rank_cols], left_on="accession",
                                 right_index=True, how="left")
            col_order = [c for c in (base_cols + rank_cols) if c in joined.columns]
            joined[col_order].to_csv(out, sep="\t", index=False,
                                     header=not wrote_header)
            wrote_header = True
            if progress is not None:
                progress()

    return out_path


# ---------- top-level: write the split ground-truth tree ----------

def build_truth_for_taxonomy(taxonomy, merged_df, mutation_results, gt_root,
                             read_sources_tsv=None, fasta_to_accession=None,
                             per_read=False, attempt_to_make_dir=os.makedirs,
                             on_per_read=None, progress=None, chunksize=500_000):
    """
    Build the ground-truth/<taxonomy>/ subtree: per-genome table, per-rank
    abundance tables, and (optionally) the per-read truth table.

    on_per_read: optional zero-arg callable invoked just before the (slow)
    per-read build, so a caller can show progress feedback around it without
    this module doing any console I/O.

    Returns {per_genome, per_rank_dir, per_read?}.
    """
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
        if on_per_read is not None:
            on_per_read()
        rt_path = os.path.join(tax_dir, "truth-per-read.tsv.gz")
        build_read_truth(read_sources_tsv, fasta_to_accession, per_genome, rt_path,
                         chunksize=chunksize, progress=progress)
        entry["per_read"] = rt_path

    return entry


def write_truth_outputs(merged_df, mutation_results, out_dir,
                        read_sources_tsv=None, fasta_to_accession=None,
                        per_read=False, attempt_to_make_dir=os.makedirs):
    """
    Write the full ground-truth/ tree, one subdirectory per taxonomy source.

    Returns a dict taxonomy -> {per_genome, per_rank_dir, per_read?} paths.
    """
    written = {}
    gt_root = os.path.join(out_dir, "ground-truth")
    _mkdir(attempt_to_make_dir, gt_root)

    for taxonomy in TAXONOMIES:
        written[taxonomy] = build_truth_for_taxonomy(
            taxonomy, merged_df, mutation_results, gt_root,
            read_sources_tsv=read_sources_tsv, fasta_to_accession=fasta_to_accession,
            per_read=per_read, attempt_to_make_dir=attempt_to_make_dir)

    return written


def _mkdir(maker, path):
    """ make a directory tolerantly whether maker is os.makedirs or bit's helper. """
    try:
        maker(path, exist_ok=True)
    except TypeError:
        if not os.path.exists(path):
            maker(path)
