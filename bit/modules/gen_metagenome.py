"""
builds mock metagenomes with ground-truth tables by composing existing bit
machinery: GTDB/NCBI selection, dl-ncbi-assemblies, mutate-seqs, and gen-reads

Phases:
  1. select    - GTDB selection and/or user accessions, normalized + merged into one genome table
  2. download  - dl-ncbi-assemblies
  3. mutate    - optional per-genome mutation
  4. abundance - assign rel-abundance/coverage, write coverage TSV
  5. reads     - gen-reads
  6. truth     - per-genome, per-rank, and optional per-read truth tables
"""

import os
from types import SimpleNamespace

import pandas as pd
from tqdm import tqdm  # type: ignore

from bit.modules.general import (color_text, report_message,
                                 attempt_to_make_dir, check_files_are_found)
from bit.modules.gtdb.get_gtdb_data import get_gtdb_data
from bit.modules.ncbi.dl_ncbi_assemblies import dl_ncbi_assemblies

from bit.modules.gen_mg import selection as SEL
from bit.modules.gen_mg import abundance as ABD
from bit.modules.gen_mg import mutation as MUT
from bit.modules.gen_mg import truth as TRU


def gen_metagenome(args):

    run = setup(args)

    section("Selecting genomes")
    run = phase_select(args, run)

    section("Downloading assemblies")
    run = phase_download(args, run)

    if args.mutation_mode != "off":
        section("Mutating genomes")
    run = phase_mutate(args, run)

    section("Assigning abundances")
    run = phase_abundance(args, run)

    section("Generating reads")
    run = phase_reads(args, run)

    section("Building truth tables")
    phase_truth(args, run)

    report_finish(run)


# ----------------------------------------------------------------------------

def section(title):
    print(color_text(f"\n  {title}", "yellow"))


def warn_all(warnings):
    for w in warnings:
        report_message(w, "orange", initial_indent="    ", subsequent_indent="    ")


def setup(args):
    if args.output_dir and not os.path.exists(args.output_dir):
        attempt_to_make_dir(args.output_dir)
    genomes_dir = os.path.join(args.output_dir, "genomes")
    attempt_to_make_dir(genomes_dir)
    return SimpleNamespace(
        out_dir=args.output_dir,
        genomes_dir=genomes_dir,
        merged=None,
        genome_paths={},      # accession -> downloaded fasta path
        working_paths={},     # accession -> fasta to read from (mutated if applicable)
        mutation_results={},
        coverage_tsv=os.path.join(args.output_dir, "per-genome-coverage.tsv"),
        reads_prefix=os.path.join(args.output_dir, args.output_prefix),
    )


# ---- phase 1: select ----

def phase_select(args, run):

    user_df = None
    if args.accessions:
        check_files_are_found([args.accessions])
        user_df = load_user_accessions(args)

    generative_df = None
    if args.num_genomes:
        gtdb_dir = get_gtdb_data(quiet=True)
        gtdb_tab = pd.read_csv(
            os.path.join(gtdb_dir, "GTDB-arc-and-bac-metadata.tsv"),
            sep="\t", low_memory=False)

        exclude = SEL.user_filled_groups(user_df, args.derep_rank) if user_df is not None else set()
        sel, sw = SEL._select_gtdb_one_per_rank(
            gtdb_tab, derep_rank=args.derep_rank, domains=args.domains,
            num_genomes=args.num_genomes, seed=args.seed, exclude_groups=exclude)
        warn_all(sw)
        generative_df = SEL._normalize_gtdb_rows(sel)

    merged, mw = SEL.merge_sources(generative_df, user_df)
    warn_all(mw)

    if len(merged) == 0:
        report_message("No genomes selected; nothing to do.", "red",
                       initial_indent="    ", subsequent_indent="    ")
        raise SystemExit(1)

    run.merged = merged
    report_message(f"Selected {len(merged)} genome(s).", "green",
                   initial_indent="    ", subsequent_indent="    ")
    return run


def load_user_accessions(args):
    """
    read user accession file. Bare one-per-line, or a TSV whose first column is
    the accession and optional named columns pin rel_abundance / coverage /
    mutation_rate. NCBI metadata (sizes, taxonomy) is fetched for these
    """
    first = open(args.accessions).readline()
    has_header = "accession" in first.lower()
    if has_header:
        udf = pd.read_csv(args.accessions, sep="\t", dtype=str)
        accs = udf["accession"].tolist()
        pins = {a: {} for a in accs}
        for col, key in (("rel_abundance", "pinned_rel_abundance"),
                         ("coverage", "pinned_coverage"),
                         ("mutation_rate", "pinned_mutation_rate")):
            if col in udf.columns:
                for a, v in zip(udf["accession"], udf[col]):
                    pins[a][key] = v
    else:
        accs = [l.strip() for l in open(args.accessions) if l.strip()]
        pins = {a: {} for a in accs}

    # fetch NCBI metadata (sizes, source) for these accessions
    ncbi_tab = fetch_ncbi_metadata_for_accessions(accs)
    udf = SEL._normalize_ncbi_rows(ncbi_tab, lineage_lookup=None, user_supplied=True)

    # apply pins
    for key in ("pinned_rel_abundance", "pinned_coverage", "pinned_mutation_rate"):
        udf[key] = [
            _num_or_na(pins.get(a, {}).get(key)) for a in udf["accession"]
        ]
    return udf


def _num_or_na(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return pd.NA


def fetch_ncbi_metadata_for_accessions(accessions):
    """
    query NCBI Datasets per-accession
    Returns a frame with at least: accession, total_sequence_length, source_database.
    """
    import subprocess
    from io import StringIO
    fields = ["accession", "assmstats-total-sequence-len", "source_database"]
    summary = subprocess.Popen(
        ["datasets", "summary", "genome", "accession", *accessions, "--as-json-lines"],
        stdout=subprocess.PIPE)
    fmt = subprocess.Popen(
        ["dataformat", "tsv", "genome", "--fields", ",".join(fields)],
        stdin=summary.stdout, stdout=subprocess.PIPE)
    summary.stdout.close()
    out, _ = fmt.communicate()
    tab = pd.read_csv(StringIO(out.decode()), sep="\t", dtype=str)
    tab.columns = ["accession", "total_sequence_length", "source_database"]
    src_map = {"SOURCE_DATABASE_GENBANK": "genbank", "SOURCE_DATABASE_REFSEQ": "refseq"}
    tab["source_database"] = tab["source_database"].map(lambda v: src_map.get(v, v))
    return tab


# ---- phase 2: download ----

def phase_download(args, run):
    acc_file = os.path.join(run.out_dir, "selected-accessions.txt")
    with open(acc_file, "w") as fh:
        for a in run.merged["accession"]:
            fh.write(a + "\n")

    dl_args = SimpleNamespace(
        wanted_accessions=acc_file, format="fasta",
        jobs=args.jobs, output_dir=run.genomes_dir)
    dl_ncbi_assemblies(dl_args)

    # resolve downloaded fasta paths per accession
    for a in run.merged["accession"]:
        run.genome_paths[a] = find_downloaded_fasta(run.genomes_dir, a)
    run.working_paths = dict(run.genome_paths)
    return run


def find_downloaded_fasta(genomes_dir, accession):
    """ locate the downloaded fasta for an accession (dl-ncbi-assemblies naming) """
    for fn in os.listdir(genomes_dir):
        if fn.startswith(accession) and (fn.endswith(".fasta") or fn.endswith(".gz")):
            return os.path.join(genomes_dir, fn)
    return None


# ---- phase 4: abundance (runs after mutate so coverage uses used_genome_size) ----

def phase_abundance(args, run):
    assigned, w = ABD.assign_abundance(
        run.merged, mode=args.abundance_mode, dist=args.abundance_dist,
        sigma=args.sigma, total_reads=args.total_reads, read_length=args.read_length,
        read_type=args.type, fragment_size=args.fragment_size,
        median_coverage=args.median_coverage, seed=args.seed)
    warn_all(w)
    run.merged = assigned
    return run


# ---- phase 3: mutate (runs before abundance so coverage uses used sizes) ----

def measure_fasta_size(path):
    """ sum sequence lengths in a fasta (handles .gz). measured, not metadata. """
    import gzip
    from Bio import SeqIO  # type: ignore
    if path is None or not os.path.exists(path):
        return None
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            total += len(rec.seq)
    return total


def phase_mutate(args, run):
    accs = list(run.merged["accession"])
    rates = MUT.assign_rates(
        accs, mode=args.mutation_mode, mutation_rate=args.mutation_rate,
        rate_min=args.mutation_rate_min, rate_max=args.mutation_rate_max, seed=args.seed)
    if "pinned_mutation_rate" in run.merged.columns:
        for a, p in zip(run.merged["accession"], run.merged["pinned_mutation_rate"]):
            if pd.notna(p):
                rates[a] = float(p)

    if args.mutation_mode == "off" and all(v == 0 for v in rates.values()):
        run.mutation_results = {a: {"rate": 0.0} for a in accs}
        run.working_paths = dict(run.genome_paths)
    else:
        mut_dir = os.path.join(run.out_dir, "mutated-genomes")
        with tqdm(total=len(accs), desc="    Mutating", ncols=70) as pbar:
            run.mutation_results = MUT.run_mutation(
                run.genome_paths, rates, mut_dir, ti_tv_ratio=args.ti_tv_ratio,
                indel_rate=args.indel_rate, seed=args.seed, progress=pbar)
        run.working_paths = {a: run.mutation_results[a]["mutated_fasta"] for a in accs}

    # populate used_genome_size: measured from the fasta reads will come from.
    # mutation agg already carries the mutated size; otherwise measure the download.
    used = []
    for a in accs:
        mr = run.mutation_results.get(a, {})
        size = mr.get("genome_size")             # present only when this genome was mutated
        if size is None:
            size = measure_fasta_size(run.working_paths.get(a))
        used.append(size if size else pd.NA)
    run.merged["used_genome_size"] = used
    return run


# ---- phase 5: reads ----

def phase_reads(args, run):
    from bit.modules.gen_reads import generate_reads

    # coverage tsv keyed by the fasta gen-reads will actually read
    run.merged["_fasta_for_reads"] = run.merged["accession"].map(run.working_paths)
    ABD.write_coverage_tsv(run.merged, "_fasta_for_reads", run.coverage_tsv)

    gr_args = TRU.build_gen_reads_args(
        run.merged["_fasta_for_reads"].dropna().tolist(), run.coverage_tsv,
        run.reads_prefix, read_type=args.type, read_length=args.read_length,
        fragment_size=args.fragment_size, fragment_size_range=args.fragment_size_range,
        long_read_length_range=args.long_read_length_range, seed=args.seed,
        include_Ns=args.include_Ns)
    generate_reads(gr_args)
    run.read_sources_tsv = f"{run.reads_prefix}-read-sources.tsv"
    return run


# ---- phase 6: truth ----

def phase_truth(args, run):
    per_genome = TRU.build_per_genome_table(run.merged, run.mutation_results)
    pg_path = os.path.join(run.out_dir, "truth-per-genome.tsv")
    per_genome.to_csv(pg_path, sep="\t", index=False)

    per_rank = TRU.build_per_rank_tables(per_genome)
    rank_dir = os.path.join(run.out_dir, "truth-per-rank")
    attempt_to_make_dir(rank_dir)
    for rank, tab in per_rank.items():
        tab.to_csv(os.path.join(rank_dir, f"truth-{rank}-abundance.tsv"),
                   sep="\t", index=False)

    run.truth_files = [pg_path, rank_dir]

    if args.per_read_tsv:
        fasta2acc = {}
        for a, p in run.working_paths.items():
            if p:
                fasta2acc[p] = a
                fasta2acc[os.path.basename(p)] = a
        rt_path = os.path.join(run.out_dir, "truth-read-level.tsv")
        TRU.build_read_truth(run.read_sources_tsv, fasta2acc, per_genome, rt_path)
        run.truth_files.append(rt_path)


def report_finish(run):
    report_message("Mock metagenome complete!", "green",
                   initial_indent="    ", subsequent_indent="    ",
                   trailing_newline=True)
    print(f"        Reads:            {run.reads_prefix}_R1.fastq.gz, _R2.fastq.gz")
    print(f"        Per-genome truth: {os.path.join(run.out_dir, 'truth-per-genome.tsv')}")
    print(f"        Per-rank truth:   {os.path.join(run.out_dir, 'truth-per-rank/')}")
    print("")
