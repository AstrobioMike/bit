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
                                 attempt_to_make_dir, check_files_are_found, spinner,
                                 log_command_run)
from bit.modules.gtdb.get_gtdb_data import get_gtdb_data
from bit.modules.ncbi.dl_ncbi_assemblies import dl_ncbi_assemblies
from bit.modules.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
from bit.modules.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data
from bit.modules.gen_mg import selection as SEL
from bit.modules.gen_mg import abundance as ABD
from bit.modules.gen_mg import mutation as MUT
from bit.modules.gen_mg import truth as TRU
from bit.modules.gen_mg import taxonomy as TAX


def gen_metagenome(args):

    run = setup(args)
    check_required_dbs(args, run)

    mutating = args.mutation_mode != "off"
    n = _phase_counter()

    section(f"Phase {n()}: Selecting genomes and resolving taxonomy...")
    run = phase_select(args, run)

    section(f"Phase {n()}: Downloading assemblies...")
    run = phase_download(args, run)

    if mutating:
        section(f"Phase {n()}: Mutating genomes...")
    run = phase_mutate(args, run)

    # when not mutating, genome sizes are measured here (needed for abundance)
    if mutating:
        section(f"Phase {n()}: Assigning abundances...")
    else:
        section(f"Phase {n()}: Getting genome sizes and assigning abundances...")
    run = phase_abundance(args, run)

    section(f"Phase {n()}: Generating reads...")
    run = phase_reads(args, run)

    section(f"Phase {n()}: Building truth tables...")
    phase_truth(args, run)

    report_finish(args, run)


def _phase_counter():
    """ returns a callable yielding 1, 2, 3, ... on each call (for phase labels). """
    state = {"i": 0}
    def nxt():
        state["i"] += 1
        return state["i"]
    return nxt


def check_required_dbs(args, run):
    need_gtdb = True
    need_ncbi_tax = True
    # assembly summary is only consulted for accessions not in GTDB; we can't
    # know that until selection, but only user-supplied accessions can ever need
    # it, so prematurely fetching it up front only when --accessions is given
    need_assembly_summary = bool(getattr(args, "accessions", None))
    if need_gtdb:
        run.gtdb_dir = get_gtdb_data(quiet=True)
    if need_ncbi_tax:
        get_ncbi_tax_data(quiet=True)
        log_data_source(run, "NCBI taxonomy (taxdump) retrieved",
                        _read_retrieved_date(os.environ.get("TAXONKIT_DB")))
    if need_assembly_summary:
        get_ncbi_assembly_data(quiet=True)


# ----------------------------------------------------------------------------

def section(title):
    print(color_text(f"\n  {title}", "yellow"))


def section_border():
    print(color_text("      " + "- " * 34, "yellow"))


def warn_all(warnings):
    for w in warnings:
        report_message(w, "orange", initial_indent="    ", subsequent_indent="    ")


def setup(args):
    if args.output_dir and not os.path.exists(args.output_dir):
        attempt_to_make_dir(args.output_dir)
    genomes_dir = os.path.join(args.output_dir, "genomes")
    attempt_to_make_dir(genomes_dir)

    # start the runlog: bit version + the rendered command (matches other subcommands)
    log_file = os.path.join(args.output_dir, "runlog.txt")
    log_command_run(getattr(args, "full_cmd_executed", "(command not captured)"),
                    args.output_dir, log_file)

    return SimpleNamespace(
        out_dir=args.output_dir,
        genomes_dir=genomes_dir,
        log_file=log_file,
        merged=None,
        genome_paths={},      # accession -> downloaded fasta path
        working_paths={},     # accession -> fasta to read from (mutated if applicable)
        mutation_results={},
        coverage_tsv=os.path.join(args.output_dir, "per-genome-coverage.tsv"),
        reads_prefix=os.path.join(args.output_dir, args.output_prefix),
    )


def log_data_source(run, label, detail):
    """ 
    append a reference-data provenance line to the runlog (e.g. GTDB version,
    NCBI assembly-summary retrieval date, NCBI taxdump date)
    """
    log_file = getattr(run, "log_file", None)
    if not log_file:
        return
    with open(log_file, "a") as fh:
        fh.write(f"{label}: {detail}\n")


def _read_retrieved_date(data_dir):
    """
    read a date-retrieved.txt (stored YYYY,MM,DD) from a data dir and return
    it formatted like 'Jun 20, 2026'. Returns '(date unknown)' if absent
    """
    from datetime import datetime
    if not data_dir:
        return "(date unknown)"
    path = os.path.join(data_dir, "date-retrieved.txt")
    if not os.path.exists(path):
        return "(date unknown)"
    with open(path) as fh:
        line = fh.readline().strip()
    if not line:
        return "(date unknown)"
    try:
        return datetime.strptime(line, "%Y,%m,%d").strftime("%b %d, %Y")
    except ValueError:
        return line.replace(",", "-")

def _read_gtdb_version(location):
    """ 
    read GTDB version + release date from the version-info file (no printing;
    the value goes to the runlog and the console shows only a spinner)
    """
    version_info = []
    with open(os.path.join(location, "GTDB-version-info.txt")) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line)
    gtdb_version = version_info[0]
    gtdb_release_date = version_info[1]
    return gtdb_version, gtdb_release_date


# ---- selection ----

def _load_gtdb_table(args, run):
    """ 
    load the GTDB metadata table once (needed for generative selection and/or
    GTDB-taxonomy resolution of user accessions). Cached on run. Silent — the
    caller wraps this in an appropriately-labeled spinner. GTDB version is
    recorded to the runlog here
    """
    if getattr(run, "gtdb_tab", None) is not None:
        return run.gtdb_tab
    gtdb_dir = get_gtdb_data(quiet=True)
    gtdb_version, gtdb_release_date = _read_gtdb_version(gtdb_dir)
    log_data_source(run, "GTDB version", f"{gtdb_version} ({gtdb_release_date})")
    run.gtdb_tab = pd.read_csv(
        os.path.join(gtdb_dir, "GTDB-arc-and-bac-metadata.tsv"),
        sep="\t", low_memory=False)
    return run.gtdb_tab


def _assembly_info_path():
    """ path to bit's stored NCBI assembly-info table (taxid source), if set. """
    base = os.environ.get("NCBI_assembly_data_dir")
    if not base:
        return None
    p = os.path.join(base, "ncbi-assembly-info.tsv")
    return p if os.path.exists(p) else None


def phase_select(args, run):

    generative = bool(args.num_genomes)

    # GTDB table is needed if selecting generatively OR resolving GTDB taxonomy
    # for user accessions. In generative mode the load is the "pool" of genomes
    # to pick from and gets its own spinner; in accessions-only mode it's loaded
    # silently as part of taxonomy resolution (no pool, no selection step).
    gtdb_tab = None
    if generative:
        with spinner("Loading pool of potential genomes...", "Loading pool of potential genomes...", indent="    "):
            gtdb_tab = _load_gtdb_table(args, run)

    user_df = None
    if args.accessions:
        check_files_are_found([args.accessions])
        user_df = load_user_accessions(args)

    if generative:
        # selection (+ user-accession GTDB fill + merge) under one spinner;
        # warnings are collected and printed AFTER the checkmark so they don't
        # collide with the spinner's line redraw.
        warnings = []
        with spinner("Selecting genomes...", "Selecting genomes...", indent="    "):
            if user_df is not None:
                # resolve GTDB taxonomy for user rows so derep exclusion can use it
                user_df = TAX.fill_gtdb_taxonomy(user_df, gtdb_tab)
            exclude = SEL.user_filled_groups(user_df, args.derep_rank) if user_df is not None else set()

            if args.domain == "both":
                domains = None
            else:
                alias = {"bacteria": "Bacteria", "archaea": "Archaea"}
                domains = [alias[args.domain]]
            sel, sw = SEL._select_gtdb_one_per_rank(
                gtdb_tab, derep_rank=args.derep_rank, domains=domains,
                num_genomes=args.num_genomes, seed=args.seed, exclude_groups=exclude)
            warnings += sw
            generative_df = SEL._normalize_gtdb_rows(sel)

            merged, mw = SEL.merge_sources(generative_df, user_df)
            warnings += mw
        warn_all(warnings)
    else:
        # accessions-only: no pool, no selection — just the user's genomes
        merged, mw = SEL.merge_sources(None, user_df)
        warn_all(mw)

    if len(merged) == 0:
        report_message("No genomes selected; nothing to do.", "red",
                    initial_indent="    ", subsequent_indent="    ")
        raise SystemExit(1)

    # resolve both taxonomies. In accessions-only mode the GTDB table hasn't been
    # loaded yet, so it loads silently inside this spinner (it's in service of
    # taxonomy resolution here, not pooling). version/date provenance -> runlog.
    ai_path = _assembly_info_path()
    if _needs_ncbi_assembly_info(merged) and ai_path:
        log_data_source(run, "NCBI assembly summary retrieved",
                        _read_retrieved_date(os.path.dirname(ai_path)))
    with spinner("Resolving taxonomy...", "Resolving taxonomy...", indent="    "):
        if gtdb_tab is None:
            gtdb_tab = _load_gtdb_table(args, run)
            if user_df is not None:
                merged = TAX.fill_gtdb_taxonomy(merged, gtdb_tab)
        merged = TAX.resolve_all(merged, gtdb_tab, ai_path)

    run.merged = merged
    print()

    return run


def _needs_ncbi_assembly_info(merged):
    """ 
    True if any genome lacks GTDB taxonomy after the early GTDB fill, meaning
    it will fall through to the NCBI assembly-summary (and beyond) for its taxid
    """
    if "gtdb_domain" not in merged.columns:
        return False
    return bool((merged["gtdb_domain"].astype(str) == "NA").any())


def load_user_accessions(args):
    """
    read user accession file. Bare one-per-line, or a TSV whose first column is
    the accession and optional named columns pin rel_abundance / coverage /
    mutation_rate. Taxonomy is resolved downstream (resolve_all) and genome sizes
    are measured from the downloaded FASTAs; no NCBI metadata query is needed.
    """
    first = open(args.accessions).readline()
    has_header = "accession" in first.lower()
    if has_header:
        udf_in = pd.read_csv(args.accessions, sep="\t", dtype=str)
        accs = udf_in["accession"].tolist()
        pins = {a: {} for a in accs}
        for col, key in (("rel_abundance", "pinned_rel_abundance"),
                         ("coverage", "pinned_coverage"),
                         ("mutation_rate", "pinned_mutation_rate")):
            if col in udf_in.columns:
                for a, v in zip(udf_in["accession"], udf_in[col]):
                    pins[a][key] = v
    else:
        accs = [l.strip() for l in open(args.accessions) if l.strip()]
        pins = {a: {} for a in accs}

    # build normalized rows directly (ranks NA here; filled by resolve_all)
    ncbi_tab = pd.DataFrame({"accession": accs})
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


# ---- downloading ----

def phase_download(args, run):
    acc_file = os.path.join(run.out_dir, "selected-accessions.txt")
    with open(acc_file, "w") as fh:
        for a in run.merged["accession"]:
            fh.write(a + "\n")

    print()
    dl_args = SimpleNamespace(
        wanted_accessions=acc_file, format="fasta",
        jobs=args.jobs, output_dir=run.genomes_dir, quiet=True)
    dl_ncbi_assemblies(dl_args)

    # resolve downloaded fasta paths per accession
    for a in run.merged["accession"]:
        run.genome_paths[a] = find_downloaded_fasta(run.genomes_dir, a)
    run.working_paths = dict(run.genome_paths)
    print()
    return run


def find_downloaded_fasta(genomes_dir, accession):
    """ locate the downloaded fasta for an accession (dl-ncbi-assemblies naming) """
    for fn in os.listdir(genomes_dir):
        if fn.startswith(accession) and (fn.endswith(".fasta") or fn.endswith(".gz")):
            return os.path.join(genomes_dir, fn)
    return None


# ---- mutating (runs before abundance so coverage uses used sizes) ----

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
        # sizes are measured in the abundance phase (parallel, with a bar), since
        # nothing has read the FASTAs yet and abundance needs the sizes.
    else:
        mut_dir = os.path.join(run.out_dir, "mutated-genomes")
        with tqdm(total=len(accs), desc="    Mutating", ncols=70) as pbar:
            run.mutation_results = MUT.run_mutation(
                run.genome_paths, rates, mut_dir, ti_tv_ratio=args.ti_tv_ratio,
                indel_rate=args.indel_rate, seed=args.seed, progress=pbar)
        run.working_paths = {a: run.mutation_results[a]["mutated_fasta"] for a in accs}

        # mutation already read every input FASTA, so both sizes come free from
        # its aggregated counts (no extra reads): input length -> downloaded,
        # written length -> used (differ only under indels).
        run.merged["downloaded_genome_size"] = [
            _size_or_na(run.mutation_results.get(a, {}).get("input_genome_size")) for a in accs
        ]
        run.merged["used_genome_size"] = [
            _size_or_na(run.mutation_results.get(a, {}).get("genome_size")) for a in accs
        ]
    return run


def _size_or_na(v):
    return v if (v is not None and not pd.isna(v) and v) else pd.NA


# ---- abundance figuring (runs after mutate so coverage uses used_genome_size) ----

def phase_abundance(args, run):
    # if sizes weren't set during mutation (mutation off), measure them now —
    # abundance needs them. Parallelized across --jobs with a progress bar.
    if "used_genome_size" not in run.merged.columns:
        print()
        sizes = measure_sizes_parallel(
            list(run.merged["accession"]), run.working_paths, jobs=args.jobs)
        run.merged["downloaded_genome_size"] = [
            _size_or_na(sizes.get(a)) for a in run.merged["accession"]]
        # no mutation -> used == downloaded
        run.merged["used_genome_size"] = run.merged["downloaded_genome_size"]

    assigned, w = ABD.assign_abundance(
        run.merged, mode=args.abundance_mode, dist=args.abundance_dist,
        sigma=args.sigma, total_reads=args.total_reads, read_length=args.read_length,
        read_type=args.type, fragment_size=args.fragment_size,
        median_coverage=args.median_coverage, seed=args.seed)
    warn_all(w)
    run.merged = assigned
    print()
    return run


def measure_sizes_parallel(accessions, path_map, jobs=10):
    """ 
    get genome size for each accession's FASTA in parallel, with a bar.
    Returns dict accession -> size (or None on failure)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    sizes = {}
    with ThreadPoolExecutor(max_workers=min(max(int(jobs), 1), 20)) as pool:
        futures = {pool.submit(measure_fasta_size, path_map.get(a)): a
                   for a in accessions}
        with tqdm(total=len(futures), desc="    Progress",
                  unit=" genome", ncols=78) as pbar:
            for fut in as_completed(futures):
                a = futures[fut]
                try:
                    sizes[a] = fut.result()
                except Exception:
                    sizes[a] = None
                pbar.update(1)
    return sizes


# ---- generating reads ----

def phase_reads(args, run):
    from bit.modules.gen_reads import generate_reads

    # coverage tsv keyed by the fasta gen-reads will actually read
    run.merged["_fasta_for_reads"] = run.merged["accession"].map(run.working_paths)
    ABD.write_coverage_tsv(run.merged, "_fasta_for_reads", run.coverage_tsv)

    # sizes already measured (mutate or abundance phase); hand them to gen-reads
    # keyed by the working fasta path so it doesn't re-read every genome to size
    # it. used_genome_size matches the working fasta (post-mutation if mutated).
    genome_sizes = {}
    for path, size in zip(run.merged["_fasta_for_reads"], run.merged["used_genome_size"]):
        if path is not None and not pd.isna(path) and pd.notna(size):
            genome_sizes[path] = int(size)

    gr_args = TRU.build_gen_reads_args(
        run.merged["_fasta_for_reads"].dropna().tolist(), run.coverage_tsv,
        run.reads_prefix, read_type=args.type, read_length=args.read_length,
        fragment_size=args.fragment_size, fragment_size_range=args.fragment_size_range,
        long_read_length_range=args.long_read_length_range, seed=args.seed,
        include_Ns=args.include_Ns, genome_sizes=genome_sizes)
    generate_reads(gr_args)
    run.read_sources_tsv = f"{run.reads_prefix}-read-sources.tsv"
    print()
    return run


# ---- making truth outputs ----

def phase_truth(args, run):
    fasta2acc = None
    if args.per_read_tsv:
        fasta2acc = {}
        for a, p in run.working_paths.items():
            if p:
                fasta2acc[p] = a
                fasta2acc[os.path.basename(p)] = a

    run.truth_written = TRU.write_truth_outputs(
        run.merged, run.mutation_results, run.out_dir,
        read_sources_tsv=run.read_sources_tsv if args.per_read_tsv else None,
        fasta_to_accession=fasta2acc,
        per_read=args.per_read_tsv,
        attempt_to_make_dir=attempt_to_make_dir)
    print()


def report_finish(args, run):
    report_message("-" * 78, "green", initial_indent="  ")
    report_message("Mock metagenome complete!", "green",
                   initial_indent=" " * 28, subsequent_indent="    ",
                   leading_newline=False, trailing_newline=False)
    report_message("-" * 78, "green", initial_indent="  ", leading_newline=False,
                   trailing_newline=True)
 
    num_genomes = len(run.merged)
    total_reads = int(run.merged["assigned_reads"].sum()) if "assigned_reads" in run.merged.columns else 0
 
    if args.type == "paired-end":
        reads_line = f"{run.reads_prefix}_R1.fastq.gz, {run.reads_prefix}_R2.fastq.gz"
        read_unit = "Read-pairs" if total_reads else "Reads"
        reported = total_reads // 2 if total_reads else 0
    else:
        reads_line = f"{run.reads_prefix}.fastq.gz"
        read_unit = "reads"
        reported = total_reads
 
    print(f"      Genomes included:      {num_genomes:,}")
    print(f"      {read_unit} generated:  {reported:,}\n")

    print(f"      Reads:                 {reads_line}")
    gt_root = os.path.join(run.out_dir, "ground-truth")
    print(f"      Ground truth:          {gt_root}{os.sep}")
    print(f"                                 gtdb{os.sep}  (GTDB-taxonomy truth tables)")
    print(f"                                 ncbi{os.sep}  (NCBI-taxonomy truth tables)")
    print("")
