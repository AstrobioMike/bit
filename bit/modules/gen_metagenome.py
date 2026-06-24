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
import pandas as pd # type: ignore
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

    section(f"Phase {n()}: Finding genomes and resolving taxonomy...")
    run = phase_select(args, run)

    section(f"Phase {n()}: Downloading assemblies...")
    run = phase_download(args, run)

    if mutating:
        section(f"Phase {n()}: Mutating genomes...")
    run = phase_mutate(args, run)

    # when not mutating, genome sizes are measured here (needed for abundance)
    if mutating:
        run = phase_abundance(args, run)
    else:
        section(f"Phase {n()}: Getting genome sizes and assigning abundances...")
        run = phase_abundance(args, run)

    section(f"Phase {n()}: Generating reads...")
    run = phase_reads(args, run)
    print()

    section(f"Phase {n()}: Generating final output tables...")
    print()
    phase_truth(args, run)

    write_reproducible_input(args, run)

    report_finish(args, run)


def _phase_counter():
    """ returns a callable yielding 1, 2, 3, ... on each call (for phase labels). """
    state = {"i": 0}
    def nxt():
        state["i"] += 1
        return state["i"]
    return nxt


def check_required_dbs(args, run):
    # all three reference datasets are needed for every run: GTDB (selection +
    # taxonomy), NCBI taxonomy (lineage resolution), and the NCBI assembly
    # summary (download phase resolves accession -> FTP path from it). Fetch them
    # all up front so no large download surprises the user mid-run.
    run.gtdb_dir = get_gtdb_data(quiet=True)
    get_ncbi_tax_data(quiet=True)
    log_data_source(run, "NCBI taxonomy (taxdump) retrieved",
                    _read_retrieved_date(os.environ.get("TAXONKIT_DB")))
    get_ncbi_assembly_data(quiet=True)
    log_data_source(run, "NCBI assembly summary retrieved",
                    _read_retrieved_date(os.environ.get("NCBI_assembly_data_dir")))


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

    # record the resolved seed (whether user-supplied or generated) so the run is
    # reproducible from the runlog even when no --seed was passed.
    seed = getattr(args, "seed", None)
    if seed is not None:
        with open(log_file, "a") as fh:
            fh.write(f"Random seed: {seed}\n")

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
    """ append a reference-data provenance line to the runlog (e.g. GTDB version,
    NCBI assembly-summary retrieval date, NCBI taxdump date). """
    log_file = getattr(run, "log_file", None)
    if not log_file:
        return
    with open(log_file, "a") as fh:
        fh.write(f"{label}: {detail}\n")


def _read_retrieved_date(data_dir):
    """ read a date-retrieved.txt (stored YYYY,MM,DD) from a data dir and return
    it formatted like 'Jun 20, 2026'. Returns '(date unknown)' if absent. """
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
    """ read GTDB version + release date from the version-info file (no printing;
    the value goes to the runlog and the console shows only a spinner). """
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

# columns gen-metagenome actually uses from the GTDB metadata table. The full
# table is ~110 columns; reading only these (via usecols) skips parsing the rest,
# which is both faster and far lighter in memory. Column names vary across GTDB
# releases (e.g. checkm2_* vs checkm_*, and some composition columns absent in
# older releases), so the read intersects this list with the columns actually
# present rather than passing it verbatim (usecols errors on missing names).
GTDB_USED_COLUMNS = [
    "accession", "ncbi_genbank_assembly_accession", "ncbi_taxid",
    "gtdb_representative", "ncbi_refseq_category",
    "domain", "phylum", "class", "order", "family", "genus", "species",
    "checkm2_completeness", "checkm2_contamination",
    "checkm_completeness", "checkm_contamination",
    "genome_size", "contig_count", "gc_count", "gc_percentage", "ambiguous_bases",
    "coding_bases", "coding_density",
]


def _load_gtdb_table(args, run):
    """
    load the GTDB metadata table once (needed for generative selection and/or
    GTDB-taxonomy resolution of user accessions). Cached on run. Silent — the
    caller wraps this in an appropriately-labeled spinner. GTDB version is
    recorded to the runlog here.

    Only the columns gen-metagenome uses are read (see GTDB_USED_COLUMNS); the
    full table has ~110 columns and skipping the unused ones at parse time is
    markedly faster and lighter in memory.
    """
    if getattr(run, "gtdb_tab", None) is not None:
        return run.gtdb_tab
    gtdb_dir = get_gtdb_data(quiet=True)
    gtdb_version, gtdb_release_date = _read_gtdb_version(gtdb_dir)
    log_data_source(run, "GTDB version", f"{gtdb_version} ({gtdb_release_date})")
    gtdb_path = os.path.join(gtdb_dir, "GTDB-arc-and-bac-metadata.tsv")
    # intersect the desired columns with those actually present (column sets vary
    # by GTDB release), so usecols never errors on an absent name.
    present = pd.read_csv(gtdb_path, sep="\t", nrows=0).columns
    usecols = [c for c in GTDB_USED_COLUMNS if c in present]
    run.gtdb_tab = pd.read_csv(
        gtdb_path, sep="\t", low_memory=False, usecols=usecols)
    return run.gtdb_tab


def _assembly_info_path():
    """ path to bit's stored NCBI assembly-info table (taxid source), if set. """
    base = os.environ.get("NCBI_assembly_data_dir")
    if not base:
        return None
    p = os.path.join(base, "ncbi-assembly-info.tsv")
    return p if os.path.exists(p) else None


def _select_with_suppression_screen(gtdb_tab, derep_rank, domains, num_genomes, seed,
                                    user_excluded_groups, assembly_info_path):
    """
    Select num_genomes GTDB genomes, screening each pick against the NCBI
    assembly-summary (absent == suppressed/removed/version-drifted) and backfill-
    ing replacements until num_genomes live genomes are found or the pool is
    exhausted. Excludes dead accessions (not their whole group), so at coarse
    derep ranks a suppressed pick is replaced by a sibling in the same group; at
    species rank (one representative per group) it degrades to a different group.

    Returns (selected_subframe, warnings).
    """
    warnings = []
    on = derep_rank != "off"

    kept_rows = []                 # list of selected GTDB rows (as 1-row frames)
    kept_accs = set()              # accessions already kept (for dedup safety)
    kept_groups = set(user_excluded_groups) if user_excluded_groups else set()
    dead_accs = set()
    seen_accs = set()              # every accession ever picked (kept or dead)

    def acc_of(row_frame):
        return row_frame["ncbi_genbank_assembly_accession"].iloc[0]

    while len(kept_rows) < num_genomes:
        need = num_genomes - len(kept_rows)
        excl_groups = kept_groups if on else None
        sel, sw = SEL._select_gtdb_one_per_rank(
            gtdb_tab, derep_rank=derep_rank, domains=domains,
            num_genomes=need, seed=seed,
            exclude_groups=excl_groups, exclude_accessions=dead_accs or None)
        # keep only the "available" warnings from the FINAL shortfall; intermediate
        # iterations' shortfall warnings are misleading, so collect at the end.
        if len(sel) == 0:
            break

        accs = list(sel["ncbi_genbank_assembly_accession"])

        # stop only when a round surfaces no genome we haven't already seen — i.e.
        # the candidate pool is genuinely exhausted. A round can legitimately add
        # zero *live* genomes (all its head(1) picks suppressed) while still making
        # progress: recording those as dead changes the next round's picks, letting
        # live siblings in the same groups surface
        if all(a in seen_accs for a in accs):
            break
        seen_accs.update(accs)

        # screen this batch
        if assembly_info_path:
            live = TAX.present_accessions(accs, assembly_info_path)
        else:
            live = set(accs)       # no screen possible -> assume all live

        for _, r in sel.iterrows():
            a = r["ncbi_genbank_assembly_accession"]
            if a in live and a not in kept_accs:
                kept_rows.append(r.to_frame().T)
                kept_accs.add(a)
                if on:
                    kept_groups.add(r[derep_rank])
            elif a not in live:
                dead_accs.add(a)

    import pandas as _pd # type: ignore
    selected = (_pd.concat(kept_rows, ignore_index=True)
                if kept_rows else gtdb_tab.iloc[0:0].copy())

    if len(selected) < num_genomes:
        n_dead = len(dead_accs)
        warnings.append(
            f"Requested {num_genomes} genomes; returning {len(selected)} after "
            f"screening out {n_dead} suppressed/unavailable accession(s) and "
            f"exhausting available candidates."
        )

    return selected, warnings


def phase_select(args, run):

    generative = bool(args.num_genomes)

    # GTDB table is needed if selecting generatively OR resolving GTDB taxonomy
    # for user accessions. In generative mode the load is the "pool" of genomes
    # to pick from and gets its own spinner; in accessions-only mode it's loaded
    # silently as part of taxonomy resolution (no pool, no selection step).
    gtdb_tab = None
    print()
    if generative:
        with spinner("Loading pool of potential genomes...",
                     "Loaded pool of potential genomes", indent="    "):
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
        with spinner("Selecting genomes...", "Selected genomes", indent="    "):
            if user_df is not None:
                # resolve GTDB taxonomy for user rows so derep exclusion can use it
                user_df = TAX.fill_gtdb_taxonomy(user_df, gtdb_tab)
            exclude = SEL.user_filled_groups(user_df, args.derep_rank) if user_df is not None else set()

            if args.domain == "both":
                domains = None
            else:
                alias = {"bacteria": "Bacteria", "archaea": "Archaea"}
                domains = [alias[args.domain]]

            sel, sw = _select_with_suppression_screen(
                gtdb_tab, derep_rank=args.derep_rank, domains=domains,
                num_genomes=args.num_genomes, seed=args.seed,
                user_excluded_groups=exclude,
                assembly_info_path=_assembly_info_path())
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
    # taxonomy resolution here, not pooling). version/date provenance -> runlog
    # (assembly-summary date is logged up front in check_required_dbs).
    ai_path = _assembly_info_path()
    with spinner("Resolving taxonomy...", "Resolved taxonomy", indent="    "):
        if gtdb_tab is None:
            gtdb_tab = _load_gtdb_table(args, run)
            if user_df is not None:
                merged = TAX.fill_gtdb_taxonomy(merged, gtdb_tab)
        merged = TAX.resolve_all(merged, gtdb_tab, ai_path)

    run.merged = merged
    print()

    return run


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
        # a single-column file is bare accessions, but tolerate a lone 'accession'
        # header line if present (skip it rather than treat it as an accession).
        if accs and accs[0].lower() == "accession":
            accs = accs[1:]
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

    # resolve downloaded fasta paths, and reconcile the accession in run.merged to
    # the version actually downloaded. NCBI may serve a newer version (or the GCF
    # twin of a requested GCA) than the accession we asked for, and the file is
    # named after that actual accession; carrying the as-downloaded accession
    # forward keeps the genome files, path maps, and truth tables all consistent.
    new_accessions = []
    new_paths = {}
    for a in run.merged["accession"]:
        path = find_downloaded_fasta(run.genomes_dir, a)
        actual = _accession_from_fasta_name(path) if path else None
        eff = actual or a                    # fall back to requested if unparsable
        new_accessions.append(eff)
        new_paths[eff] = path
    run.merged["accession"] = new_accessions
    run.genome_paths = new_paths

    # drop genomes that didn't download (e.g. suppressed/removed accessions) so
    # later phases don't try to size/mutate/generate reads for absent files. The
    # download step already reported which accessions weren't found.
    missing = [a for a in run.merged["accession"] if run.genome_paths.get(a) is None]
    if missing:
        run.merged = run.merged[~run.merged["accession"].isin(missing)].reset_index(drop=True)
        for a in missing:
            run.genome_paths.pop(a, None)

    if len(run.merged) == 0:
        report_message("No genomes downloaded successfully; nothing to do.", "red",
                       initial_indent="    ", subsequent_indent="    ")
        raise SystemExit(1)

    run.working_paths = dict(run.genome_paths)

    # rename dl-ncbi-assemblies' info table for consistency with the gtdb summary
    _rename_if_exists(
        os.path.join(run.genomes_dir, "wanted-ncbi-accessions-info.tsv"),
        os.path.join(run.genomes_dir, "selected-genomes-ncbi-info.tsv"))

    # write the GTDB info summary for the final community
    write_gtdb_summary(run)

    # the full GTDB metadata table (~GBs) is no longer needed after this point —
    # selection, taxonomy resolution, and the GTDB summary are all done. Drop it
    # so it isn't held through mutate/abundance/reads/truth (the per-read truth
    # build is memory-heavy and shouldn't stack on top of the GTDB table).
    run.gtdb_tab = None

    print()
    return run


def _rename_if_exists(src, dst):
    if os.path.exists(src):
        os.replace(src, dst)


# columns for the per-genome GTDB info summary written to genomes/
GTDB_SUMMARY_EXTRA = ["genome_size", "contig_count", "gc_count", "gc_percentage", "ambiguous_bases", "checkm2_completeness", "checkm2_contamination",
                      "coding_bases", "coding_density"]
GTDB_SUMMARY_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def write_gtdb_summary(run):
    """
    Write genomes/selected-genomes-gtdb-info.tsv: GTDB taxonomy + quality/
    composition columns for each genome in the final community, keyed by
    run.merged['accession']. Genomes not in the GTDB table (e.g. user-supplied
    accessions absent from GTDB) get NA for all GTDB columns. checkm2 columns are
    emitted strictly (NA if the source lacks checkm2_*, never back-filled from
    an older checkm_* column).
    """
    gtdb_tab = getattr(run, "gtdb_tab", None)
    out_cols = ["accession"] + GTDB_SUMMARY_EXTRA + GTDB_SUMMARY_RANKS
    accs = list(run.merged["accession"])

    # build a lookup from the GTDB table keyed by numeric accession core so user
    # GCF/GCA accessions match GTDB's GCA-form rows.
    lookup = {}
    if gtdb_tab is not None and "ncbi_genbank_assembly_accession" in gtdb_tab.columns:
        cols_present = [c for c in (GTDB_SUMMARY_EXTRA + GTDB_SUMMARY_RANKS)
                        if c in gtdb_tab.columns]
        gcol = gtdb_tab["ncbi_genbank_assembly_accession"].map(_acc_digits_core)
        for key, row in zip(gcol.values, gtdb_tab[cols_present].to_dict("records")):
            lookup.setdefault(key, row)

    rows = []
    for a in accs:
        rec = {"accession": a}
        src = lookup.get(_acc_digits_core(a), {})
        for c in GTDB_SUMMARY_EXTRA + GTDB_SUMMARY_RANKS:
            rec[c] = src.get(c, pd.NA)
        rows.append(rec)

    out_path = os.path.join(run.genomes_dir, "selected-genomes-gtdb-info.tsv")
    pd.DataFrame(rows, columns=out_cols).to_csv(out_path, sep="\t", index=False)
    run.gtdb_summary_path = out_path


def _acc_digits_core(acc):
    """
    numeric core of an accession (drop RS_/GB_ prefix, GCA_/GCF_ prefix, and
    version), so GCA/GCF twins and version variants collapse to one key
    """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    s = s.split(".")[0]
    if "_" in s:
        s = s.split("_", 1)[1]
    return s


def _accession_from_fasta_name(path):
    """
    recover the assembly accession from a dl-ncbi-assemblies filename, which is
    '<accession>.fasta.gz' (e.g. 'GCF_000005845.2.fasta.gz' -> 'GCF_000005845.2').
    Returns None if path is falsy
    """
    if not path:
        return None
    fn = os.path.basename(path)
    for ext in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".gz"):
        if fn.endswith(ext):
            return fn[: -len(ext)]
    return fn


def find_downloaded_fasta(genomes_dir, accession):
    """
    locate the downloaded fasta for an accession (dl-ncbi-assemblies naming).

    Match on the versionless numeric accession core rather than a literal prefix:
    NCBI names the file after the assembly's CURRENT version and accession form,
    which can differ from the accession string carried in run.merged (GTDB may
    reference an older version, and a GCA pick can download as its GCF twin or
    vice versa). A strict startswith(accession) check misses those and would
    wrongly treat a present genome as not-downloaded.

    dl-ncbi-assemblies writes files as '<dl_acc>.fasta.gz' (e.g.
    'GCF_000005845.2.fasta.gz'), so the accession core is the filename with its
    extension stripped
    """
    want = _acc_digits_core(accession)
    for fn in os.listdir(genomes_dir):
        if not (fn.endswith(".fasta") or fn.endswith(".gz")):
            continue
        # strip .fasta.gz / .fasta / .gz to recover the leading accession token
        file_acc = fn
        for ext in (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".gz"):
            if file_acc.endswith(ext):
                file_acc = file_acc[: -len(ext)]
                break
        if _acc_digits_core(file_acc) == want:
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
        print()
        mut_dir = os.path.join(run.out_dir, "mutated-genomes")
        with tqdm(total=len(accs), desc="    Mutating", ncols=78, unit=" genome") as pbar:
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
    with ThreadPoolExecutor(max_workers=max(int(jobs), 1)) as pool:
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
        include_Ns=args.include_Ns, genome_sizes=genome_sizes, jobs=args.jobs,
        per_read_tsv=args.per_read_tsv)
    generate_reads(gr_args)
    run.read_sources_tsv = f"{run.reads_prefix}-read-sources.tsv.gz"
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

    gt_root = os.path.join(run.out_dir, "ground-truth")
    attempt_to_make_dir(gt_root)

    run.truth_written = {}

    # chunk size for the per-read build: drives both the actual chunked read
    # (memory bounded by this many rows) and the progress-bar granularity. Smaller
    # chunks update the bar more often at negligible speed cost and lower peak
    # memory; this value is the single source of truth for both.
    chunksize = 500_000
    # total chunks per taxonomy, for a determinate progress bar over the per-read
    # build (rows == total reads, available from the assigned_reads column).
    total_reads = int(run.merged["assigned_reads"].sum()) if "assigned_reads" in run.merged.columns else 0
    n_chunks = max(1, -(-total_reads // chunksize)) if total_reads else None  # ceil-div

    jobs = max(1, int(getattr(args, "jobs", 1) or 1))

    if not args.per_read_tsv:
        # no per-read table -> just the quick per-genome/per-rank tables, serial.
        for taxonomy in TRU.TAXONOMIES:
            label = f"Building {taxonomy.upper()} truth tables..."
            with spinner(label, label, indent="    "):
                run.truth_written[taxonomy] = TRU.build_truth_for_taxonomy(
                    taxonomy, run.merged, run.mutation_results, gt_root,
                    read_sources_tsv=None, fasta_to_accession=fasta2acc,
                    per_read=False, attempt_to_make_dir=attempt_to_make_dir)
    elif jobs <= 1 or len(TRU.TAXONOMIES) < 2:
        # serial per-read build with a determinate per-taxonomy bar.
        for taxonomy in TRU.TAXONOMIES:
            print(f"    Building {taxonomy.upper()} per-read table...")
            pbar = tqdm(total=n_chunks, desc="      Progress", unit="chunk", ncols=76)
            run.truth_written[taxonomy] = TRU.build_truth_for_taxonomy(
                taxonomy, run.merged, run.mutation_results, gt_root,
                read_sources_tsv=run.read_sources_tsv, fasta_to_accession=fasta2acc,
                per_read=True, attempt_to_make_dir=attempt_to_make_dir,
                progress=pbar.update, chunksize=chunksize)
            if pbar.total is not None:
                pbar.update(pbar.total - pbar.n)
            pbar.close()
            print()
    else:
        # parallel per-read build: the two taxonomies are independent (different
        # output, shared read-only source), built in separate spawned processes
        # for true multi-core parallelism on the CPU-bound merge+gzip. A shared
        # counter drives one aggregate progress bar over both taxonomies' chunks.
        total = (n_chunks * len(TRU.TAXONOMIES)) if n_chunks else None
        print("    Building per-read tables (gtdb, ncbi)...")
        pbar = tqdm(total=total, desc="    Progress", unit=" chunk", ncols=78)
        written = _build_truth_processes(run, gt_root, fasta2acc, chunksize, pbar)
        run.truth_written.update(written)
        if pbar.total is not None:
            pbar.update(pbar.total - pbar.n)
        pbar.close()
        print()

    # remove intermediates now that the truth tables hold their information:
    #   - read-sources TSV: consumed by build_read_truth (per-read tables)
    #   - per-genome-coverage TSV: consumed by gen-reads; its columns live in
    #     each truth-per-genome.tsv (assigned_coverage/reads/rel_abundance/sizes)
    if args.per_read_tsv and os.path.exists(run.read_sources_tsv):
        os.remove(run.read_sources_tsv)
    if os.path.exists(run.coverage_tsv):
        os.remove(run.coverage_tsv)


def _truth_process_worker(spec, counter, result_q):
    """ picklable worker for the process-parallel path: build one taxonomy's
    per-read table, incrementing the shared counter once per chunk so the parent
    can render a single aggregate progress bar, and putting the result (or an
    error marker) on the queue. counter and result_q are passed as direct Process
    args (not nested in spec) so spawn can share them. """
    try:
        import sys as _sys
        _sys.path.insert(0, os.getcwd())
        from bit.modules.gen_mg import truth as _TRU

        def tick():
            with counter.get_lock():
                counter.value += 1

        entry = _TRU.build_truth_for_taxonomy(
            spec["taxonomy"], spec["merged"], spec["mutation_results"], spec["gt_root"],
            read_sources_tsv=spec["read_sources_tsv"],
            fasta_to_accession=spec["fasta_to_accession"],
            per_read=True, attempt_to_make_dir=os.makedirs,
            progress=tick, chunksize=spec["chunksize"])
        result_q.put((spec["taxonomy"], entry))
    except Exception as e:                           # surface worker failures
        result_q.put(("__error__", f"{spec['taxonomy']}: {type(e).__name__}: {e}"))


def _build_truth_processes(run, gt_root, fasta2acc, chunksize, pbar):
    """ build the two taxonomies' per-read tables in separate (spawned) processes
    for true multi-core parallelism on the CPU-bound merge+gzip. A shared counter
    (incremented per chunk by each worker) drives the single aggregate bar; results
    come back over a queue. Raw Process (not a pool) so the shared Value/Queue can
    be passed as direct args at construction rather than pickled through submit. """
    import time
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    counter = ctx.Value("i", 0)
    result_q = ctx.Queue()

    specs = [{
        "taxonomy": taxonomy, "merged": run.merged,
        "mutation_results": run.mutation_results, "gt_root": gt_root,
        "read_sources_tsv": run.read_sources_tsv, "fasta_to_accession": fasta2acc,
        "chunksize": chunksize,
    } for taxonomy in TRU.TAXONOMIES]

    procs = [ctx.Process(target=_truth_process_worker, args=(s, counter, result_q))
             for s in specs]
    for p in procs:
        p.start()

    # advance the bar from the shared counter while workers run
    while any(p.is_alive() for p in procs):
        with counter.get_lock():
            done = counter.value
        if pbar.n < done:
            pbar.update(done - pbar.n)
        time.sleep(0.2)
    for p in procs:
        p.join()

    written = {}
    while not result_q.empty():
        taxonomy, entry = result_q.get()
        if taxonomy == "__error__":
            raise RuntimeError(f"per-read truth build failed ({entry})")
        written[taxonomy] = entry
    if len(written) != len(TRU.TAXONOMIES):
        missing = set(TRU.TAXONOMIES) - set(written)
        raise RuntimeError(f"per-read truth build failed for: {', '.join(missing)}")
    return written


def write_reproducible_input(args, run):
    """
    Write reproducible-input.tsv in the run's root: every genome in the final
    community (user-supplied and randomly-selected) with both the rel_abundance
    and coverage columns, plus a mutation_rate column when mutation was on. Fed
    back as a pure --accessions file (no --num-genomes), it recreates this
    community's composition; the abundance mode chosen at re-run time selects
    which abundance column is honored. Identical reads additionally require
    the same --seed (recorded in runlog.txt).
    """
    merged = run.merged
    if merged is None or len(merged) == 0:
        return

    cols = {"accession": list(merged["accession"])}
    if "assigned_rel_abundance" in merged.columns:
        cols["rel_abundance"] = list(merged["assigned_rel_abundance"])
    if "assigned_coverage" in merged.columns:
        cols["coverage"] = list(merged["assigned_coverage"])

    if args.mutation_mode != "off" and run.mutation_results:
        cols["mutation_rate"] = [
            run.mutation_results.get(a, {}).get("rate", pd.NA)
            for a in merged["accession"]
        ]

    out_path = os.path.join(run.out_dir, "reproducible-input.tsv")
    pd.DataFrame(cols).to_csv(out_path, sep="\t", index=False)
    run.reproducible_input_path = out_path


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
    print(f"      Ground-truth files:    {gt_root}{os.sep}")
    print(f"                                 gtdb{os.sep}  (GTDB-taxonomy truth tables)")
    print(f"                                 ncbi{os.sep}  (NCBI-taxonomy truth tables)")
    print("")
