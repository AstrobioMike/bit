import os
import shutil
import tempfile
import pandas as pd # type: ignore
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from Bio import SeqIO # type: ignore
from tqdm import tqdm # type: ignore
from dataclasses import dataclass, field
from pathlib import Path
import json
from subprocess import run
import pysam # type: ignore
from collections import defaultdict
import numpy as np # type: ignore
from natsort import natsort_keygen # type: ignore
from bit.modules.input_parsing import get_input_reads_dict_from_dir
from bit.modules.seqs import identify_seq_type
from bit.modules.general import (report_message,
                       report_failure,
                       get_package_path,
                       color_text,
                       log_command_run,
                       transient_spinner)


def run_assembly(args, full_cmd_executed):
    outputs_dir = args.output_prefix + "-outputs"
    assembly_path_dict, targets_plan = assembly_preflight(args, outputs_dir)
    run_assembly_screen(args, assembly_path_dict, outputs_dir, targets_plan)
    report_assembly_screen_finished(args, outputs_dir)
    log_command_run(full_cmd_executed, outputs_dir)


def run_reads(args, full_cmd_executed):
    reads_dict = get_input_reads_dict_from_dir(args.reads_dir)
    reads_config = ReadsRunConfiguration.from_args(args)
    run_reads_snakemake(reads_config, reads_dict)
    report_read_screen_finished(args)
    log_command_run(full_cmd_executed, reads_config.log_files_dir)


def assembly_preflight(args, outputs_dir):

    check_assembly_inputs(args.assemblies)

    if os.path.exists(outputs_dir):
        shutil.rmtree(outputs_dir)
    os.makedirs(outputs_dir)

    targets_plan = build_targets_plan(args.targets, outputs_dir)

    assembly_basenames = [os.path.splitext(os.path.basename(assembly))[0] for assembly in args.assemblies]
    basename_counts = {basename: assembly_basenames.count(basename) for basename in assembly_basenames}

    if any(count > 1 for count in basename_counts.values()):
        assembly_path_dict = {assembly: assembly for assembly in args.assemblies}
    else:
        assembly_path_dict = {assembly: os.path.splitext(os.path.basename(assembly))[0] for assembly in args.assemblies}

    return assembly_path_dict, targets_plan


def check_assembly_inputs(assemblies):

    for assembly in assemblies:
        if not os.path.exists(assembly):
            report_failure(f"Specified input assembly file not found: {assembly}")
        if os.path.isdir(assembly):
            report_failure(f"Specified input assembly is a directory, but needs to be a file or files: {assembly}")


def parses_as_fasta(path):
    """
    checks if at least a fasta record can be read from path. This is
    used to classify --targets entries
    """
    if not os.path.isfile(path):
        return False
    try:
        for _ in SeqIO.parse(path, "fasta"):
            return True
    except Exception:
        return False
    return False


def blast_db_is_readable(db_prefix, dbtype):
    """
    returns True if blastdbcmd can open the db at db_prefix as dbtype
    ('nucl' or 'prot')
    """
    try:
        subprocess.run(
            ["blastdbcmd", "-db", db_prefix, "-dbtype", dbtype, "-info"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def diamond_db_is_readable(path):
    """
    returns True if `diamond dbinfo` can open path as a DIAMOND database.
    handles the path given with or without the .dmnd extension
    """
    candidates = [path]
    if path.endswith(".dmnd"):
        candidates.append(path[:-5])
    else:
        candidates.append(path + ".dmnd")
    for candidate in candidates:
        try:
            subprocess.run(
                ["diamond", "dbinfo", "--db", candidate],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    return False


def classify_targets_entry(entry):
    """
    classify one --targets entry into (modality, seqtype):
        ("fasta", "nt")  nucleotide fasta   -> build blast db -> blastn -db
        ("fasta", "aa")  protein fasta      -> build .dmnd -> diamond blastx
        ("db",    "nt")  nucleotide BLAST db-> blastn -db
        ("db",    "aa")  DIAMOND db (.dmnd) -> diamond blastx -db

        protein BLAST dbs and unrecognizable inputs are rejected
    """

    if parses_as_fasta(entry):
        return ("fasta", identify_seq_type(entry))

    if diamond_db_is_readable(entry):
        return ("db", "aa")

    if blast_db_is_readable(entry, "nucl"):
        return ("db", "nt")

    if blast_db_is_readable(entry, "prot"):
        report_message(
            f"The follwing `--targets` entry looks like a protein BLAST database, which isn't "
            f"supported. Pass a protein fasta or a DIAMOND database (built with `diamond makedb`) instead.")
        report_message(f"    {entry}", join = False, color = "none")
        report_failure()

    report_message(
        f"The following `--targets` entry is neither a readable fasta nor a recognized BLAST or "
        f"DIAMOND database:")
    report_message(f"    {entry}", join = False, color = "none")
    report_failure()


def build_targets_plan(targets_entries, outputs_dir):
    """
    classifies all --targets entries, pull each one's {target: length}, and
    merge into a single plan used for the whole run

    applies __nt/__aa suffixes to target names that collide across
    sequence types. within-type duplicate names are left
    as-is

    returns a dict:
        {
        "entries": [ {"entry","modality","seqtype","names": set(...)}, ... ],
        "lengths": {target: length},          # merged, suffixed
        "types":   {target: "nt"|"aa"},        # merged, suffixed
        }
    where each entry's 'names' are the SUFFIXED names for that entry, so the
    search step can remap its raw sseqids to them
    """

    classified = []
    for entry in targets_entries:
        modality, seqtype = classify_targets_entry(entry)
        raw_lengths = get_targets(entry, modality, seqtype)
        classified.append({
            "entry": entry,
            "modality": modality,
            "seqtype": seqtype,
            "raw_lengths": raw_lengths,
        })

    # find cross-type name collisions
    names_by_type = {"nt": set(), "aa": set()}
    for c in classified:
        names_by_type[c["seqtype"]].update(c["raw_lengths"].keys())
    colliding = names_by_type["nt"] & names_by_type["aa"]

    merged_lengths = {}
    merged_types = {}
    entries_out = []

    for c in classified:
        st = c["seqtype"]
        suffix = f"__{st}"
        suffixed_names = set()
        for name, length in c["raw_lengths"].items():
            key = name + suffix if name in colliding else name
            merged_lengths[key] = length
            merged_types[key] = st
            suffixed_names.add(key)
        entries_out.append({
            "entry": c["entry"],
            "modality": c["modality"],
            "seqtype": st,
            "names": suffixed_names,
            "colliding_raw": {n for n in c["raw_lengths"] if n in colliding},
        })

    return {"entries": entries_out, "lengths": merged_lengths, "types": merged_types}


def run_assembly_screen(args, assembly_path_dict, outputs_dir, targets_plan):

    targets_dict = targets_plan["lengths"]
    type_map = targets_plan["types"]
    summary_df = pd.DataFrame()
    total_assemblies = len(args.assemblies)
    task = "megablast" if args.use_megablast else "blastn"

    print("")

    for i, assembly in enumerate(tqdm(args.assemblies, total = total_assemblies,
                                      desc = "Processing assemblies", unit = "assembly", ncols=76,
                                      bar_format="    {l_bar}{bar}| {n_fmt}/{total_fmt} processed [{elapsed}]")):

        unique_assembly_name = assembly_path_dict[assembly]
        out_base = safe_name(unique_assembly_name)

        per_entry_filtered = []
        per_entry_loci = []

        for entry_info in targets_plan["entries"]:
            entry = entry_info["entry"]
            seqtype = entry_info["seqtype"]
            min_perc_id = args.min_aa_perc_id if seqtype == "aa" else args.min_nt_perc_id

            blast_df = run_targets_search(assembly, entry_info, outputs_dir,
                                          out_base, args.num_jobs, args.num_threads, task)

            with transient_spinner("resolving hits", indent="        "):
                blast_df = apply_collision_suffix(blast_df, entry_info)

                # guard against id mismatch (e.g. a db built without -parse_seqids):
                # check against THIS entry's target names only
                check_sseqid_target_match(blast_df, entry_info["names"])

                filtered = filter_blast_results(blast_df, min_perc_id)
                per_entry_filtered.append(filtered)

                loci = resolve_assembly_loci(filtered, targets_dict,
                                             gap_tolerance=args.hit_merge_gap,
                                             min_perc_cov=args.min_perc_cov,
                                             min_edge_perc_cov=args.min_edge_perc_cov,
                                             edge_tolerance=args.edge_tolerance)
                if not loci.empty:
                    loci = loci.assign(target_type=loci["target"].map(type_map))
                else:
                    loci = loci.assign(target_type=pd.Series(dtype=object))
                per_entry_loci.append(loci)

        # merge across entries
        with transient_spinner("finalizing results", indent="        "):
            nonempty_filtered = [d for d in per_entry_filtered if not d.empty]
            nonempty_loci = [d for d in per_entry_loci if not d.empty]
            filtered_hits_df = (pd.concat(nonempty_filtered, ignore_index=True)
                                if nonempty_filtered else pd.DataFrame())
            loci_df = (pd.concat(nonempty_loci, ignore_index=True)
                       if nonempty_loci else pd.DataFrame())

            # per-assembly filtered BLAST table (pident-passing HSPs, all entries)
            filtered_hits_df.to_csv(f"{outputs_dir}/{out_base}-filtered-search-results.tsv",
                                    sep="\t", index=False)

            # contig lengths (not carried on loci_df) for the contig summary
            if filtered_hits_df.empty:
                contig_lengths = {}
            else:
                contig_lengths = (filtered_hits_df.drop_duplicates("qseqid")
                                  .set_index("qseqid")["qlen"].to_dict())

            # optional across-target region resolution: collapse loci from different
            # targets (including across nt/aa types) that pile onto the same contig
            # locus into a single called region, keeping the best and folding the
            # rest into 'other_targets'. on by default; disabled with --dont-resolve-regions.
            contig_island_map = {}
            if args.resolve_regions:
                region_df, num_regions_by_contig = gen_region_calls_table(
                    loci_df, overlap_frac=args.region_overlap_frac,
                    contig_lengths=contig_lengths)

                # optional island extraction: cluster regions into islands and cut
                # out those buried in large contigs (needs resolved regions). adds
                # cross-reference ids tying region-calls <-> contig-summary <-> the
                # extracted FASTAs + manifest. off via --no-island-extraction.
                if args.island_extraction:
                    region_island_map, contig_island_map = extract_islands(
                        region_df, contig_lengths, assembly, outputs_dir, out_base,
                        island_gap=args.island_gap, contig_floor=args.island_contig_floor,
                        ratio=args.island_contig_ratio, min_span=args.island_min_span,
                        min_regions=args.island_min_regions, buffer=args.island_buffer)
                    region_df["island_id"] = region_df.index.map(
                        lambda idx: region_island_map.get(idx, "none"))

                region_df.to_csv(f"{outputs_dir}/{out_base}-region-calls.tsv",
                                 sep="\t", index=False)
            else:
                num_regions_by_contig = None

            # per-assembly contig summary table
            contig_df = gen_contig_summary_table(loci_df, contig_lengths,
                                                 num_regions_by_contig,
                                                 contig_island_map=contig_island_map)
            contig_df.to_csv(f"{outputs_dir}/{out_base}-hit-contig-summary.tsv",
                             sep="\t", index=False)

            summary_df = update_assembly_summary_table(loci_df,
                                              targets_dict,
                                              unique_assembly_name,
                                              summary_df)

    if not args.report_all_targets:
        summary_df = filter_undetected_assembly_targets(summary_df)

    output_tsv = args.output_prefix + "-summary.tsv"
    if args.assemblies_as_rows:
        summary_df.to_csv(output_tsv, sep = "\t", index_label = "input-assembly")
    else:
        summary_df = summary_df.T
        summary_df.to_csv(output_tsv, sep = "\t", index_label = "target")


def apply_collision_suffix(blast_df, entry_info):
    """
    remap an entry's raw sseqids onto the merged ones. any sseqid whose
    raw name collided across types gets this entry's __nt/__aa suffix, so it
    matches the suffixed key in targets_dict. non-colliding ids are
    untouched
    """

    colliding_raw = entry_info.get("colliding_raw", set())
    if blast_df.empty or not colliding_raw:
        return blast_df
    suffix = f"__{entry_info['seqtype']}"
    blast_df = blast_df.copy()
    blast_df["sseqid"] = blast_df["sseqid"].map(
        lambda s: s + suffix if s in colliding_raw else s)
    return blast_df


def get_targets(targets, modality, seqtype):
    """
    returns {target_id: length} for one --targets entry

    fasta -> parse the file (length in residues: nt for nucleotide, aa for protein)
    nt db -> blastdbcmd (nucleotide)
    aa db -> diamond getseq / dbinfo to recover accessions + lengths
    """

    if modality == "fasta":
        return {record.id: len(record.seq) for record in SeqIO.parse(targets, "fasta")}

    if seqtype == "nt":
        # nucleotide BLAST db: one "<accession> <length>" line per sequence
        result = subprocess.run(
            ["blastdbcmd", "-db", targets, "-dbtype", "nucl",
             "-entry", "all", "-outfmt", "%a %l"],
            check=True, capture_output=True, text=True)
        targets_dict = {}
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            acc, length = line.rsplit(" ", 1)
            targets_dict[acc] = int(length)
        return targets_dict

    # aa DIAMOND db: dump sequences and measure lengths. diamond getseq writes
    # the db's sequences as fasta to stdout; we parse ids + lengths from it.
    result = subprocess.run(
        ["diamond", "getseq", "--db", targets, "--out", "/dev/stdout"],
        check=True, capture_output=True, text=True)
    targets_dict = {}
    for record in SeqIO.parse(StringIO(result.stdout), "fasta"):
        targets_dict[record.id] = len(record.seq)
    return targets_dict


def run_targets_search(assembly, entry_info, outputs_dir, out_base, num_jobs, num_threads, task):
    """
    dispatches one targets entry to the right search engine and return a
    BLAST-outfmt-6-shaped DataFrame. nucleotide -> blastn; protein -> DIAMOND
    blastx (building a blast db or .dmnd in outputs_dir first if the entry is
    a fasta, so searches always run in db mode)

    nucleotide searches are parallelized as num_jobs blastn processes with
    num_threads threads each; DIAMOND gets num_jobs * num_threads threads in one process
    """

    entry = entry_info["entry"]
    modality = entry_info["modality"]
    seqtype = entry_info["seqtype"]
    # tag outputs per-entry so multiple targets files don't clobber each other
    entry_tag = safe_name(os.path.splitext(os.path.basename(entry))[0])

    if seqtype == "nt":
        if modality == "fasta":
            db_path = f"{outputs_dir}/{entry_tag}-blast-db"
            if not os.path.exists(f"{db_path}.ndb"):
                build_blast_db(entry, db_path)
        else:
            db_path = entry  # already a blast db
        return run_blastn(assembly, db_path, outputs_dir,
                          f"{out_base}-{entry_tag}", num_jobs, num_threads, task)

    # protein
    if modality == "fasta":
        db_path = f"{outputs_dir}/{entry_tag}.dmnd"
        if not os.path.exists(db_path):
            build_diamond_db(entry, db_path)
    else:
        db_path = entry  # already a .dmnd db
    return run_diamond_blastx(assembly, db_path, outputs_dir,
                              f"{out_base}-{entry_tag}", num_jobs * num_threads)


BLASTN_OUTFMT_FIELDS = ["qseqid", "qlen", "sseqid", "slen", "qstart", "qend",
                        "sstart", "send", "length", "qcovs", "qcovhsp",
                        "qcovus", "pident", "evalue", "bitscore"]


def build_blast_db(nt_fasta, db_path):

    cmd = ["makeblastdb", "-in", nt_fasta, "-dbtype", "nucl", "-out", db_path]
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        print()
        report_message("makeblastdb failed with the following error:")
        report_message(f"{e.output}", join = False, color = "none", initial_indent = "      ", subsequent_indent = "      ")
        report_failure("")


def find_next_record_start(handle, offset, filesize):
    """
    returns the byte offset of the first fasta record header at or after
    offset (a line starting with '>'), or filesize if none remains
    """
    if offset == 0:
        return 0
    handle.seek(offset - 1)
    prev = handle.read(1)
    while True:
        chunk = handle.read(1024 * 1024)
        if not chunk:
            return filesize
        search_space = prev + chunk
        idx = search_space.find(b"\n>")
        if idx != -1:
            return handle.tell() - len(chunk) + (idx - len(prev)) + 1
        prev = search_space[-1:]


def split_fasta_by_bytes(fasta_path, num_chunks, dest_dir):
    """
    splits a plain-text fasta into up to num_chunks piece-files on record
    boundaries, targeting roughly equal byte sizes
    """
    filesize = os.path.getsize(fasta_path)
    if num_chunks <= 1 or filesize == 0:
        return [fasta_path]

    with open(fasta_path, "rb") as f:
        boundaries = [find_next_record_start(f, (filesize * i) // num_chunks, filesize)
                      for i in range(num_chunks)]
    boundaries.append(filesize)

    spans = [(start, end) for start, end in zip(boundaries, boundaries[1:])
             if end > start]

    chunk_paths = []
    with open(fasta_path, "rb") as f:
        for i, (start, end) in enumerate(spans):
            chunk_path = os.path.join(dest_dir, f"query-chunk-{i:04d}.fasta")
            f.seek(start)
            remaining = end - start
            with open(chunk_path, "wb") as out:
                while remaining > 0:
                    buf = f.read(min(4 * 1024 * 1024, remaining))
                    if not buf:
                        break
                    out.write(buf)
                    remaining -= len(buf)
            chunk_paths.append(chunk_path)

    return chunk_paths


def run_single_blastn(query_path, db_path, out_path, num_threads, task,
                      remove_query=False):
    """
    runs one blastn process writing outfmt-6 to out_path

    'remove_query' is about the temp files when chunking is happening, not
    the initial user input

    returns (returncode, captured_output)
    """
    outfmt = "6 " + " ".join(BLASTN_OUTFMT_FIELDS)
    cmd = ["blastn", "-task", task, "-query", query_path, "-db", db_path,
           "-num_threads", str(num_threads), "-outfmt", outfmt, "-out", out_path]
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    if remove_query:
        try:
            os.remove(query_path)
        except OSError:
            pass
    return result.returncode, result.stdout


def report_blast_failure(captured_output):
    print()
    report_message("BLAST failed with the following error:")
    report_message(f"{captured_output}", join = False, color = "none",
                   initial_indent = "      ", subsequent_indent = "      ")
    report_failure("")


TARGET_CHUNK_BYTES = 50 * 1024 * 1024  # ~50MB query per blastn job (floor on chunk size)
MAX_CHUNKS = 1000
CHUNK_MULTIPLIER = 2 # chunks are capped at CHUNK_MULTIPLIER * num_jobs


def run_blastn(assembly, db_path, outputs_dir, out_base, num_jobs, num_threads, task):
    """
    blastn of assembly against db_path

    with num_jobs > 1 the query is split into ~TARGET_CHUNK_BYTES chunks (with
    a floor of num_jobs and a cap of MAX_CHUNKS) and run as num_jobs parallel
    blastn processes with num_threads threads each, so large inputs yield many
    short-lived jobs and the progress bar (tracked in bytes of query processed)
    advances steadily instead of sitting on a few huge chunks. chunk query and
    output files are deleted as they're consumed to bound scratch disk.

    results stream to disk via -out rather than being captured in memory
    """

    cols = BLASTN_OUTFMT_FIELDS
    filesize = os.path.getsize(assembly)

    with tempfile.TemporaryDirectory(dir=outputs_dir, prefix="blast-tmp-") as tmp_dir:

        raw_results_path = os.path.join(tmp_dir, "all-results.tsv")

        if num_jobs <= 1:
            with transient_spinner("running blastn", indent="        "):
                returncode, output = run_single_blastn(assembly, db_path,
                                                       raw_results_path,
                                                       num_threads, task)
            if returncode != 0:
                report_blast_failure(output)
        else:
            size_based = -(-filesize // TARGET_CHUNK_BYTES)
            num_chunks = min(num_jobs * CHUNK_MULTIPLIER, size_based, MAX_CHUNKS)
            num_chunks = max(num_chunks, num_jobs)  # floor so workers aren't idle
            chunk_paths = split_fasta_by_bytes(assembly, num_chunks, tmp_dir)
            chunk_outs = [os.path.join(tmp_dir, f"chunk-{i:04d}-results.tsv")
                          for i in range(len(chunk_paths))]

            if len(chunk_paths) == 1:
                # single record / unsplittable: one process, full thread budget
                with transient_spinner("running blastn", indent="        "):
                    returncode, output = run_single_blastn(chunk_paths[0], db_path,
                                                           raw_results_path,
                                                           num_jobs * num_threads,
                                                           task, remove_query=True)
                if returncode != 0:
                    report_blast_failure(output)
            else:
                chunk_sizes = {i: os.path.getsize(chunk_path)
                               for i, chunk_path in enumerate(chunk_paths)}
                with ThreadPoolExecutor(max_workers=num_jobs) as pool:
                    futures = {pool.submit(run_single_blastn, chunk_path, db_path,
                                           chunk_out, num_threads, task,
                                           remove_query=True): i
                               for i, (chunk_path, chunk_out)
                               in enumerate(zip(chunk_paths, chunk_outs))}
                    progress = tqdm(total=filesize, desc="        blastn progress",
                                    unit="B", unit_scale=True, unit_divisor=1024,
                                    ncols=76, leave=False)
                    for future in as_completed(futures):
                        returncode, output = future.result()
                        if returncode != 0:
                            report_blast_failure(output)
                        progress.update(chunk_sizes[futures[future]])
                    progress.close()

                with transient_spinner("combining results", indent="        "):
                    with open(raw_results_path, "wb") as final:
                        for chunk_out in chunk_outs:
                            with open(chunk_out, "rb") as f:
                                shutil.copyfileobj(f, final, 4 * 1024 * 1024)
                            try:
                                os.remove(chunk_out)
                            except OSError:
                                pass

        if os.path.getsize(raw_results_path) == 0:
            blast_df = pd.DataFrame(columns=cols)
        else:
            with transient_spinner("loading results", indent="        "):
                blast_df = pd.read_csv(raw_results_path, sep="\t",
                                       header=None, names=cols)

    with transient_spinner("writing search results", indent="        "):
        blast_df["perc-subj-cov"] = (round((blast_df["length"] / blast_df["slen"]) * 100, 1)
                                     if not blast_df.empty else pd.Series(dtype=float))

        blast_df.to_csv(f"{outputs_dir}/{out_base}-blast-results.tsv", sep="\t", index=False)

    return blast_df


def build_diamond_db(protein_fasta, db_path):

    cmd = ["diamond", "makedb", "--in", protein_fasta, "--db", db_path, "--quiet"]
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        print()
        report_message("DIAMOND makedb failed with the following error:")
        report_message(f"{e.output}", join = False, color = "none", initial_indent = "      ", subsequent_indent = "      ")
        report_failure("")


def run_diamond_blastx(assembly, db_path, outputs_dir, out_base, num_threads):

    fields = ["qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart",
              "send", "length", "pident", "evalue", "bitscore"]

    with tempfile.TemporaryDirectory(dir=outputs_dir, prefix="diamond-tmp-") as tmp_dir:

        raw_results_path = os.path.join(tmp_dir, "all-results.tsv")

        cmd = [
            "diamond", "blastx",
            "--query", assembly,
            "--db", db_path,
            "--outfmt", "6", *fields,
            "--out", raw_results_path,
            "--threads", str(num_threads),
            "--quiet",
        ]

        # results stream to --out on disk rather than being captured in memory;
        # only stderr is captured, for error reporting on failure
        with transient_spinner("running DIAMOND blastx", indent="        "):
            result = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print()
            report_message("DIAMOND blastx failed with the following error:")
            report_message(f"{result.stderr}", join = False, color = "none",
                           initial_indent = "      ", subsequent_indent = "      ")
            report_failure("")

        if os.path.getsize(raw_results_path) == 0:
            blast_df = pd.DataFrame(columns=fields)
        else:
            with transient_spinner("loading results", indent="        "):
                blast_df = pd.read_csv(raw_results_path, sep="\t",
                                       header=None, names=fields)

    with transient_spinner("writing search results", indent="        "):
        # DIAMOND doesn't emit the qcov* fields blastn includes; add placeholders
        # so the combined frame has a consistent shape across entry types
        for missing in ("qcovs", "qcovhsp", "qcovus"):
            blast_df[missing] = 0

        blast_df["perc-subj-cov"] = (round((blast_df["length"] / blast_df["slen"]) * 100, 1)
                                     if not blast_df.empty else pd.Series(dtype=float))

        blast_df.to_csv(f"{outputs_dir}/{out_base}-diamond-results.tsv", sep="\t", index=False)

    return blast_df


def filter_blast_results(blast_df, min_perc_id):

    return blast_df[blast_df["pident"] >= min_perc_id].copy()


def update_assembly_summary_table(loci_df, targets_dict, unique_assembly_name, summary_df):

    if loci_df.empty:
        target_counts = {}
    else:
        target_counts = loci_df["target"].value_counts().to_dict()

    new_row = pd.DataFrame(columns=targets_dict.keys(), index=[unique_assembly_name])

    for target in targets_dict.keys():
        new_row.at[unique_assembly_name, target] = target_counts.get(target, 0)

    summary_df = pd.concat([summary_df, new_row])

    return summary_df


def filter_undetected_assembly_targets(summary_df):
    """ filters out targets that were not found (count 0) in any assembly """

    cols_to_drop = summary_df.columns[
        summary_df.apply(lambda col: set(col.astype(str)) <= {"0"})
    ]

    return summary_df.drop(columns=cols_to_drop)


border = "-" * 80
def report_assembly_screen_finished(args, outputs_dir):
    print(f"\n{border}")
    report_message("DONE!", color = "green", trailing_newline = True)
    out_file = f"{args.output_prefix}-summary.tsv"
    out_dir = f"{outputs_dir}/"
    print(f"    Summary table written to: {color_text(out_file, 'green')}")
    print(f"    Additional outputs written in subdirectory: {color_text(out_dir, 'green')}\n")
    print(f"{border}\n")


def report_read_screen_finished(args):
    print(f"\n{border}")
    report_message("DONE!", color = "green", trailing_newline = True)
    out_file = f"{args.output_prefix}-reads-summary.tsv"
    out_dir = f"{args.output_prefix}-mapping/"
    print(f"    Summary table written to: {color_text(out_file, 'green')}")
    print(f"    Mapping info and logs written in subdirectory: {color_text(out_dir, 'green')}\n")
    print(f"{border}\n")


@dataclass
class ReadsRunConfiguration:
    base_output_prefix: str = field(init=False)
    base_output_dir: Path = field(init=False)
    mapping_output_dir: Path = field(init=False)
    log_files_dir: Path = field(init=False)
    targets: str = field(init=False)
    reads_dir: str = field(init=False)
    min_perc_id: float = field(init=False)
    min_perc_cov: float = field(init=False)
    num_cores: int = field(init=False)
    rerun_incomplete: bool = field(init=False)
    dry_run: bool = field(init=False)

    @classmethod
    def from_args(cls, args):
        reads_run_data = cls()
        reads_run_data.populate_read_run_data(args)
        return reads_run_data

    def populate_read_run_data(self, args):
        self.base_output_prefix = Path(args.output_prefix).resolve().name
        self.base_output_dir = Path(args.output_prefix).resolve().parent
        self.mapping_output_dir = self.base_output_dir / f"{args.output_prefix}-mapping"
        self.log_files_dir = self.mapping_output_dir / f"log-files"
        self.targets = Path(args.targets).absolute()
        self.reads_dir = Path(args.reads_dir).absolute()
        self.min_perc_id = args.min_perc_id
        self.min_perc_cov = args.min_perc_cov
        self.num_cores = args.jobs
        self.rerun_incomplete = args.rerun_incomplete
        self.dry_run = args.dry_run

    @property
    def key_value_pairs(self):
        return [f"{key}={str(value)}" for key, value in vars(self).items()]


def run_reads_snakemake(config, reads_dict):
    reads_json = json.dumps(reads_dict)

    cmd = [
        "snakemake",
        "--snakefile", str(get_package_path("smk/ez-screen-reads.smk")),
        "--cores", str(config.num_cores),
        "--printshellcmds",
        "--directory", config.mapping_output_dir,
        "--config", f'reads_json={reads_json}',
        *config.key_value_pairs,
    ]

    if config.dry_run:
        cmd.append("--dry-run")
    if config.rerun_incomplete:
        cmd.append("--rerun-incomplete")

    process = run(cmd)
    if process.returncode != 0:
        message = "Snakemake failed. Hopefully its output above can help you spot why."
        report_failure(message)


def gen_reads_summary_table(input_bam, input_global_dist_tab, outpath,
                            min_perc_id, min_perc_cov):

    ref_read_pids = gen_ref_read_pids(input_bam)

    filtered_df = gen_coverage_filtered_reads_df(input_global_dist_tab, min_perc_cov)

    # making dictionary of ref_name: mean-of-aligned-read-percent-IDs
    mean_pid_dict = {
        ref: np.mean(pids)
        for ref, pids in ref_read_pids.items()
    }
    # making dictionary of ref_name: number-of-aligned-reads
    read_counts_dict = {
        ref: len(pids)
        for ref, pids in ref_read_pids.items()
    }

    # adding mean percent identity column
    filtered_df = filtered_df.assign(mean_perc_id = filtered_df['target'].map(mean_pid_dict))

    # filtering based on mean percent identity
    filtered_df = filtered_df[filtered_df['mean_perc_id'] >= min_perc_id]

    filtered_df['mean_perc_id'] = filtered_df['mean_perc_id'].map("{:.2f}".format)

    # adding number of reads recruited column
    filtered_df = filtered_df.assign(num_reads_recruited = filtered_df['target'].map(read_counts_dict))

    # re-ordering columns
    filtered_df = filtered_df[
        ['target', 'num_reads_recruited', 'detection', 'mean_perc_id']
    ]

    if filtered_df.empty:
        with open(outpath, 'w') as f:
            f.write("No reads successfully mapped to any targets above the set thresholds.\n")
    else:
        filtered_df.to_csv(outpath, sep='\t', index=False)


def gen_ref_read_pids(input_bam):

    with pysam.AlignmentFile(input_bam, "rb") as bam:

        # store read % identity by reference
        ref_read_pids = defaultdict(list)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            try:
                nm = read.get_tag("NM")
            except KeyError:
                continue


            full_aligned_length = sum(
                length for op, length in read.cigartuples
                if op in {0, 1, 2, 7, 8}   # M, I, D, =, X
            )

            pid = (full_aligned_length - nm) / full_aligned_length * 100

            refname = bam.get_reference_name(read.reference_id)
            ref_read_pids[refname].append(pid)

        return ref_read_pids


def gen_coverage_filtered_reads_df(input_global_dist_tab, min_perc_cov):

    detection_df = pd.read_csv(input_global_dist_tab, sep='\t')
    detection_df.columns = ['target', 'depth', 'detection']

    mask = (detection_df['target'] != "total") & (detection_df['depth'] == 1) & (detection_df['detection'] >= min_perc_cov / 100)
    filtered_df = detection_df[mask]
    filtered_df = filtered_df.drop('depth', axis=1)

    return filtered_df


def combine_reads_summary_outputs(samples_output_summaries_dict, output_tsv):

    expected_cols = {"target", "num_reads_recruited", "detection", "mean_perc_id"}
    dfs = []

    for sample, path in samples_output_summaries_dict.items():
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            # skipping files that don't parse
            continue

        # skipping empty tables
        if df.empty:
            continue

        # skipping if columns aren't exactly the expected ones
        if set(df.columns) != expected_cols:
            continue

        # inserting "sample" as first column
        df.insert(0, "sample", sample)
        dfs.append(df)

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        combined.to_csv(output_tsv, sep="\t", index=False)
    else:
        # no valid data
        with open(output_tsv, 'w') as f:
            f.write("No valid read-mappings were found to any targets.\n")


def safe_name(name):
    return name.replace("/", "_").replace("\\", "_")


def gen_contig_summary_table(loci_df, contig_lengths, num_regions_by_contig=None,
                             contig_island_map=None):

    cols = ["contig", "length", "bases_aligned_to_targets",
            "perc_contig_aligned_to_targets", "num_unique_target_hits",
            "num_total_target_hits"]

    if loci_df.empty:
        if num_regions_by_contig is not None:
            cols.insert(5, "num_regions")
        if contig_island_map:
            cols.append("island_ids")
        return pd.DataFrame(columns=cols)

    grouped = loci_df.groupby("contig")

    bases_aligned = {
        contig: merge_intervals(list(zip(g["q_low"], g["q_high"])))
        for contig, g in grouped
    }
    bases_aligned = pd.Series(bases_aligned)

    lengths = pd.Series({c: contig_lengths[c] for c in bases_aligned.index})

    contig_df = pd.DataFrame({
        "length": lengths,
        "bases_aligned_to_targets": bases_aligned,
        "num_unique_target_hits": grouped["target"].nunique(),
        "num_total_target_hits": grouped.size(),
    })

    for int_col in ("length", "bases_aligned_to_targets", "num_unique_target_hits",
                    "num_total_target_hits"):
        contig_df[int_col] = contig_df[int_col].astype(int)

    contig_df["perc_contig_aligned_to_targets"] = round(
        contig_df["bases_aligned_to_targets"] / contig_df["length"] * 100, 1)

    ordered_cols = ["length", "bases_aligned_to_targets", "perc_contig_aligned_to_targets",
                    "num_unique_target_hits", "num_total_target_hits"]

    if num_regions_by_contig is not None:
        contig_df["num_regions"] = pd.Series(num_regions_by_contig)
        contig_df["num_regions"] = contig_df["num_regions"].fillna(0).astype(int)
        ordered_cols = ["length", "bases_aligned_to_targets", "perc_contig_aligned_to_targets",
                        "num_regions", "num_unique_target_hits", "num_total_target_hits"]

    # island cross-reference: which extracted island(s) came from each contig.
    # 'none' for contigs that produced no extracted island.
    if contig_island_map:
        contig_df["island_ids"] = contig_df.index.map(
            lambda c: ";".join(contig_island_map[c]) if c in contig_island_map else "none")
        ordered_cols = ordered_cols + ["island_ids"]

    contig_df = contig_df[ordered_cols]

    contig_df = contig_df.sort_values("perc_contig_aligned_to_targets", ascending=False)
    contig_df.index.name = "contig"
    contig_df = contig_df.reset_index()

    return contig_df


def merge_intervals(intervals):
    """
    given a list of (start, end) tuples, returns total length covered
    after merging overlaps (1-based inclusive coordinates)
    """

    if not intervals:
        return 0

    intervals = sorted(intervals)
    merged = []

    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return sum(end - start + 1 for start, end in merged)


def locus_abuts_contig_edge(q_low, q_high, contig_length, edge_tolerance=100):
    """ True if the locus reaches within edge_tolerance bp of either contig
        terminus -- i.e. it starts at/near position 1 or ends at/near
        contig_length. used to relax the coverage gate for targets that may be
        truncated by a contig boundary. """
    if contig_length is None:
        return False
    abuts_left = q_low <= edge_tolerance
    abuts_right = q_high >= contig_length - edge_tolerance + 1
    return abuts_left or abuts_right


def resolve_assembly_loci(filtered_hits_df, targets_dict, gap_tolerance=200,
                          min_perc_cov=80.0, contig_lengths=None,
                          min_edge_perc_cov=50.0, edge_tolerance=100):
    """
    given pident-passing blast HSPs for one assembly,
    resolve them into distinct target occurrences ("loci") and
    keep only those that clear the coverage gate. every downstream output
    (summary counts, contig table, region-calls) draws from this
    resolution

    per (target, contig):
        1. gap-join HSPs into loci in contig space (two HSPs share a locus if
            they are within gap_tolerance bp on the contig)
        2. for each locus, merge the subject intervals (sstart..send) and
            divide by the target length -> perc_target_cov at that locus
        3. keep the locus if perc_target_cov >= min_perc_cov; OR, for a locus
            abutting a contig edge (within edge_tolerance bp of either end), keep
            it under the relaxed min_edge_perc_cov threshold -- this rescues
            targets truncated by a contig boundary on fragmented assemblies.
            edge rescue needs contig_lengths; without it, only the normal gate
            applies.

    returns a DataFrame with one row per passing locus, columns:
        contig, target, q_low, q_high, perc_target_cov, pident, bitscore,
        length, slen
    where q_low/q_high are the locus's span on the contig, and pident /
    bitscore / length are taken from the locus's single best HSP (highest
    bitscore) so the region-resolution downstream has a representative hit
    to rank
    """

    cols = ["contig", "target", "q_low", "q_high", "perc_target_cov",
            "pident", "bitscore", "length", "slen"]

    if filtered_hits_df.empty:
        return pd.DataFrame(columns=cols)

    df = filtered_hits_df.copy()
    df["q_low"] = df[["qstart", "qend"]].min(axis=1)
    df["q_high"] = df[["qstart", "qend"]].max(axis=1)
    df["s_low"] = df[["sstart", "send"]].min(axis=1)
    df["s_high"] = df[["sstart", "send"]].max(axis=1)

    # contig lengths for the edge-rescue test. the HSP frame already carries
    # qlen per row, so derive them here if not supplied by the caller -- keeps
    # the edge gate self-contained.
    if contig_lengths is None:
        contig_lengths = df.drop_duplicates("qseqid").set_index("qseqid")["qlen"].to_dict()

    # sort once so rows for each (target, contig) are contiguous and ordered
    # along the contig; then walk the arrays in a single linear pass. pandas
    # groupby over tens of thousands of (sseqid, qseqid) pairs -- most of them a
    # single HSP -- is dominated by per-group overhead and .to_dict("records")
    # materialization, so we avoid it entirely and index raw numpy columns.
    df = df.sort_values(["sseqid", "qseqid", "q_low", "q_high"])

    sseqid = df["sseqid"].to_numpy()
    qseqid = df["qseqid"].to_numpy()
    q_low = df["q_low"].to_numpy()
    q_high = df["q_high"].to_numpy()
    s_low = df["s_low"].to_numpy()
    s_high = df["s_high"].to_numpy()
    pident = df["pident"].to_numpy()
    bitscore = df["bitscore"].to_numpy()
    length = df["length"].to_numpy()
    slen = df["slen"].to_numpy()

    n = len(df)
    locus_rows = []

    def flush_locus(idxs, target, contig, target_length):
        # merged subject coverage of the target at this locus
        subj_bases = merge_intervals([(s_low[j], s_high[j]) for j in idxs])
        perc_target_cov = round(subj_bases / target_length * 100, 1)

        locus_q_low = min(q_low[j] for j in idxs)
        locus_q_high = max(q_high[j] for j in idxs)

        # two-tier coverage gate: normal threshold everywhere, OR the relaxed
        # edge threshold for loci abutting a contig terminus (which may be
        # truncated by the boundary rather than genuinely partial)
        if perc_target_cov < min_perc_cov:
            contig_length = contig_lengths.get(contig) if contig_lengths else None
            rescued = (
                perc_target_cov >= min_edge_perc_cov and
                locus_abuts_contig_edge(locus_q_low, locus_q_high,
                                        contig_length, edge_tolerance)
            )
            if not rescued:
                return  # fails both the normal and edge coverage gates

        # representative HSP for the locus = highest bitscore
        best = max(idxs, key=lambda j: bitscore[j])

        locus_rows.append({
            "contig": contig,
            "target": target,
            "q_low": locus_q_low,
            "q_high": locus_q_high,
            "perc_target_cov": perc_target_cov,
            "pident": pident[best],
            "bitscore": bitscore[best],
            "length": length[best],
            "slen": slen[best],
        })

    i = 0
    while i < n:
        target = sseqid[i]
        contig = qseqid[i]

        # advance to the end of this (target, contig) block
        j = i + 1
        while j < n and sseqid[j] == target and qseqid[j] == contig:
            j += 1

        target_length = targets_dict.get(target)
        if target_length:
            # gap-join the block's HSPs (already q_low-ordered by the sort) into
            # loci: a new locus starts when a hit's q_low exceeds the running
            # max q_high by more than gap_tolerance
            cur = [i]
            cur_qhigh = q_high[i]
            for k in range(i + 1, j):
                if q_low[k] > cur_qhigh + gap_tolerance:
                    flush_locus(cur, target, contig, target_length)
                    cur = [k]
                    cur_qhigh = q_high[k]
                else:
                    cur.append(k)
                    if q_high[k] > cur_qhigh:
                        cur_qhigh = q_high[k]
            flush_locus(cur, target, contig, target_length)

        i = j

    if not locus_rows:
        return pd.DataFrame(columns=cols)

    return pd.DataFrame(locus_rows, columns=cols)


def _reciprocal_overlap_frac(a_low, a_high, b_low, b_high):
    """
    fraction of the SHORTER of the two spans that is covered by their
    overlap (1-based inclusive coordinates). returns 0.0 if they don't
    overlap
    """

    overlap = min(a_high, b_high) - max(a_low, b_low) + 1
    if overlap <= 0:
        return 0.0

    len_a = a_high - a_low + 1
    len_b = b_high - b_low + 1

    return overlap / min(len_a, len_b)


def resolve_contig_regions(contig_loci_df, overlap_frac=0.5):
    """
    given the passing loci for an individual contig (a DataFrame with at least
    q_low, q_high, target, bitscore, pident, length, slen, perc_target_cov
    columns -- as produced by resolve_assembly_loci), cluster loci that
    reciprocal-overlap by >= overlap_frac of the shorter span, pick a winner
    per cluster (by bitscore, tie-broken by pident*length), and return one
    region record per cluster

    a locus joins the current cluster if it reciprocal-overlaps ANY locus
    already in that cluster by >= overlap_frac (single left-to-right sweep
    over loci sorted by q_low, then q_high). region coordinates are the
    winner's span; target_length is the winner's subject length (slen);
    perc_target_cov is the winner's coverage

    returns a list of dicts, one per resolved region, ordered by region_start

    this is a thin DataFrame adapter; the resolution logic lives in
    _resolve_contig_regions_from_rows, which gen_region_calls_table calls
    directly with tuple lists to stay off the pandas per-group path
    """

    if contig_loci_df.empty:
        return []

    # tuple layout: (q_low, q_high, target, bitscore, pident, length, slen, perc_target_cov[, target_type])
    has_type = "target_type" in contig_loci_df.columns
    pull_cols = ["q_low", "q_high", "target", "bitscore", "pident",
                 "length", "slen", "perc_target_cov"]
    if has_type:
        pull_cols.append("target_type")
    rows = list(contig_loci_df[pull_cols].itertuples(index=False, name=None))

    return _resolve_contig_regions_from_rows(rows, overlap_frac=overlap_frac,
                                             has_type=has_type)


def _resolve_contig_regions_from_rows(rows, overlap_frac=0.5, has_type=True):
    """
    core region resolution operating on a list of tuples in the layout
        (q_low, q_high, target, bitscore, pident, length, slen,
         perc_target_cov[, target_type])
    see resolve_contig_regions for the algorithm; kept tuple-based so callers
    with many contigs can avoid per-contig DataFrame construction
    """

    if not rows:
        return []

    # sort by start, then end
    rows = sorted(rows, key=lambda r: (r[0], r[1]))

    def is_better(r, best):
        # winner = highest bitscore, tie-break by pident*length
        r_key = (r[3], r[4] * r[5])
        best_key = (best[3], best[4] * best[5])
        return r_key > best_key

    clusters = []
    current = [rows[0]]

    for r in rows[1:]:
        q_low, q_high = r[0], r[1]
        joins = any(
            _reciprocal_overlap_frac(q_low, q_high, m[0], m[1]) >= overlap_frac
            for m in current
        )
        if joins:
            current.append(r)
        else:
            clusters.append(current)
            current = [r]
    clusters.append(current)

    regions = []
    for cluster in clusters:
        winner = cluster[0]
        for r in cluster[1:]:
            if is_better(r, winner):
                winner = r

        # alternates, best-first, excluding the winner (by identity)
        alternates = [r for r in cluster if r is not winner]
        alternates.sort(key=lambda r: (r[3], r[4] * r[5]), reverse=True)

        region_start = int(winner[0])
        region_end = int(winner[1])

        regions.append({
            "region_start": region_start,
            "region_end": region_end,
            "region_length": region_end - region_start + 1,
            "aligned_target": winner[2],
            "target_type": winner[8] if has_type else "nt",
            "target_length": int(winner[6]),
            "pident": float(winner[4]),
            "perc_target_cov": float(winner[7]),
            "bitscore": float(winner[3]),
            "n_overlapping_targets": len(cluster),
            "other_targets": [r[2] for r in alternates],
        })

    regions.sort(key=lambda d: d["region_start"])

    return regions


natural_sort_key = natsort_keygen()

def gen_region_calls_table(loci_df, overlap_frac=0.5, contig_lengths=None):
    """
    builds a per-resolved-region table across all contigs from the resolved
    loci (resolve_assembly_loci): one row per region, with the chosen (best)
    target reported as 'aligned_target', the region's coordinates and length
    on the contig, the aligned target's length, its pident, its per-locus
    coverage (perc_target_cov), bitscore, and any other targets that
    overlapped the same region folded into 'other_targets'
    (semicolon-delimited; "NA" if none)

    perc_target_cov is the winner's coverage only

    columns: contig, region_start, region_end, region_length,
    aligned_target, target_type, target_length, pident, perc_target_cov,
    bitscore, n_overlapping_targets, other_targets. target_type ('nt'/'aa')
    is the units key for target_length -- for 'aa' targets the length is in
    amino acids

    also returns a {contig: num_regions} dict so the contig summary can
    report a region count sourced from the exact same resolution, keeping
    the two outputs from drifting
    """

    cols = ["contig", "contig_length", "region_start", "region_end", "region_length",
            "aligned_target", "target_type", "target_length", "pident",
            "perc_target_cov", "bitscore", "n_overlapping_targets", "other_targets"]

    if loci_df.empty:
        return pd.DataFrame(columns=cols), {}

    region_rows = []
    num_regions_by_contig = {}

    has_type = "target_type" in loci_df.columns
    pull_cols = ["contig", "q_low", "q_high", "target", "bitscore", "pident",
                 "length", "slen", "perc_target_cov"]
    if has_type:
        pull_cols.append("target_type")

    sorted_loci = loci_df[pull_cols].sort_values("contig", kind="stable")
    records = list(sorted_loci.itertuples(index=False, name=None))
    # records tuple layout: (contig, q_low, q_high, target, bitscore, pident,
    #                        length, slen, perc_target_cov[, target_type])

    n = len(records)
    i = 0
    while i < n:
        contig = records[i][0]
        j = i + 1
        while j < n and records[j][0] == contig:
            j += 1

        # drop the leading contig field to match the row layout
        contig_rows = [rec[1:] for rec in records[i:j]]
        regions = _resolve_contig_regions_from_rows(contig_rows,
                                                    overlap_frac=overlap_frac,
                                                    has_type=has_type)
        num_regions_by_contig[contig] = len(regions)
        for r in regions:
            region_rows.append({
                "contig": contig,
                "contig_length": (contig_lengths or {}).get(contig, ""),
                "region_start": r["region_start"],
                "region_end": r["region_end"],
                "region_length": r["region_length"],
                "aligned_target": r["aligned_target"],
                "target_type": r["target_type"],
                "target_length": r["target_length"],
                "pident": round(r["pident"], 1),
                "perc_target_cov": r["perc_target_cov"],
                "bitscore": r["bitscore"],
                "n_overlapping_targets": r["n_overlapping_targets"],
                "other_targets": ";".join(r["other_targets"]) if r["other_targets"] else "none",
            })

        i = j

    if not region_rows:
        return pd.DataFrame(columns=cols), num_regions_by_contig

    region_df = pd.DataFrame(region_rows, columns=cols)
    region_df = region_df.sort_values(
        ["contig", "region_start"],
        key=lambda col: col.map(natural_sort_key) if col.name == "contig" else col
    ).reset_index(drop=True)

    return region_df, num_regions_by_contig


def cluster_regions_into_islands(region_df_for_contig, island_gap=2500):
    """
    given the resolved regions for an individual contig (rows with region_start /
    region_end and a DataFrame index), this chains regions within island_gap bp of
    each other into islands

    returns a list of dicts:
        {island_start, island_end, num_regions, region_indices}
    where region_indices are the DataFrame index labels of the member
    regions (used to tag region-calls with the island id)
    """

    if region_df_for_contig.empty:
        return []

    rows = sorted(
        ((int(r.region_start), int(r.region_end), idx)
         for idx, r in region_df_for_contig.iterrows()),
        key=lambda t: (t[0], t[1]))

    islands = []
    cur_start, cur_end, idxs, cur_count = rows[0][0], rows[0][1], [rows[0][2]], 1

    for start, end, idx in rows[1:]:
        if start <= cur_end + island_gap:
            cur_end = max(cur_end, end)
            cur_count += 1
            idxs.append(idx)
        else:
            islands.append({"island_start": cur_start, "island_end": cur_end,
                            "num_regions": cur_count, "region_indices": idxs})
            cur_start, cur_end, idxs, cur_count = start, end, [idx], 1
    islands.append({"island_start": cur_start, "island_end": cur_end,
                    "num_regions": cur_count, "region_indices": idxs})

    return islands


def island_passes_extraction(island, contig_length, contig_floor=20000,
                             ratio=3.0, min_span=2000, min_regions=3):
    """
    four-part extraction trigger on one island (base span (no buffer)):
        contig_length >= contig_floor          (big enough to bother cutting)
        contig_length >= ratio * island_span   (island is a minority of contig)
        island_span   >= min_span              (long enough to cut)
        num_regions   >= min_regions           (a real cluster, not a lone hit)
    """
    span = island["island_end"] - island["island_start"] + 1
    return (
        contig_length >= contig_floor and
        contig_length >= ratio * span and
        span >= min_span and
        island["num_regions"] >= min_regions
    )


def buffered_island_coords(island, contig_length, buffer=500):
    """
    adds the flanking buffer to an island's base coords, up
    to any contig ends at least
    """
    bstart = max(1, island["island_start"] - buffer)
    bend = min(contig_length, island["island_end"] + buffer)
    return bstart, bend


def make_island_id(contig, bstart, bend):
    """
    makes the identifier shared by the island's fasta filename, its
    manifest row, and the cross-reference columns in the other tables
    """
    return f"{safe_name(contig)}_island_{bstart}-{bend}"


def extract_islands(region_df, contig_lengths, assembly_path,
                    outputs_dir, out_base,
                    island_gap=5000, contig_floor=20000,
                    ratio=3.0, min_span=2000, min_regions=3,
                    buffer=500):
    """
    clusters resolved regions into islands, keeping those passing the extraction
    trigger, slices each (buffered) out of its contig, and writes one
    fasta per island into a per-assembly subdir plus a manifest tsv

    returns:
        region_island_map: {region_df_index: island_id} for passing islands
        contig_island_map: {contig: [island_id, ...]}
    so we can add island_id / island_ids columns to the other tables

    sequences are read lazily: only contigs that have a passing island are
    pulled from the input assembly(ies). islands require resolved regions, so if
    region_df is empty (e.g., --dont-resolve-regions) nothing is extracted
    """

    manifest_cols = ["island_id", "contig", "island_start", "island_end",
                     "island_span", "buffered_start", "buffered_end",
                     "num_regions", "contig_length"]
    manifest_path = f"{outputs_dir}/{out_base}-extracted-islands.tsv"

    region_island_map = {}
    contig_island_map = {}

    if region_df is None or region_df.empty:
        pd.DataFrame(columns=manifest_cols).to_csv(manifest_path, sep="\t", index=False)
        return region_island_map, contig_island_map

    # 1. cluster per contig and collect passing islands (no sequence reading yet)
    passing = []  # (contig, island dict, bstart, bend)
    for contig, g in region_df.groupby("contig"):
        contig_length = contig_lengths.get(contig)
        if not contig_length:
            continue
        for island in cluster_regions_into_islands(g, island_gap=island_gap):
            if island_passes_extraction(island, contig_length,
                                        contig_floor=contig_floor, ratio=ratio,
                                        min_span=min_span, min_regions=min_regions):
                bstart, bend = buffered_island_coords(island, contig_length, buffer=buffer)
                passing.append((contig, island, bstart, bend))

    if not passing:
        pd.DataFrame(columns=manifest_cols).to_csv(manifest_path, sep="\t", index=False)
        return region_island_map, contig_island_map

    # 2. lazily read only the contigs we need from the assembly
    needed_contigs = {contig for contig, *_ in passing}
    contig_seqs = {}
    for record in SeqIO.parse(assembly_path, "fasta"):
        if record.id in needed_contigs:
            contig_seqs[record.id] = record.seq
            if len(contig_seqs) == len(needed_contigs):
                break

    # 3. slice, write per-island FASTAs into the subdir, build manifest + maps
    island_subdir = f"{outputs_dir}/{out_base}-extracted-islands"
    os.makedirs(island_subdir, exist_ok=True)

    manifest_rows = []
    for contig, island, bstart, bend in passing:
        seq = contig_seqs.get(contig)
        if seq is None:
            continue  # contig not found in assembly (shouldn't happen); skip safely
        island_id = make_island_id(contig, bstart, bend)
        subseq = seq[bstart - 1:bend]  # 1-based inclusive -> python slice

        with open(f"{island_subdir}/{island_id}.fasta", "w") as f:
            f.write(f">{island_id}\n{str(subseq)}\n")

        for idx in island["region_indices"]:
            region_island_map[idx] = island_id
        contig_island_map.setdefault(contig, []).append(island_id)

        manifest_rows.append({
            "island_id": island_id,
            "contig": contig,
            "island_start": island["island_start"],
            "island_end": island["island_end"],
            "island_span": island["island_end"] - island["island_start"] + 1,
            "buffered_start": bstart,
            "buffered_end": bend,
            "num_regions": island["num_regions"],
            "contig_length": contig_lengths[contig],
        })

    pd.DataFrame(manifest_rows, columns=manifest_cols).to_csv(
        manifest_path, sep="\t", index=False)

    return region_island_map, contig_island_map


def check_sseqid_target_match(blast_df, target_names):
    """
    guards against the silent-empty-output failure where BLAST/DIAMOND sseqids
    don't match any target name (e.g. a db built without -parse_seqids). In that
    case every coverage lookup would miss and all loci would get dropped,
    producing an empty run with no error

    target_names is the set of target names for the entry being checked (already
    collision-suffixed to match the remapped sseqids), so per-entry checking
    doesn't false-trip against other entries' names
    """

    if blast_df.empty:
        return  # genuine no-hits; nothing to diagnose

    blast_sseqids = set(blast_df["sseqid"].unique())
    target_keys = set(target_names)

    if blast_sseqids & target_keys:
        return  # at least some overlap; ids line up

    # no overlap at all -> guaranteed every locus would be dropped
    example_sseqid = next(iter(blast_sseqids))
    example_target = next(iter(target_keys)) if target_keys else "(none)"

    msg = (
        "A target search returned hits, but none of the subject IDs don't match the headers "
        "pulled out in fasta form from the database. (This won't work for us...)")

    print()
    report_message(msg)

    report_message(f"Example subject ID from search:  {example_sseqid}",
                   initial_indent = "      ", subsequent_indent = "      ",
                   color = "none")
    report_message(f"Example header pulled out:       {example_target}",
                   initial_indent = "      ", subsequent_indent = "      ",
                   color = "none", leading_newline = False, join = False)

    msg = (
        "For a BLAST db this usually means it was built without "
        "'-parse_seqids'. Rebuild with '-parse_seqids', or pass the original "
        "fasta to `--targets` instead."
    )

    report_message(msg)

    report_failure()
