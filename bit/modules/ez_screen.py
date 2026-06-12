import os
import shutil
import pandas as pd # type: ignore
import subprocess
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
                       log_command_run)


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
        ("fasta", "nt")  nucleotide fasta   -> blastn -subject
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
        report_failure(
            f"--targets entry looks like a protein BLAST database, which isn't "
            f"supported:\n    {entry}\n\n  Pass a protein fasta or a DIAMOND "
            f"database (built with `diamond makedb`) instead.")

    report_failure(
        f"--targets entry is neither a readable fasta nor a recognized BLAST or "
        f"DIAMOND database:\n    {entry}")


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

    print("")

    for i, assembly in enumerate(tqdm(args.assemblies, total = total_assemblies,
                                      desc = "Processing assemblies", unit = "assembly", ncols=76,
                                      bar_format="    {l_bar}{bar} | {n_fmt}/{total_fmt} processed, {remaining} remaining")):

        unique_assembly_name = assembly_path_dict[assembly]
        out_base = safe_name(unique_assembly_name)

        # ---- per-targets-entry search + per-type loci resolution ----
        # each entry is homogeneous (one modality, one seqtype), so it gets the
        # right search engine and the right pident threshold; resulting loci are
        # tagged with target_type and concatenated, so all downstream outputs
        # see one merged, type-aware loci frame.
        per_entry_filtered = []
        per_entry_loci = []

        for entry_info in targets_plan["entries"]:
            entry = entry_info["entry"]
            seqtype = entry_info["seqtype"]
            min_perc_id = args.min_aa_perc_id if seqtype == "aa" else args.min_nt_perc_id

            blast_df = run_targets_search(assembly, entry_info, outputs_dir,
                                          out_base, args.num_threads)

            # remap this entry's sseqids into the merged (possibly suffixed)
            # keyspace, so coverage lookups and downstream joins line up
            blast_df = apply_collision_suffix(blast_df, entry_info)

            # guard against id mismatch (e.g. a db built without -parse_seqids):
            # check against THIS entry's target names only
            check_sseqid_target_match(blast_df, entry_info["names"])

            filtered = filter_blast_results(blast_df, min_perc_id)
            per_entry_filtered.append(filtered)

            loci = resolve_assembly_loci(filtered, targets_dict,
                                         gap_tolerance=args.hit_merge_gap,
                                         min_perc_cov=args.min_perc_cov)
            if not loci.empty:
                loci = loci.assign(target_type=loci["target"].map(type_map))
            else:
                loci = loci.assign(target_type=pd.Series(dtype=object))
            per_entry_loci.append(loci)

        # merge across entries
        filtered_hits_df = (pd.concat(per_entry_filtered, ignore_index=True)
                            if per_entry_filtered else pd.DataFrame())
        loci_df = (pd.concat(per_entry_loci, ignore_index=True)
                   if per_entry_loci else pd.DataFrame())

        # per-assembly filtered BLAST table (pident-passing HSPs, all entries)
        filtered_hits_df.to_csv(f"{outputs_dir}/{out_base}-filtered-blast-results.tsv",
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
        if args.resolve_regions:
            region_df, num_regions_by_contig = gen_region_calls_table(
                loci_df, overlap_frac=args.region_overlap_frac)
            region_df.to_csv(f"{outputs_dir}/{out_base}-region-calls.tsv",
                             sep="\t", index=False)
        else:
            num_regions_by_contig = None

        # per-assembly contig summary table
        contig_df = gen_contig_summary_table(loci_df, contig_lengths,
                                             num_regions_by_contig)
        contig_df.to_csv(f"{outputs_dir}/{out_base}-hit-contig-summary.tsv",
                         sep="\t", index=False)

        summary_df = update_assembly_summary_table(loci_df,
                                          targets_dict,
                                          unique_assembly_name,
                                          summary_df)

    if not args.report_all_targets:
        summary_df = filter_undetected_assembly_targets(summary_df)

    output_tsv = args.output_prefix + "-summary.tsv"
    # internally summary_df is built assemblies-as-rows (one row appended per
    # assembly). the default OUTPUT orientation is targets-as-rows, which stays
    # readable when targets greatly outnumber assemblies
    # and matches the feature-as-row convention of presence/absence tables.
    # --assemblies-as-rows restores the older assemblies-as-rows orientation.
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
        ["diamond", "getseq", "--db", targets],
        check=True, capture_output=True, text=True)
    targets_dict = {}
    for record in SeqIO.parse(StringIO(result.stdout), "fasta"):
        targets_dict[record.id] = len(record.seq)
    return targets_dict


def run_targets_search(assembly, entry_info, outputs_dir, out_base, num_threads):
    """
    dispatches one targets entry to the right search engine and return a
    BLAST-outfmt-6-shaped DataFrame. nucleotide -> blastn; protein -> DIAMOND
    blastx (building a .dmnd in outputs_dir first if the entry is a fasta)
    """

    entry = entry_info["entry"]
    modality = entry_info["modality"]
    seqtype = entry_info["seqtype"]
    # tag outputs per-entry so multiple targets files don't clobber each other
    entry_tag = safe_name(os.path.splitext(os.path.basename(entry))[0])

    if seqtype == "nt":
        return run_blastn(assembly, entry, modality, outputs_dir,
                          f"{out_base}-{entry_tag}", num_threads)

    # protein
    if modality == "fasta":
        db_path = f"{outputs_dir}/{entry_tag}.dmnd"
        if not os.path.exists(db_path):
            build_diamond_db(entry, db_path)
    else:
        db_path = entry  # already a .dmnd db
    return run_diamond_blastx(assembly, db_path, outputs_dir,
                              f"{out_base}-{entry_tag}", num_threads)


def run_blastn(assembly, targets, modality, outputs_dir, out_base, num_threads):

    outfmt = ("6 qseqid qlen sseqid slen qstart qend sstart send length "
              "qcovs qcovhsp qcovus pident evalue bitscore")

    blast_command = ["blastn", "-task", "blastn", "-query", assembly, "-outfmt", outfmt]

    if modality == "db":
        blast_command += ["-db", targets, "-num_threads", str(num_threads)]
    else:
        blast_command += ["-subject", targets]

    try:
        blast_results = subprocess.check_output(blast_command, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        report_failure("BLAST failed with the following error:  \n" + e.output)

    cols = ["qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart",
            "send", "length", "qcovs", "qcovhsp", "qcovus", "pident", "evalue", "bitscore"]
    if blast_results.strip() == "":
        blast_df = pd.DataFrame(columns=cols)
    else:
        blast_df = pd.read_csv(StringIO(blast_results), sep="\t", header=None, names=cols)

    blast_df["perc-subj-cov"] = (round((blast_df["length"] / blast_df["slen"]) * 100, 1)
                                 if not blast_df.empty else pd.Series(dtype=float))

    blast_df.to_csv(f"{outputs_dir}/{out_base}-blast-results.tsv", sep="\t", index=False)

    return blast_df


def build_diamond_db(protein_fasta, db_path):

    cmd = ["diamond", "makedb", "--in", protein_fasta, "--db", db_path, "--quiet"]
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        report_failure("DIAMOND makedb failed with the following error:  \n" + e.output)


def run_diamond_blastx(assembly, db_path, outputs_dir, out_base, num_threads):

    fields = ["qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart",
              "send", "length", "pident", "evalue", "bitscore"]

    cmd = [
        "diamond", "blastx",
        "--query", assembly,
        "--db", db_path,
        "--outfmt", "6", *fields,
        "--threads", str(num_threads),
        "--quiet",
    ]

    try:
        results = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        report_failure("DIAMOND blastx failed with the following error:  \n" + e.output)

    if results.strip() == "":
        blast_df = pd.DataFrame(columns=fields)
    else:
        blast_df = pd.read_csv(StringIO(results), sep="\t", header=None, names=fields)

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
    report_message("DONE!", color = "green")
    out_file = f"{args.output_prefix}-summary.tsv"
    out_dir = f"{outputs_dir}/"
    print(f"    Summary table written to: {color_text(out_file, 'green')}")
    print(f"    Additional outputs written in subdirectory: {color_text(out_dir, 'green')}\n")
    print(f"{border}\n")


def report_read_screen_finished(args):
    print(f"\n{border}")
    report_message("DONE!", color = "green")
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


def gen_contig_summary_table(loci_df, contig_lengths, num_regions_by_contig=None):

    cols = ["contig", "length", "bases_aligned_to_targets",
            "perc_contig_aligned_to_targets", "num_unique_target_hits",
            "num_total_target_hits"]

    if loci_df.empty:
        if num_regions_by_contig is not None:
            cols.insert(5, "num_regions")
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


def resolve_assembly_loci(filtered_hits_df, targets_dict, gap_tolerance=200,
                          min_perc_cov=80.0):
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
        3. keep the locus only if perc_target_cov >= min_perc_cov

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

    locus_rows = []

    for (target, contig), g in df.groupby(["sseqid", "qseqid"]):

        target_length = targets_dict.get(target)
        if not target_length:
            continue

        # order HSPs along the contig, then gap-join into loci
        hsps = g.sort_values(["q_low", "q_high"]).to_dict("records")

        current = [hsps[0]]
        cur_qhigh = hsps[0]["q_high"]
        groups = []

        for h in hsps[1:]:
            if h["q_low"] > cur_qhigh + gap_tolerance:
                groups.append(current)
                current = [h]
                cur_qhigh = h["q_high"]
            else:
                current.append(h)
                cur_qhigh = max(cur_qhigh, h["q_high"])
        groups.append(current)

        for locus in groups:
            # merged subject coverage of the target at this locus
            subj_bases = merge_intervals([(h["s_low"], h["s_high"]) for h in locus])
            perc_target_cov = round(subj_bases / target_length * 100, 1)

            if perc_target_cov < min_perc_cov:
                continue  # hard coverage gate, same treatment as pident filter

            # representative HSP for the locus = highest bitscore
            best = max(locus, key=lambda h: h["bitscore"])

            locus_rows.append({
                "contig": contig,
                "target": target,
                "q_low": min(h["q_low"] for h in locus),
                "q_high": max(h["q_high"] for h in locus),
                "perc_target_cov": perc_target_cov,
                "pident": best["pident"],
                "bitscore": best["bitscore"],
                "length": best["length"],
                "slen": best["slen"],
            })

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

    returns a list of dicts, one per resolved region, ordered by
    region_start
    """

    if contig_loci_df.empty:
        return []

    # tuple layout: (q_low, q_high, target, bitscore, pident, length, slen, perc_target_cov, target_type)
    has_type = "target_type" in contig_loci_df.columns
    pull_cols = ["q_low", "q_high", "target", "bitscore", "pident",
                 "length", "slen", "perc_target_cov"]
    if has_type:
        pull_cols.append("target_type")
    rows = list(contig_loci_df[pull_cols].itertuples(index=False, name=None))

    # sort by start, then end
    rows.sort(key=lambda r: (r[0], r[1]))

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

def gen_region_calls_table(loci_df, overlap_frac=0.5):
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

    cols = ["contig", "region_start", "region_end", "region_length",
            "aligned_target", "target_type", "target_length", "pident",
            "perc_target_cov", "bitscore", "n_overlapping_targets", "other_targets"]

    if loci_df.empty:
        return pd.DataFrame(columns=cols), {}

    region_rows = []
    num_regions_by_contig = {}

    for contig, g in loci_df.groupby("contig"):
        regions = resolve_contig_regions(g, overlap_frac=overlap_frac)
        num_regions_by_contig[contig] = len(regions)
        for r in regions:
            region_rows.append({
                "contig": contig,
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

    if not region_rows:
        return pd.DataFrame(columns=cols), num_regions_by_contig

    region_df = pd.DataFrame(region_rows, columns=cols)
    region_df = region_df.sort_values(
        ["contig", "region_start"],
        key=lambda col: col.map(natural_sort_key) if col.name == "contig" else col
    ).reset_index(drop=True)

    return region_df, num_regions_by_contig


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
        "A targets-search returned hits, but none of the subject IDs match the "
        "target IDs used for coverage gating, so every hit would be silently "
        "dropped.\n\n"
        f"    Example subject ID from search:  {example_sseqid}\n"
        f"    Example target ID expected:      {example_target}\n"
    )

    msg += (
        "\n  For a BLAST db this usually means it was built without "
        "'-parse_seqids'. Rebuild with '-parse_seqids', or pass the original "
        "fasta to --targets instead."
    )

    report_failure(msg)
