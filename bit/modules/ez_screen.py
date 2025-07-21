import os
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
from tqdm import tqdm
from dataclasses import dataclass, field
from pathlib import Path
import json
from subprocess import run
import pysam
from collections import defaultdict
import numpy as np
from bit.modules.input_parsing import get_input_reads_dict_from_dir
from bit.modules.general import (report_message,
                       report_failure,
                       get_package_path,
                       color_text,
                       log_command_run)


def run_assembly(args, full_cmd_executed):
    blast_results_dir = args.output_prefix + "-blast-results"
    assembly_path_dict = assembly_preflight(args, blast_results_dir)
    run_assembly_screen(args, assembly_path_dict, blast_results_dir)
    report_assembly_screen_finished(args, blast_results_dir)
    log_command_run(full_cmd_executed, blast_results_dir)


def run_reads(args, full_cmd_executed):
    reads_dict = get_input_reads_dict_from_dir(args.reads_dir)
    reads_config = ReadsRunConfiguration.from_args(args)
    run_reads_snakemake(reads_config, reads_dict)
    report_read_screen_finished(args)
    log_command_run(full_cmd_executed, reads_config.log_files_dir)


def assembly_preflight(args, blast_results_dir):

    check_assembly_inputs(args.assemblies, args.targets)

    if not os.path.exists(blast_results_dir):
        os.makedirs(blast_results_dir)

    # checking if there'd be any duplicates with just basenames (like same filename from different path),
    # and retaining input path info if so
    assembly_basenames = [os.path.splitext(os.path.basename(assembly))[0] for assembly in args.assemblies]
    basename_counts = {basename: assembly_basenames.count(basename) for basename in assembly_basenames}

    if any(count > 1 for count in basename_counts.values()):
        assembly_path_dict = {assembly: assembly for assembly in args.assemblies}
    else:
        assembly_path_dict = {assembly: os.path.splitext(os.path.basename(assembly))[0] for assembly in args.assemblies}

    return assembly_path_dict


def check_assembly_inputs(assemblies, targets):
    """ checks that input files exist """

    for assembly in assemblies:
        if not os.path.exists(assembly):
            report_failure(f"Specified input assembly file not found: {assembly}")
        if os.path.isdir(assembly):
            report_failure(f"Specified input assembly is a directory, but needs to be a file or files: {assembly}")

    if not os.path.exists(targets):
        report_failure(f"Specified input targets file not found: {targets}")
    if os.path.isdir(targets):
        report_failure(f"Specified input targets is a directory, but needs to be a file: {assembly}")



def run_assembly_screen(args, assembly_path_dict, blast_results_dir):

    targets_dict = get_targets(args.targets)
    summary_df = pd.DataFrame()
    total_assemblies = len(args.assemblies)

    print("")

    for i, assembly in enumerate(tqdm(args.assemblies, total = total_assemblies,
                                      desc = "Processing assemblies", unit = "assembly",
                                      bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} processed, {remaining} remaining")):

        blast_df = run_blast(assembly, args.targets, blast_results_dir)
        filtered_blast_df, long_targets_results_dict = filter_blast_results(blast_df, targets_dict,
                                                                            args.min_perc_id,
                                                                            args.min_perc_cov)

        unique_assembly_name = assembly_path_dict[assembly]
        summary_df = update_assembly_summary_table(filtered_blast_df,
                                          long_targets_results_dict,
                                          targets_dict,
                                          unique_assembly_name,
                                          summary_df)

    if args.filter_if_not_detected:
        summary_df = filter_undetected_assembly_targets(summary_df)

    output_tsv = args.output_prefix + "-assembly-summary.tsv"
    if args.transpose_output_tsv:
        summary_df = summary_df.T
        summary_df.to_csv(output_tsv, sep = "\t", index_label = "target")
    else:
        summary_df.to_csv(output_tsv, sep = "\t", index_label = "input-assembly")


def get_targets(targets):
    """ creates a dictionary of the targets based on input fasta headers as keys and lengths as values"""

    targets_dict = {}

    for record in SeqIO.parse(targets, "fasta"):
        targets_dict[record.id] = len(record.seq)

    return targets_dict


def run_blast(assembly, targets, blast_results_dir):
    """ runs BLAST to search for targets in assembly """

    blast_command = [
        "blastn",
        "-task", "blastn",
        "-query", assembly,
        "-subject", targets,
        "-outfmt", "6 qseqid qlen sseqid slen qstart qend sstart send length qcovs qcovhsp qcovus pident evalue bitscore",
    ]

    try:
        blast_results = subprocess.check_output(blast_command, stderr = subprocess.STDOUT, text = True)
    except subprocess.CalledProcessError as e:
        report_failure("BLAST failed with the following error:  \n" + e.output)

    blast_results = StringIO(blast_results)

    blast_df = pd.read_csv(blast_results, sep = "\t", header = None, names = [
        "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart",
        "send", "length", "qcovs", "qcovhsp", "qcovus", "pident", "evalue", "bitscore"
    ])

    assembly_base = os.path.splitext(os.path.basename(assembly))[0]

    # adding percent of subject covered by alignment
    blast_df["perc-subj-cov"] = round((blast_df["length"] / blast_df["slen"]) * 100, 1)

    blast_df.to_csv(f"{blast_results_dir}/{assembly_base}-blast-results.tsv", sep = "\t", index = False)

    return blast_df


def filter_blast_results(blast_df, targets_dict, min_perc_id, min_perc_cov):
    """ Filter BLAST results based on user-specified thresholds """

    long_target_length_cutoff = 10000

    # if query is targets and subject is assembly
    # short_targets_df, long_targets_df = blast_df[blast_df["qlen"] < long_target_length_cutoff], blast_df[blast_df["qlen"] >= long_target_length_cutoff]

    # if query is assembly and subject is targets
    short_targets_df, long_targets_df = blast_df[blast_df["slen"] < long_target_length_cutoff], blast_df[blast_df["slen"] >= long_target_length_cutoff]

    filtered_short_blast_df = short_targets_df[(short_targets_df["pident"] >= min_perc_id) & (short_targets_df["perc-subj-cov"] >= min_perc_cov)]

    if not long_targets_df.empty:

        filtered_long_blast_df = long_targets_df[(long_targets_df["pident"] >= min_perc_id)]
        long_targets_results_dict = gen_long_targets_assembly_results_dict(filtered_long_blast_df,
                                                                           targets_dict,
                                                                           long_target_length_cutoff,
                                                                           min_perc_cov)

    else:
        long_targets_results_dict = None

    return filtered_short_blast_df, long_targets_results_dict


def update_assembly_summary_table(filtered_short_blast_df, long_targets_detected_dict, targets_dict, unique_assembly_name, summary_df):
    """ Update the summary table with results from the current assembly """

    target_counts = filtered_short_blast_df["sseqid"].value_counts()

    new_row = pd.DataFrame(columns=targets_dict.keys(), index=[unique_assembly_name])

    for target in targets_dict.keys():
        new_row.at[unique_assembly_name, target] = target_counts.get(target, 0)

    if long_targets_detected_dict:
        for col, new_value in long_targets_detected_dict.items():
            new_row.at[unique_assembly_name, col] = new_value

    summary_df = pd.concat([summary_df, new_row])

    return summary_df


def filter_undetected_assembly_targets(summary_df):
    """ filters out targets that weren't detected in any input assemblies """

    cols_to_drop = summary_df.columns[
        summary_df.apply(lambda col: set(col.astype(str)) <= {"0", "NOT-DETECTED"})
    ]

    return summary_df.drop(columns=cols_to_drop)


def gen_long_targets_assembly_results_dict(filtered_long_blast_df, targets_dict, long_target_length_cutoff, min_perc_cov):
    """ remove overlapping alignments and calculate total percent of subject covered by alignment """

    # the sstart and sstop positions aren't necessarily in order from lowest to highest,
    # which is needed for how we merge the covered regions below
    # so here was are making new columns that are the min and max of the two
    filtered_long_blast_df = filtered_long_blast_df.copy() # make a copy to avoid SettingWithCopyWarning
    filtered_long_blast_df["mod_sstart"] = filtered_long_blast_df[["sstart", "send"]].min(axis = 1)
    filtered_long_blast_df["mod_send"] = filtered_long_blast_df[["sstart", "send"]].max(axis = 1)

    long_subject_ids = [key for key, value in targets_dict.items() if value >= long_target_length_cutoff]

    long_targets_detected_dict = {}

    for ID in long_subject_ids:

        sub_df = filtered_long_blast_df[filtered_long_blast_df["sseqid"] == ID]

        if sub_df.empty:
            long_targets_detected_dict[ID] = "NOT-DETECTED"
            continue

        sub_df = sub_df.sort_values(by = ["mod_sstart"])

        merged_regions = []

        for index, row in sub_df.iterrows():

            start, end = row['mod_sstart'], row['mod_send']

            # if this is the first region or no overlap, adding the region
            if not merged_regions or start > merged_regions[-1][1]:
                merged_regions.append([start, end])
            else:
                # otherwise, merging the current region with the previous one
                merged_regions[-1][1] = max(merged_regions[-1][1], end)

        total_subject_bases_covered = sum(end - start + 1 for start, end in merged_regions)

        if total_subject_bases_covered >= min_perc_cov * targets_dict[ID] / 100:
            long_targets_detected_dict[ID] = "DETECTED"
        else:
            long_targets_detected_dict[ID] = "NOT-DETECTED"

    return long_targets_detected_dict


border = "-" * 80
def report_assembly_screen_finished(args, blast_results_dir):
    print(f"\n{border}")
    report_message("DONE!", color = "green")
    out_file = f"{args.output_prefix}-assembly-summary.tsv"
    out_dir = f"{blast_results_dir}/"
    print(f"    Summary table written to: {color_text(out_file, 'green')}")
    print(f"    Full and filtered BLAST results written in subdirectory: {color_text(out_dir, 'green')}\n")
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
                nm = read.get_tag("NM")  # edit distance
            except KeyError:
                continue

            matches = read.query_length - nm
            pid = matches / read.query_length * 100
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
