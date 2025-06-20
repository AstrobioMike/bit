import os
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
from tqdm import tqdm
from bit.utils import (report_message,
                       report_failure,
                       get_input_reads_dict)


def run_assembly(args):
    blast_results_dir = args.output_prefix + "-blast-results"
    assembly_path_dict = assembly_preflight(args, blast_results_dir)
    run_assembly_screen(args, assembly_path_dict, blast_results_dir)
    report_assembly_screen_finished(args, blast_results_dir)


def run_reads(args):
    reads_dict = get_input_reads_dict(args.reads_dir)
    print(reads_dict)
    print(reads_dict['perfect-reads']['R1'])


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

    if args.transpose_output_tsv:
        summary_df = summary_df.T
        summary_df.to_csv(args.output_prefix + "-summary-table.tsv", sep = "\t", index_label = "target")
    else:
        summary_df.to_csv(args.output_prefix + "-summary-table.tsv", sep = "\t", index_label = "input-assembly")


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


def report_assembly_screen_finished(args, blast_results_dir):

    report_message("Done!", color = "green")
    print(f"    Summary table written to: {args.output_prefix}-summary-table.tsv")
    print(f"    Full and filtered BLAST results written in subdirectory: {blast_results_dir}/\n")
