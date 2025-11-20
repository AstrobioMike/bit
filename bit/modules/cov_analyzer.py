import os
import subprocess
from pathlib import Path
import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from colorama import Fore, init
from bit.modules.general import (check_files_are_found,
                                 notify_premature_exit,
                                 log_command_run,
                                 tee)

init(autoreset=True)


def run_cov_analyzer(
    reference_fasta: str,
    bam_file: str,
    output_dir: str,
    high_threshold: float,
    low_threshold: float,
    min_region_length: int,
    exclude_contigs: list,
    per_contig: bool,
    sliding_window_size: int,
    step_size: int,
    allowed_gap: int,
    buffer: int,
    no_window_stats: bool,
    log_file = str,
    full_cmd_executed: str = None
):

    preflight_checks(reference_fasta, bam_file, output_dir, exclude_contigs,
                     log_file, full_cmd_executed)

    contig_lengths = get_contig_lengths(reference_fasta)

    window_bed_path = generate_sliding_bed_file(reference_fasta,
                                                sliding_window_size,
                                                step_size,
                                                output_dir)

    mosdepth_regions_file = run_mosdepth(bam_file, window_bed_path, output_dir)

    cov_df = read_mosdepth_regions_file(mosdepth_regions_file)
    if exclude_contigs:
        cov_df = cov_df.loc[~cov_df["contig"].isin(exclude_contigs)].reset_index(drop=True)

    cov_stats = CoverageStats(cov_df, per_contig)

    high_merged_regions = cov_stats.merge_windows(type="high", threshold = high_threshold, contig_lengths = contig_lengths, allowed_gap = allowed_gap)
    low_merged_regions = cov_stats.merge_windows(type="low", threshold = 1/low_threshold, contig_lengths = contig_lengths, allowed_gap = allowed_gap)
    low_merged_regions = filter_nonATGCs_from_low_merged_regions(reference_fasta, low_merged_regions)

    if min_region_length > 0:
        high_merged_regions = high_merged_regions.loc[high_merged_regions["length"] >= min_region_length].reset_index(drop=True)
        low_merged_regions = low_merged_regions.loc[low_merged_regions["length"] >= min_region_length].reset_index(drop=True)

    high_merged_regions = annotate_zero_cov_bases(high_merged_regions, bam_file)
    low_merged_regions = annotate_zero_cov_bases(low_merged_regions, bam_file)

    generate_outputs(reference_fasta, high_merged_regions, low_merged_regions,
                     cov_df, cov_stats, buffer, output_dir, contig_lengths, no_window_stats,
                     log_file)


def preflight_checks(reference_fasta, bam_file, output_dir, exclude_contigs, log_file, full_cmd_executed):

    paths_list = [reference_fasta, bam_file]
    check_files_are_found(paths_list)
    check_bam_file_is_indexed(bam_file)
    check_fasta_file_is_indexed(reference_fasta)
    check_excluded_contigs(reference_fasta, exclude_contigs)
    os.makedirs(f"{output_dir}/mosdepth-files", exist_ok=True)

    log_command_run(full_cmd_executed, output_dir, log_file)


def check_bam_file_is_indexed(bam_file):
    if not Path(bam_file + ".bai").is_file():
        cmd = f"samtools index {bam_file}"
        subprocess.run(cmd, shell=True)


def check_fasta_file_is_indexed(reference_fasta):
    if not Path(reference_fasta + ".fai").is_file():
        cmd = f"samtools faidx {reference_fasta}"
        subprocess.run(cmd, shell=True)


def check_excluded_contigs(reference_fasta, exclude_contigs):
    if exclude_contigs:
        with pysam.FastaFile(reference_fasta) as fasta:
            ref_contigs = set(fasta.references)
        missing_contigs = [contig for contig in exclude_contigs if contig not in ref_contigs]
        if missing_contigs:
            print(f"\n    The following contigs were specified for exclusion but were not found in the reference fasta:\n")
            for contig in missing_contigs:
                print(f"        {contig}")
            print("\n    Please check the contig names and try again.")
            notify_premature_exit()


def get_contig_lengths(reference_fasta):
    with pysam.FastaFile(reference_fasta) as fasta:
        return dict(zip(fasta.references, fasta.lengths))


def generate_sliding_bed_file(reference_fasta, sliding_window_size, step_size, output_dir):
    bed_outfile = f"{output_dir}/mosdepth-files/sliding-windows.bed"
    with pysam.FastaFile(reference_fasta) as fasta, open(bed_outfile, "w") as bed_out:
        for contig, length in zip(fasta.references, fasta.lengths):
            for start in range(0, length - sliding_window_size + 1, step_size):
                end = start + sliding_window_size
                bed_out.write(f"{contig}\t{start}\t{end}\n")

    return bed_outfile


def run_mosdepth(bam_file, window_bed_path, output_dir):
    mosdepth_out_prefix = f"{output_dir}/mosdepth-files/{bam_file[:-4]}"
    cmd = f"mosdepth --by {window_bed_path} {mosdepth_out_prefix} {bam_file}"
    subprocess.run(cmd, shell=True)
    mosdepth_regions_file = f"{mosdepth_out_prefix}.regions.bed.gz"

    return mosdepth_regions_file


def read_mosdepth_regions_file(mosdepth_regions_file):
    cov_df = pd.read_csv(
        mosdepth_regions_file,
        sep="\t",
        header=None,
        names=["contig", "start", "end", "cov"],
        compression="gzip",
    ).astype({"contig": str, "start": int, "end": int, "cov": float})

    return cov_df


class CoverageStats:
    def __init__(self, df, per_contig = False):
        self.df = df.copy()
        self.per_contig = per_contig

        # global stats
        self.global_mean = self.df["cov"].mean()
        self.global_std  = self.df["cov"].std()
        self.global_median = self.df["cov"].median()
        self.global_min = self.df["cov"].min()
        self.global_max = self.df["cov"].max()
        self.df['percentile'] = self.df['cov'].rank(pct=True) * 100

        # per-contig stats
        group = self.df.groupby("contig")["cov"]
        self.df["contig_mean"] = group.transform("mean")
        self.df["contig_std"] = group.transform("std")
        self.df["contig_median"] = group.transform("median")
        self.df["contig_min"] = group.transform("min")
        self.df["contig_max"] = group.transform("max")
        self.df["contig_percentile"] = group.transform(lambda x: x.rank(pct=True) * 100)

        if per_contig:
            self.df["baseline_mean"] = self.df["contig_mean"]
            self.df["baseline_std"] = self.df["contig_std"]
            self.df["percentile"] = self.df["contig_percentile"]
        else:
            self.df["baseline_mean"] = self.global_mean
            self.df["baseline_std"] = self.global_std
            self.df["percentile"] = self.df["percentile"]

        # fold diff, zscore, log2-fold diff
        nobody_likes_zero = 1e-6
        self.df["fold_diff"] = self.df["cov"] / self.df["baseline_mean"]
        self.df["signed_fold_diff"] = self.df["fold_diff"].where(
            self.df["fold_diff"] >= 1,
            -1.0 / self.df["fold_diff"]
        )
        self.df["zscore"] = (self.df["cov"] - self.df["baseline_mean"]) / self.df["baseline_std"]
        self.df["log2_fold_diff"] = np.log2((self.df["cov"] + nobody_likes_zero) / (self.df["baseline_mean"] + nobody_likes_zero))


    def global_percentile(self, p):
        return np.percentile(self.df["cov"], p)


    def get_windows(self, type = "high", method = "fold", threshold = 10):
        """
        Filter for windows above or below a given threshold.

        methods: "percentile", "zscore", "fold", or "log2fold"
        """
        if method == "percentile":
            metric = self.df["percentile"]
        elif method == "zscore":
            metric = self.df["zscore"]
        elif method == "fold":
            metric = self.df["fold_diff"]
        elif method == "log2fold":
            metric = self.df["log2_fold_diff"]
        else:
            raise ValueError(f"Nope nope nope: {method!r}")

        # filter and return a new dataframe
        if type == "low":
            mask = (metric <= threshold) | (self.df["cov"] == 0)
        elif type == "high":
            mask = (metric >= threshold)
        else:
            raise ValueError(f"Nope nope nope: {type!r}")

        return self.df.loc[mask].reset_index(drop=True)


    def merge_windows(self, *, type = "high", method = "fold", threshold = 10, edge_pad = 500, contig_lengths = None, allowed_gap = 0):

        df = self.get_windows(type=type, method=method, threshold=threshold)

        # here, if type == "low", we are filtering out windows that are within <edge_pad> bases of the start/end of a contig
        if type =="low":
            df["contig_length"] = df["contig"].map(contig_lengths)
            df = df[(df.start >= edge_pad) & (df.end <= df.contig_length - edge_pad)].drop(columns="contig_length")

        regions = []

        for contig, group in df.sort_values(["contig", "start"]).groupby("contig", sort=False):
            curr_start = curr_end = None
            covs = []
            # iterating through each window
            for row in group.itertuples(index = False):
                start, end, cov = row.start, row.end, row.cov
                if curr_start is None:
                    # starting a new region
                    curr_start, curr_end = start, end
                    covs = [cov]
                elif start <= curr_end + allowed_gap:
                    # adjacent, or within allowed_gap range, so extending region
                    curr_end = max(curr_end, end)
                    covs.append(cov)
                else:
                    # not adjacent, so saving the region and starting a new one
                    regions.append(self._summarize_region(contig, curr_start, curr_end, covs))
                    curr_start, curr_end = start, end
                    covs = [cov]

            if curr_start is not None:
                # saving the last region
                regions.append(self._summarize_region(contig, curr_start, curr_end, covs))

        return pd.DataFrame(regions,
                            columns=["contig","start","end", "length",
                                     "cov", "percentile", "zscore", "signed_fold_diff",
                                     "fold_diff", "log2_fold_diff"])


    def _summarize_region(self, contig, start, end, covs):
        "helper to compute the summary metrics given a list of window coverages"
        region_length = end - start
        region_mean_cov = float(np.mean(covs))

        if self.per_contig:
            subset = self.df.loc[self.df["contig"] == contig, "cov"]
            percentile = float((subset <= region_mean_cov).mean() * 100)
            baseline_mean = float(self.df.loc[self.df["contig"] == contig, "contig_mean"].iloc[0])
            baseline_std = float(self.df.loc[self.df["contig"] == contig, "contig_std"].iloc[0])
        else:
            percentile = float((self.df["cov"] <= region_mean_cov).mean() * 100)
            baseline_mean = self.global_mean
            baseline_std = self.global_std

        # computing metrics based on appropriate baseline
        fold_diff = region_mean_cov / baseline_mean
        if fold_diff == 0:
            signed_fold_diff = 0.0
        else:
            signed_fold_diff = fold_diff if fold_diff >= 1 else -1.0 / fold_diff
        zscore = (region_mean_cov - baseline_mean) / baseline_std
        nobody_likes_zero = 1e-6
        log2_fold_diff = float(np.log2((region_mean_cov + nobody_likes_zero) / (baseline_mean + nobody_likes_zero)))

        return (contig, start, end, region_length, region_mean_cov, percentile, zscore, signed_fold_diff, fold_diff, log2_fold_diff)


def filter_nonATGCs_from_low_merged_regions(reference_fasta, low_merged_regions):
    # filtering out regions comprised of >= 5% non-ATGCs
        # originally was doing just Ns, but some seqs had odd characters
            # e.g., NC_003070.9:15345880-15346580 in GCF_000001735.4
    max_nonATGCs_fraction = 0.05
    regions_to_keep = []
    with pysam.FastaFile(reference_fasta) as fasta:
        for idx, row in low_merged_regions.iterrows():
            seq = fasta.fetch(row.contig, row.start, row.end).upper()
            total_len = len(seq)
            non_atgc = total_len - (
                seq.count("A") +
                seq.count("T") +
                seq.count("C") +
                seq.count("G")
            )
            non_atgc_frac = non_atgc / total_len

            if non_atgc_frac < max_nonATGCs_fraction:
                regions_to_keep.append(idx)

    low_merged_regions = low_merged_regions.loc[regions_to_keep].reset_index(drop=True)

    return low_merged_regions


def annotate_zero_cov_bases(merged_regions, bam_file):

    zero_counts = []

    with pysam.AlignmentFile(bam_file, "rb") as aln:
        for _, row in merged_regions.iterrows():

            # count_coverage returns 4 arrays (A,C,G,T) of per-base counts
            cov_arrays = aln.count_coverage(
                row.contig,
                start = int(row.start),
                stop = int(row.end),
                quality_threshold = 0
            )

            # sum across nucleotides → total depth at each position
            depth = np.zeros(len(cov_arrays[0]), dtype=int)

            for arr in cov_arrays:
                depth += np.array(arr, dtype=int)
            zero_counts.append(int((depth == 0).sum()))

    merged_regions = merged_regions.copy()
    merged_regions['zero_cov_bases'] = zero_counts

    # reordering columns
    cols = merged_regions.columns.tolist()
    cols.insert(4, cols.pop(cols.index('zero_cov_bases')))

    return merged_regions[cols]


def generate_outputs(reference_fasta, high_merged_regions, low_merged_regions,
                     cov_df, cov_stats, buffer, output_dir, contig_lengths, no_window_stats,
                     log_file):

    print()
    write_window_cov_stats(cov_stats, output_dir, contig_lengths, log_file)
    if not no_window_stats:
        write_windows_table(cov_stats, output_dir, log_file)
    write_window_plot_cov_histogram(cov_df, output_dir, log_file)

    write_regions_of_interest_table(high_merged_regions, output_dir, "high", log_file)
    write_regions_fasta(reference_fasta, high_merged_regions, buffer, output_dir, "high", contig_lengths)

    write_regions_of_interest_table(low_merged_regions, output_dir, "low", log_file)
    write_regions_fasta(reference_fasta, low_merged_regions, buffer, output_dir, "low", contig_lengths)
    print()


def write_window_cov_stats(cov_stats, output_dir, contig_lengths, log_file):
    out_txt = f"{output_dir}/window-coverage-overview.txt"
    out_tsv = f"{output_dir}/window-coverage-overview.tsv"
    percentiles = [0.01, 0.1, 1, 5, 10, 25, 50, 75, 90, 95, 99, 99.9, 99.99]

    rows = []

    # assembly‐level
    assembly_length = sum(contig_lengths.values())
    rows.append({
        "name":        "full-assembly",
        "length_bp":   assembly_length,
        "num_windows": len(cov_stats.df),
        "mean":        cov_stats.global_mean,
        "median":      cov_stats.global_median,
        "std":         cov_stats.global_std,
        "min":         cov_stats.global_min,
        "max":         cov_stats.global_max,
        **{f"p{p}": cov_stats.df["cov"].quantile(p/100) for p in percentiles}
    })

    # per-contig
    for contig, group in cov_stats.df.groupby("contig", sort=False):
        contig_len = contig_lengths.get(contig, 0)
        rows.append({
            "name":        contig,
            "length_bp":   contig_len,
            "num_windows": len(group),
            "mean":        group["cov"].mean(),
            "median":      group["cov"].median(),
            "std":         group["cov"].std(),
            "min":         group["cov"].min(),
            "max":         group["cov"].max(),
            **{f"p{p}": group["cov"].quantile(p/100) for p in percentiles}
        })

    # writing assembly level
    with open(out_txt, "w") as out:
        out.write("====================================\n")
        out.write("=== GLOBAL WINDOW COVERAGE STATS ===\n")
        out.write("====================================\n")
        out.write(f"Assembly size: {assembly_length:,} bp\n")
        out.write(f"Number of windows: {len(cov_stats.df)}\n")
        out.write(f"Mean:   {cov_stats.global_mean:.2f}\n")
        out.write(f"Median: {cov_stats.global_median:.2f}\n")
        out.write(f"Std:    {cov_stats.global_std:.2f}\n")
        out.write(f"Min:    {cov_stats.global_min:.2f}\n")
        out.write(f"Max:    {cov_stats.global_max:.2f}\n\n")

        out.write("Percentiles (global):\n")
        for p in percentiles:
            val = cov_stats.df["cov"].quantile(p/100)
            out.write(f"  {p:>6.2f}% : {val:.2f}\n")

        # now writing per-contig
        out.write("\n========================================\n")
        out.write("=== PER-CONTIG WINDOW COVERAGE STATS ===\n")
        out.write("========================================\n")
        for contig, group in cov_stats.df.groupby("contig", sort=False):
            out.write(f"\n-- {contig} --\n")
            contig_length = contig_lengths.get(contig, 0)
            gm = group["cov"].mean()
            md = group["cov"].median()
            sd = group["cov"].std()
            mn = group["cov"].min()
            mx = group["cov"].max()
            out.write(f"Contig size: {contig_length:,} bp\n")
            out.write(f"Number of windows: {len(group)}\n")
            out.write(f"Mean:   {gm:.2f}\n")
            out.write(f"Median: {md:.2f}\n")
            out.write(f"Std:    {sd:.2f}\n")
            out.write(f"Min:    {mn:.2f}\n")
            out.write(f"Max:    {mx:.2f}\n")
            out.write("\nPercentiles:\n")
            for p in percentiles:
                val = group["cov"].quantile(p/100)
                out.write(f"  {p:>6.2f}% : {val:.2f}\n")
            out.write("\n")

    tee(f"  Window-coverage summary written to:\n      {Fore.YELLOW}{out_txt}", log_file)

    # writing out as tsv
    df_summary = pd.DataFrame(rows)
    cols = ["name","length_bp","num_windows","mean","median","std","min","max"] \
           + [f"p{p}" for p in percentiles]
    df_summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.2f", columns=cols)

    tee(f"  Window-coverage summary table written to:\n      {Fore.YELLOW}{out_tsv}", log_file)


def write_windows_table(cov_stats, output_dir, log_file):
    out_path = f"{output_dir}/window-coverage-stats.tsv.gz"
    cov_stats.df.to_csv(out_path, sep="\t", index=False, float_format="%.2f", compression="gzip",
                        columns=["contig", "start", "end",
                                 "cov", "percentile", "zscore", "signed_fold_diff",
                                 "fold_diff", "log2_fold_diff"])
    tee(f"  Window-coverage stats table written to:\n      {Fore.YELLOW + out_path}", log_file)


def write_window_plot_cov_histogram(cov_df, output_dir, log_file):
    plt.figure()
    plt.hist(np.log10(cov_df["cov"] + 1), bins=100, density=True)
    plt.xlabel("log₁₀(Coverage + 1)")
    plt.ylabel("Frequency")
    plt.title("Distribution of Window Coverages (log₁₀)")
    mean_cov = cov_df["cov"].mean()
    median_cov = cov_df["cov"].median()
    num_windows = len(cov_df)
    ax = plt.gca()
    ax.text(0.95, 0.95,
            f"N: {num_windows:,}\nMean: {mean_cov:.2f}\nMedian: {median_cov:.2f}",
            transform=ax.transAxes,
            fontsize=10,
            va="top",
            ha="right",
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))
    out_path = f"{output_dir}/window-coverage-histogram.png"
    plt.savefig(out_path, dpi=125, bbox_inches="tight")
    plt.close()
    tee(f"  Window-coverage histogram written to:\n      {Fore.YELLOW + out_path}", log_file)


def write_regions_of_interest_table(merged_regions, output_dir, type, log_file):
    out_path = f"{output_dir}/{type}-coverage-regions-of-interest.tsv"

    sorted_df = merged_regions.sort_values(by=["contig", "start"], ascending=[True, True])

    sorted_df.to_csv(out_path, sep="\t", index=False, float_format="%.2f")

    if len(sorted_df) > 0:
        tee(f"\n  Number of {type}-coverage regions identified: {Fore.YELLOW + str(len(sorted_df))}", log_file)
        tee(f"    {type.capitalize()}-coverage regions-of-interest table written to:\n      {Fore.YELLOW + out_path}", log_file)
    else:
        tee(f"\n  No {type}-coverage regions-of-interest identified.", log_file)


def write_regions_fasta(reference_fasta, regions_df, buffer, output_dir, type, contig_lengths):
    output_fasta = f"{output_dir}/{type}-coverage-regions-of-interest.fasta"
    fasta = pysam.FastaFile(reference_fasta)
    records = []

    for _, row in regions_df.iterrows():
        contig = row["contig"]
        orig_start = row["start"]
        orig_end = row["end"]
        cov = row["cov"]

        # adding the buffer to each side of the region
        start_buf = max(0, orig_start - buffer)
        end_buf = min(contig_lengths[contig], orig_end + buffer)

        # getting the seq
        seq = fasta.fetch(contig, start_buf, end_buf)

        # having header include original and buffered coordinates
        header = (f"{contig}:{start_buf}-{end_buf}"
            f"_orig_{orig_start}-{orig_end}"
            f"_cov_{cov:.2f}")

        rec = SeqRecord(Seq(seq),
                        id=header,
                        description="")
        records.append(rec)
    fasta.close()

    if len(records) > 0:
        with open(output_fasta, "w") as out:
            SeqIO.write(records, out, "fasta")

        print(f"    {type.capitalize()}-coverage regions-of-interest fasta written to:\n      {Fore.YELLOW + output_fasta}")
