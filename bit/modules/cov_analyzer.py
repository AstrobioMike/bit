import io
import os
import sys
import subprocess
from pathlib import Path
import pysam # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from Bio import SeqIO # type: ignore
from Bio.Seq import Seq # type: ignore
from Bio.SeqRecord import SeqRecord # type: ignore
import matplotlib.pyplot as plt # type: ignore
from tqdm import tqdm # type: ignore
from colorama import Fore, init # type: ignore
from bit.modules.general import (check_if_output_dir_exists,
                                 check_files_are_found,
                                 notify_premature_exit,
                                 report_message,
                                 log_command_run,
                                 tee, spinner)

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
    window_size: int,
    step_size: int,
    allowed_gap: int,
    buffer: int,
    write_window_stats: bool,
    force_overwrite: bool,
    log_file = str,
    full_cmd_executed: str = None
):

    preflight_checks(reference_fasta, bam_file, output_dir, exclude_contigs,
                     force_overwrite, log_file, full_cmd_executed)

    contig_lengths = get_contig_lengths(reference_fasta)

    report_message(f"Generating windows for the input fasta (window size: {window_size}; step size: {step_size})...")
    window_bed_path, total_windows = generate_windows_bed_file(reference_fasta,
                                                               window_size,
                                                               step_size,
                                                               output_dir)

    report_message("Running mosdepth to generate window-coverage data...")
    mosdepth_regions_file = run_mosdepth(bam_file, window_bed_path, output_dir)

    report_message("Reading in coverage data...")
    cov_df = read_mosdepth_regions_file(mosdepth_regions_file, total_windows)
    if exclude_contigs:
        cov_df = cov_df.loc[~cov_df["contig"].isin(exclude_contigs)].reset_index(drop=True)

    report_message("Computing per-window coverage stats...")
    with spinner("", ""):
        cov_stats = CoverageStats(cov_df, per_contig)

    report_message("Identifying regions of interest...")
    with spinner("", ""):
        high_merged_regions = cov_stats.merge_windows(type="high", threshold = high_threshold, contig_lengths = contig_lengths, allowed_gap = allowed_gap)
        low_merged_regions = cov_stats.merge_windows(type="low", threshold = 1/low_threshold, contig_lengths = contig_lengths, allowed_gap = allowed_gap)
        low_merged_regions = filter_nonATGCs_from_low_merged_regions(reference_fasta, low_merged_regions)

        zero_merged_regions = find_zero_coverage_regions(cov_df)
        zero_merged_regions = filter_nonATGCs_from_low_merged_regions(reference_fasta, zero_merged_regions)

        if min_region_length > 0:
            high_merged_regions = high_merged_regions.loc[high_merged_regions["length"] >= min_region_length].reset_index(drop=True)
            low_merged_regions = low_merged_regions.loc[low_merged_regions["length"] >= min_region_length].reset_index(drop=True)
            zero_merged_regions = zero_merged_regions.loc[zero_merged_regions["length"] >= min_region_length].reset_index(drop=True)

        high_merged_regions = annotate_zero_cov_bases(high_merged_regions, bam_file)
        low_merged_regions = annotate_zero_cov_bases(low_merged_regions, bam_file)

        high_merged_regions = annotate_low_complexity(reference_fasta, high_merged_regions)
        low_merged_regions = annotate_low_complexity(reference_fasta, low_merged_regions)
        zero_merged_regions = annotate_low_complexity(reference_fasta, zero_merged_regions)

        for df in [high_merged_regions, low_merged_regions, zero_merged_regions]:
            df.insert(1, "contig_length", df["contig"].map(contig_lengths))

        # high_merged_regions["contig_length"] = high_merged_regions["contig"].map(contig_lengths)
        # low_merged_regions["contig_length"] = low_merged_regions["contig"].map(contig_lengths)
        # zero_merged_regions["contig_length"] = zero_merged_regions["contig"].map(contig_lengths)

    report_message("Generating outputs...")
    with spinner("", ""):
        old_stdout = sys.stdout
        sys.stdout = captured = io.StringIO()
        try:
            generate_outputs(reference_fasta, high_merged_regions, low_merged_regions,
                             zero_merged_regions,
                             cov_df, cov_stats, buffer, output_dir, contig_lengths, write_window_stats,
                             log_file)
        finally:
            sys.stdout = old_stdout
    print()
    print(captured.getvalue(), end="")


def preflight_checks(reference_fasta, bam_file, output_dir, exclude_contigs, force_overwrite, log_file, full_cmd_executed):

    check_if_output_dir_exists(output_dir, force_overwrite)
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


def generate_windows_bed_file(reference_fasta, window_size, step_size, output_dir):
    bed_outfile = f"{output_dir}/mosdepth-files/windows.bed"
    flush_size = 100000
    with pysam.FastaFile(reference_fasta) as fasta, open(bed_outfile, "w") as bed_out:
        contigs = [(c, l) for c, l in zip(fasta.references, fasta.lengths) if l >= window_size]
        total_windows = sum((l - window_size) // step_size + 1 for _, l in contigs)
        with tqdm(total=total_windows, ncols=80, unit=" windows", leave=True, bar_format="    {l_bar}{bar}| of {total:,} total{unit} [{elapsed}, {rate_fmt}]") as pbar:
            for contig, length in contigs:
                buffered_lines = []
                for start in range(0, length - window_size + 1, step_size):
                    end = start + window_size
                    buffered_lines.append(f"{contig}\t{start}\t{end}\n")

                    if len(buffered_lines) >= flush_size:
                        bed_out.write("".join(buffered_lines))
                        pbar.update(len(buffered_lines))
                        buffered_lines.clear()

                if buffered_lines:
                    bed_out.write("".join(buffered_lines))
                    pbar.update(len(buffered_lines))

    return bed_outfile, total_windows


def run_mosdepth(bam_file, window_bed_path, output_dir):
    mosdepth_out_prefix = str(Path(output_dir) / "mosdepth-files" / Path(bam_file).stem)
    cmd = f"mosdepth --by {window_bed_path} {mosdepth_out_prefix} {bam_file}"

    with spinner("mosdepthinitely in progress...", "mosdepthinitely done "):
        subprocess.run(cmd, shell=True)

    mosdepth_regions_file = f"{mosdepth_out_prefix}.regions.bed.gz"

    return mosdepth_regions_file


def read_mosdepth_regions_file(mosdepth_regions_file, total_rows=None):

    chunks = []
    reader = pd.read_csv(
        mosdepth_regions_file,
        sep="\t",
        header=None,
        names=["contig", "start", "end", "cov"],
        compression="gzip",
        dtype={"contig": str, "start": np.int64, "end": np.int64, "cov": np.float64},
        chunksize=500_000,
    )
    with tqdm(total=total_rows, ncols=74, unit=" rows", leave=True, bar_format="    {l_bar}{bar}| of {total:,} total{unit} [{elapsed}, {rate_fmt}]") as pbar:
        for chunk in reader:
            chunks.append(chunk)
            pbar.update(len(chunk))

    cov_df = pd.concat(chunks, ignore_index=True)

    return cov_df


class CoverageStats:
    def __init__(self, df, per_contig = False):
        self.df = df
        self.per_contig = per_contig

        cov = self.df["cov"]

        # global stats
        self.global_mean   = float(cov.mean())
        self.global_std    = float(cov.std())
        self.global_median = float(cov.median())
        self.global_min    = float(cov.min())
        self.global_max    = float(cov.max())

        # percentile + per-contig / global baseline for derived columns
        if per_contig:
            group = self.df.groupby("contig")["cov"]
            baseline_mean = group.transform("mean")
            baseline_std  = group.transform("std")
            self.df["percentile"] = group.rank(pct=True) * 100
        else:
            baseline_mean = self.global_mean
            baseline_std  = self.global_std
            self.df["percentile"] = cov.rank(pct=True) * 100

        # derived columns
        nobody_likes_zero = 1e-6
        fold = cov / baseline_mean
        self.df["fold_diff"] = fold
        self.df["signed_fold_diff"] = fold.where(fold >= 1, -1.0 / fold)
        self.df["zscore"] = (cov - baseline_mean) / baseline_std
        self.df["log2_fold_diff"] = np.log2(
            (cov + nobody_likes_zero) / (baseline_mean + nobody_likes_zero)
        )

        # pre-sorted arrays for fast percentile lookup in _summarize_region
        self._global_cov_sorted = np.sort(self.df["cov"].values)
        self._contig_stats = {}
        for contig, grp in self.df.groupby("contig"):
            vals = grp["cov"].values
            self._contig_stats[contig] = {
                "mean": float(np.mean(vals)),
                "std":  float(np.std(vals, ddof=1)),
                "cov_sorted": np.sort(vals),
            }


    def global_percentile(self, p):
        return np.percentile(self.df["cov"], p)


    def get_windows(self, type = "high", method = "fold", threshold = 10):
        """
        filters for windows above or below a given threshold

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

        # inline filtering (same logic as get_windows) selecting only needed columns
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

        if type == "low":
            mask = (metric <= threshold) | (self.df["cov"] == 0)
        elif type == "high":
            mask = (metric >= threshold)
        else:
            raise ValueError(f"Nope nope nope: {type!r}")

        df = self.df.loc[mask, ["contig", "start", "end", "cov"]]

        # here, if type == "low", we are filtering out windows that are within <edge_pad> bases of the start/end of a contig
        if type == "low":
            contig_len = df["contig"].map(contig_lengths)
            df = df[(df["start"] >= edge_pad) & (df["end"] <= contig_len - edge_pad)]

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
        length = end - start
        region_mean_cov = float(np.mean(covs))

        if self.per_contig:
            stats = self._contig_stats[contig]
            sorted_cov = stats["cov_sorted"]
            percentile = float(np.searchsorted(sorted_cov, region_mean_cov, side="right") / len(sorted_cov) * 100)
            baseline_mean = stats["mean"]
            baseline_std = stats["std"]
        else:
            percentile = float(np.searchsorted(self._global_cov_sorted, region_mean_cov, side="right") / len(self._global_cov_sorted) * 100)
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

        return (contig, start, end, length, region_mean_cov, percentile, zscore, signed_fold_diff, fold_diff, log2_fold_diff)


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


def linguistic_complexity(seq, k=3):
    """
    we get the ratio of observed unique k-mers to the max possible unique k-mers in the seq
    """
    seq = seq.upper()
    n = len(seq)
    if n < k:
        return 0.0
    observed = len({seq[i:i + k] for i in range(n - k + 1)})
    max_possible = min(n - k + 1, 4 ** k)
    return observed / max_possible


# def filter_low_complexity_regions(reference_fasta, merged_regions, min_complexity=0.3, k=3):
#     """
#     this is here to help catch simple tandem repeats that pass
#     the non-ATGC filter but are produce un-informative or false-positive
#     zero/low-coverage calls
#     """
#     if merged_regions.empty:
#         return merged_regions

#     regions_to_keep = []
#     with pysam.FastaFile(reference_fasta) as fasta:
#         for idx, row in merged_regions.iterrows():
#             seq = fasta.fetch(row.contig, row.start, row.end)
#             if linguistic_complexity(seq, k) >= min_complexity:
#                 regions_to_keep.append(idx)

#     return merged_regions.loc[regions_to_keep].reset_index(drop=True)


def annotate_low_complexity(reference_fasta, merged_regions, k=3, min_complexity=0.4):
    """
    rather than removing them, this just flags the low-complexity seqs in the output tsvs
    """
    if merged_regions.empty:
        merged_regions = merged_regions.copy()
        merged_regions["low_complexity"] = pd.Series(dtype=bool)
        return merged_regions

    flags = []
    with pysam.FastaFile(reference_fasta) as fasta:
        for _, row in merged_regions.iterrows():
            seq = fasta.fetch(row.contig, row.start, row.end)
            flags.append(linguistic_complexity(seq, k) <= min_complexity)

    merged_regions = merged_regions.copy()
    merged_regions["low_complexity"] = flags
    return merged_regions


def find_zero_coverage_regions(cov_df):

    "merge adjacent windows with zero coverage into contiguous regions"

    zero_df = cov_df.loc[cov_df["cov"] == 0, ["contig", "start", "end"]].sort_values(["contig", "start"])

    if zero_df.empty:
        return pd.DataFrame(columns=["contig", "start", "end", "length", "cov"])

    regions = []
    for contig, group in zero_df.groupby("contig", sort=False):
        curr_start = curr_end = None
        for row in group.itertuples(index=False):
            if curr_start is None:
                curr_start, curr_end = row.start, row.end
            elif row.start <= curr_end:
                curr_end = max(curr_end, row.end)
            else:
                regions.append((contig, curr_start, curr_end, curr_end - curr_start, 0.0))
                curr_start, curr_end = row.start, row.end
        if curr_start is not None:
            regions.append((contig, curr_start, curr_end, curr_end - curr_start, 0.0))

    return pd.DataFrame(regions, columns=["contig", "start", "end", "length", "cov"])


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
                     zero_merged_regions,
                     cov_df, cov_stats, buffer, output_dir, contig_lengths, write_window_stats,
                     log_file):

    write_window_cov_stats(cov_stats, output_dir, contig_lengths, log_file)
    if write_window_stats:
        write_windows_table(cov_stats, output_dir, log_file)
        write_window_plot_cov_histogram(cov_df, output_dir, log_file)

    write_regions_of_interest_table(high_merged_regions, output_dir, "high", log_file)
    write_regions_fasta(reference_fasta, high_merged_regions, buffer, output_dir, "high", contig_lengths)

    write_regions_of_interest_table(low_merged_regions, output_dir, "low", log_file)
    write_regions_fasta(reference_fasta, low_merged_regions, buffer, output_dir, "low", contig_lengths)

    write_regions_of_interest_table(zero_merged_regions, output_dir, "zero", log_file)
    write_regions_fasta(reference_fasta, zero_merged_regions, buffer, output_dir, "zero", contig_lengths)
    print()


def write_window_cov_stats(cov_stats, output_dir, contig_lengths, log_file):
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

    df_summary = pd.DataFrame(rows)
    cols = ["name","length_bp","num_windows","mean","median","std","min","max"] \
           + [f"p{p}" for p in percentiles]
    df_summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.2f", columns=cols)

    tee(f"\n  Window-coverage overview written to:\n      {Fore.YELLOW}{out_tsv}{Fore.RESET}", log_file)


def write_windows_table(cov_stats, output_dir, log_file):
    out_path = f"{output_dir}/window-coverage-stats.tsv.gz"
    cov_stats.df.to_csv(out_path, sep="\t", index=False, float_format="%.2f", compression="gzip",
                        columns=["contig", "start", "end",
                                 "cov", "percentile", "zscore", "signed_fold_diff",
                                 "fold_diff", "log2_fold_diff"])
    tee(f"  Window-coverage stats table written to:\n      {Fore.YELLOW}{out_path}{Fore.RESET}", log_file)


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
    tee(f"  Window-coverage histogram written to:\n      {Fore.YELLOW}{out_path}{Fore.RESET}", log_file)


def write_regions_of_interest_table(merged_regions, output_dir, type, log_file):
    out_path = f"{output_dir}/{type}-coverage-regions-of-interest.tsv"

    sorted_df = merged_regions.sort_values(by=["contig", "start"], ascending=[True, True])

    if type == "zero":
        sorted_df = sorted_df.drop(columns="cov")

    sorted_df.to_csv(out_path, sep="\t", index=False, float_format="%.2f")

    if len(sorted_df) > 0:
        n_low_complexity = int(sorted_df["low_complexity"].sum()) if "low_complexity" in sorted_df.columns else 0
        lc_note = f" ({n_low_complexity} flagged as low-complexity)" if n_low_complexity > 0 else ""
        tee(f"\n\n  Number of {type}-coverage regions identified: {Fore.YELLOW}{len(sorted_df)}{Fore.RESET}{lc_note}", log_file)
        tee(f"    {type.capitalize()}-coverage regions-of-interest table written to:\n      {Fore.YELLOW}{out_path}{Fore.RESET}", log_file)
    else:
        tee(f"\n\n  No {type}-coverage regions-of-interest identified.", log_file)


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
            for record in records:
                out.write(f">{record.id}\n")
                out.write(f"{str(record.seq)}\n")

        print(f"    {type.capitalize()}-coverage regions-of-interest fasta written to:\n      {Fore.YELLOW}{output_fasta}{Fore.RESET}")
