#!/usr/bin/env python

import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import textwrap
from bit.modules.general import check_files_are_found


RANK_DICT = {
    "D": "Domain",
    "P": "Phylum",
    "C": "Class",
    "O": "Order",
    "F": "Family",
    "G": "Genus",
    "S": "Species",
}


def gen_kraken2_tax_plots(input_kraken, output_prefix, max_taxa, min_percent, no_annots):

    preflight_checks(input_kraken)

    df = read_in_kraken_report(input_kraken)

    total_reads, unclassified_percent, percent_reads_stuck_at_root = get_stats(df)

    df = reduce_kraken_report(df)

    df = add_domain_letter_to_taxon_name(df)

    for code, rank_name in RANK_DICT.items():

        sub_df = df[df["rank_code"] == code].copy()
        sub_df = df[(df["rank_code"] == code) | ((code == "D") & (df["rank_code"] == "R1"))].copy()


        if sub_df.empty:
            continue

        sub_df = sub_df.groupby("taxon_name", as_index = False)["percent"].sum() # in case there are duplicates (which i don't think should happen)

        if min_percent > 0:

            above_thresh_df = sub_df[sub_df["percent"] >= min_percent].copy()
            below_thresh_df = sub_df[sub_df["percent"] < min_percent].copy()

            sum_others = below_thresh_df["percent"].sum()

            if above_thresh_df.empty and sum_others == 0:

                print(f"No taxa pass the {min_percent}% threshold at rank {rank_name}. Skipping plot.")
                continue

            if sum_others > 0:

                other_df = pd.DataFrame({"taxon_name": ["Other"], "percent": [sum_others]})
                above_thresh_df = pd.concat([above_thresh_df, other_df], ignore_index=True)

            make_barplot_threshold(above_thresh_df, rank_name, output_prefix, total_reads,
                                   unclassified_percent, percent_reads_stuck_at_root, no_annots)

        else:

            make_barplot_top_N(sub_df, rank_name, output_prefix, max_taxa, total_reads,
                               unclassified_percent, percent_reads_stuck_at_root, no_annots)


def preflight_checks(input_report):

    check_files_are_found([input_report])


def read_in_kraken_report(filepath: str) -> pd.DataFrame:
    """
    Expected columns per line (in one typical format):

        1. Percentage of reads covered by the clade rooted at this taxon
        2. Number of reads covered by the clade rooted at this taxon
        3. Number of reads assigned directly to this taxon
        4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
        5. NCBI taxonomy ID
        6. indented scientific name
    """

    df = pd.read_csv(filepath, sep = "\t", header = None,
        names=["percent", "clade_reads", "this_level_reads",
               "rank_code", "tax_id", "taxon_name"],
        skip_blank_lines = True,
        dtype={
            "percent": float,
            "clade_reads": int,
            "this_level_reads": int,
            "rank_code": str,
            "tax_id": str,
            "taxon_name": str
        }
    )

    # removing leading/trailing whitespace from the 'taxon_name' column, as it's indented by rank
    df["taxon_name"] = df["taxon_name"].str.strip()

    return df


def get_stats(df):

    total_reads = df.loc[df["rank_code"] == "U", "clade_reads"].sum() + df.loc[df["taxon_name"] == "root", "clade_reads"].sum()
    unclassified_percent = df.loc[df["rank_code"] == "U", "percent"].sum()
    total_reads_stuck_at_root = df.loc[df["taxon_name"] == "root", "this_level_reads"].sum()
    percent_reads_stuck_at_root = total_reads_stuck_at_root / total_reads * 100
    percent_reads_stuck_at_root = f"{percent_reads_stuck_at_root:.2f}".rstrip('0').rstrip('.')

    return total_reads, unclassified_percent, percent_reads_stuck_at_root


def reduce_kraken_report(df):
    """
    Simplify the Kraken report by removing unnecessary columns and rows.

    We'll keep the taxon name, rank code, and number of reads assigned to that taxon.
    """

    target_ranks = set(RANK_DICT.keys())

    # adding unclassied and root, as i might do something with them at some point,
    # and keeping R1 if present, like it GTDB tax
    target_ranks.update(["U", "R", "R1"])
    filtered_df = df.loc[df["rank_code"].isin(target_ranks)].copy()

    return filtered_df[["rank_code", "taxon_name", "percent", "clade_reads"]]


def add_domain_letter_to_taxon_name(df):
    """ Adding a domain letter to each taxon name, for easier interpretation. """

    domain_map = {
        "Bacteria": "B",
        "Archaea": "A",
        "Eukaryota": "E",
        "Viruses": "V",
        "root": "R",
        "unclassified": "U",
    }

    current_domain_letter = "?"

    for i, row in df.iterrows():

        # "R1" is here for GTDB tax reports that use it for domain
        if row["rank_code"] in {"D", "R1"} and row["taxon_name"] in domain_map:

            this_domain_name = row["taxon_name"]
            current_domain_letter = domain_map.get(this_domain_name, "?")

        else:

            df.at[i, "taxon_name"] = f"({current_domain_letter}) {row['taxon_name']}"

    return df


def make_barplot_threshold(sub_df, rank_name, out_prefix,
                           total_reads, unclassified_percent,
                           percent_reads_stuck_at_root, no_annots):
    """
    Make a barplot for all taxa that meet the user-specified min-percent threshold.
    """
    df_sorted = sub_df.sort_values("percent", ascending=False).reset_index(drop=True)

    # No 'Other' slice is automatically created here, because the user wants
    # to strictly see only those that exceed the threshold.

    plot_barplot_and_save(df_sorted, rank_name, out_prefix, total_reads,
                          unclassified_percent, percent_reads_stuck_at_root, no_annots)


def make_barplot_top_N(sub_df, rank_name, out_prefix, max_taxa, total_reads,
                       unclassified_percent, percent_reads_stuck_at_root, no_annots):
    """
    Make a barplot for the top N taxa, combining the remainder into 'Other'.
    """
    df_sorted = sub_df.sort_values("percent", ascending=False).reset_index(drop=True)
    top_df = df_sorted.head(max_taxa).copy()

    sum_others = df_sorted["percent"][max_taxa:].sum()
    if sum_others > 0:
        other_df = pd.DataFrame({"taxon_name": ["Other"], "percent": [sum_others]})
        top_df = pd.concat([top_df, other_df], ignore_index=True)

    plot_barplot_and_save(top_df, rank_name, out_prefix, total_reads,
                          unclassified_percent, percent_reads_stuck_at_root, no_annots)


def plot_barplot_and_save(df_to_plot, rank_name, out_prefix, total_reads,
                          unclassified_percent, percent_reads_stuck_at_root, no_annots):
    """
    Given a (filtered) dataframe, plot and save the barplot.
    """

    df_sorted = df_to_plot.copy()

    # making sure "Other" is the last entry
    if "Other" in df_sorted["taxon_name"].values:

        other_row = df_sorted[df_sorted["taxon_name"] == "Other"]
        df_sorted = df_sorted[df_sorted["taxon_name"] != "Other"]
        df_sorted = pd.concat([df_sorted, other_row], ignore_index=True)

    df_to_plot = df_sorted.copy()

    def wrap_name(x):
        return "\n".join(textwrap.wrap(x, width=30))

    df_to_plot["taxon_name_wrapped"] = df_to_plot["taxon_name"].apply(wrap_name)

    x_positions = range(len(df_to_plot))
    main_title = f"{rank_name.capitalize()}-level dominant taxa"

    fig, ax = plt.subplots(figsize=(6, 4.5))

    ax.bar(x_positions, df_to_plot["percent"], color="steelblue")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(df_to_plot["taxon_name_wrapped"], rotation=45, ha='right')

    ax.set_title(f"{main_title}")
    ax.set_ylabel("% of total reads", fontsize=10)

    if not no_annots:
        additional_text = (f"Sample total reads: {total_reads:,}\n"
                            f"Percent unclassified: {unclassified_percent:.2f}%\n"
                            f"Percent reads solely at Root: {percent_reads_stuck_at_root}%")
        ax.text(0.97, 0.97, additional_text, transform=ax.transAxes,
                ha='right', va='top', fontsize=9, color="black",
                bbox=dict(facecolor='white', alpha=0.6, linewidth=0.2, boxstyle="round,pad=0.2"))


    fig.tight_layout()

    if out_prefix == "":
        out_png = f"{rank_name.lower()}-barplot.png"
    else:
        out_png = f"{out_prefix}-{rank_name.lower()}-barplot.png"

    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"Generated: {out_png}")
