"""
Expects rows to be units (e.g. genes/KOs/etc.), and columns to be samples.

This script normalizes a table by either coverage per million (CPM) or based on the median-ratio
method as performed in DESeq2. But unlike DESeq2, we don't care here if there are floats in there.

I initially wrote this for normalizing metagenomic coverage data, like gene-level coverage, or summed KO coverages.
These are normalized for gene-length already because they are "coverages", but they are not yet normalized
for sampling depth - which is where this script comes in.

I also found myself wanting this because I wanted to do differential abundance testing of coverages
of KO terms. DESeq2 doesn't require normalizing for gene-length because it is the same unit being analyzed
across all samples - the same gene, so the same size. However, after grouping genes into their KO annotations,
(which we may need to compare across samples that don't all share the same underlying assembly or genes),
they no longer all represent the same units across all samples. It is because of this I decided to stick with
gene-level coverages (which are normalized for gene-length), and then sum those values based on KO annotations.

The CPM (coverage per million) normalization is just like a percent, except scaled to 1 million instead of 100.
So each row's entry (e.g. gene/KO/etc.) is the proportion out of 1 million for that column (sample),
and each column will sum to 1 million.

The median-ration normalization method (MR) was initially described in this paper
(http://dx.doi.org/10.1186/gb-2010-11-10-r106; e.q. 5), and this site is super-informative in general
about the DESeq2 process overall, and helped me understand the normalizaiton process better to implement it:
https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html. Columns will not sum to
the same amount when the median-ratio method is applied.
"""

import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help



def build_parser():

    desc = """
        This script normalizes a table by either counts- or coverage-per-million (CPM) or with the median-ratio
        method as performed in DESeq2. See note at top of module for more info. It expects a
        tab-delimited table with samples as columns and units (e.g. genes/KOs/OTUs/etc.) as rows.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-normalize-table -i input-table.tsv -n CPM -o output-table.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-table",
        metavar="<FILE>",
        help="Input tab-delimited table",
        required=True,
    )

    optional.add_argument(
        "-n",
        "--normalization",
        help='Desired normalization method of either "CPM" for counts- or coverage-per-million or "MR" for median-ratio (default: "CPM")',
        choices=["CPM", "MR"],
        default="CPM",
    )

    optional.add_argument(
        "-o",
        "--output-table",
        metavar="<FILE>",
        help='Output filename (default: "normalized.tsv")',
        default="normalized.tsv",
    )

    add_help(optional)

    return parser


def main(args=None):

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args(args)

    from bit.modules.general import check_files_are_found

    check_files_are_found([args.input_table])

    import pandas as pd # type: ignore

    tab = pd.read_csv(args.input_table, sep="\t", index_col=0, low_memory=False)

    # removing columns with all zeroes prior to normalization (will be restored after)
    tab, zero_column_names, ordered_columns = remove_zero_columns(tab)

    if args.normalization == "CPM":
        norm_tab = normalize_cpm(tab)
    else:
        norm_tab = normalize_median_ratio(tab)

    # restoring zero columns and original column order
    norm_tab = restore_zero_columns(norm_tab, zero_column_names, ordered_columns)

    norm_tab.to_csv(args.output_table, sep="\t")


def remove_zero_columns(tab):
    ordered_columns = tab.columns.tolist()
    column_sums = tab.sum()
    zero_column_names = column_sums[column_sums == 0].index.tolist()
    tab = tab.drop(zero_column_names, axis=1)
    return tab, zero_column_names, ordered_columns


def restore_zero_columns(tab, zero_column_names, ordered_columns):
    for col in zero_column_names:
        tab[col] = 0.0
    return tab[ordered_columns]


def normalize_cpm(tab):
    return tab / tab.sum() * 1000000


def normalize_median_ratio(tab):

    import numpy as np # type: ignore
    from scipy.stats.mstats import gmean # type: ignore

    # getting geometric means for each row
    with np.errstate(divide='ignore'):
        geomeans = gmean(tab, axis=1)

    # getting ratios of values to geometric means
    ratios_tab = (tab.T / geomeans).T

    # calculating size factors from rows with non-zero geometric means
    size_factors = ratios_tab[geomeans > 0].median().to_list()

    return tab / size_factors

if __name__ == "__main__":
    main()
