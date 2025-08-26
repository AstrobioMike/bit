#!/usr/bin/env python
import sys
import argparse
from bit.modules.gen_kraken2_tax_plots import gen_kraken2_tax_plots
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)

def main():

    desc = """
        This script takes as input a kraken2 report and generates bar plots for the most abundant taxa, or those
        above a specified % threshold, at each rank. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-gen-kraken2-tax-plots -i kraken2.report`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-kraken2-report",
        metavar="<FILE>",
        help="Input kraken2 report",
        required=True,
        dest="input_kraken"
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help='Prefix for output plots if wanted (default: "")',
        default=""
    )

    optional.add_argument(
        "-m",
        "--maximum-taxa-per-plot",
        metavar="<INT>",
        help='The maximum number of unique taxa included in a plot (the rest are represented in "Other"; default: 10)',
        type=int,
        default=10,
        dest="max_taxa"
    )

    optional.add_argument(
        "-p",
        "--minimum-percent-threshold",
        metavar="<FLOAT>",
        help='The minimum percent of reads a taxon must have in order to be included (the rest are represented in "Other"). \
            This overrides --maximum-taxa-per-plot if set > 0. (default: 0.0 and to use --maximum-taxa-per-plot instead)',
        type=float,
        default=0.0,
        dest="min_percent"
    )
    optional.add_argument(
        "--no-annots",
        metavar="",
        help="Don't add annotations (total num. reads; perc. unclassified; perc. stuck at root) as overlain text on plots",
    )

    add_help(optional)

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    gen_kraken2_tax_plots(args.input_kraken, args.output_prefix, args.max_taxa, args.min_percent, args.no_annots)
