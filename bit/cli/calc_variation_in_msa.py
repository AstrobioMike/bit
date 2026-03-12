#!/usr/bin/env python
import sys
import argparse
from bit.modules.seqs import calc_variation_in_msa
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)

def main():

    desc = """
        This script takes an alignment in fasta format as input and returns the Shannon uncertainty values for each column
        (using: https://scikit.bio/docs/dev/generated/skbio.alignment.TabularMSA.conservation.html). In output, a "variation" value of 0 would
        mean the same character in all sequences for that position (highest conservation); 1 would mean equal probability of any character
        (greatest variability). "Conservation" column is inverse. As written, any ambiguous bases or residues are converted to gap characters.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-calc-variation-in-msa -i alignment.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-alignment-fasta",
        metavar="<FILE>",
        help="Input alignment fasta file",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tab-separated file (default: "variation.tsv")',
        action="store",
        default="variation.tsv"
    )

    optional.add_argument(
        "-t",
        "--type",
        help='Molecule type (default: "Protein")',
        choices=["Protein", "DNA", "3Di"],
        action="store",
        default="Protein"
    )

    optional.add_argument(
        "-g",
        "--gap-treatment",
        help='How to treat gaps (default: "include")',
        choices=["include", "nan", "ignore", "error"],
        action="store",
        default="include"
    )

    add_help(optional)

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    df = calc_variation_in_msa(args)
    df.to_csv(args.output_tsv, sep="\t", index=False)
