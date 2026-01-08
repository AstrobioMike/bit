import os
import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_help)
from bit.modules.lineage_to_tsv import convert_lineage_to_tsv
from bit.modules.general import check_files_are_found

def build_parser():

    desc = """
        This script converts lineages (in, e.g., this format: "d__Bacteria;p__Campylobacterota;c__Campylobacteria",
        only including up to the standard 7 ranks d,p,c,o,f,g,s) into tsv format (e.g., the previous would be output
        tab-separated as "Bacteria\tCampylobacterota\tCampylobacteria\tNA\tNA\tNA\tNA"). It expects as input a 2-column
        tab-delimited file with column 1 holding an identifier and column 2 holding the lineage. For version info,
        run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-lineage-to-tsv -i input-lineages.tsv -o formatted-tax.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-lineage",
        metavar="<FILE>",
        help="Input lineage file",
        required=True
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output standard taxonomy tsv file (default: formatted-tax.tsv)',
        default="formatted-tax.tsv"
    )

    optional.add_argument(
        "--make-taxid",
        help="Provide this flag to make a unique taxid (string of all rank fields) for each lineage (" \
             "will be added as second column of output).",
        action="store_true"
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_files_are_found([args.input_lineage])

    convert_lineage_to_tsv(args.input_lineage, args.output_tsv, args.make_taxid)
