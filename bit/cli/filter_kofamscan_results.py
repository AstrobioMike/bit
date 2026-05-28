import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from bit.modules.filter_kofamscan_results import filter_kofamscan_results
from bit.modules.general import check_files_are_found


def build_parser(parent_subparsers=None):

    desc = """
        This program filters the "detail-tsv"-formatted output file from KOFamScan to retain
        only those above the KO-specific score threshold, and retains only the hit with the
        lowest e-value for each gene if there are multiple (thereby simplifying things to one
        KO annotation per input gene). It outputs a 3-column tab-delimited file
        with: gene_ID, KO_ID, and KO_annotation.
    """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "filter-ko-results",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit filter-kofamscan-results -i initial-KOFamScan-results.txt -o KO-annotations.tsv`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-file",
        metavar="<FILE>",
        help="Input KOFamScan detail-tsv annotation table",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output table filename (default: "output.tsv")',
        default="output.tsv",
    )

    add_help(optional)

    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_files_are_found([args.input_file])

    filter_kofamscan_results(args.input_file, args.output_file)
