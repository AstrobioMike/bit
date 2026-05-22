import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.ncbi.get_lineage_from_taxids import get_lineage_from_taxids
from bit.modules.general import check_files_are_found


def build_parser():

    desc = """
        This program uses taxonkit to get NCBI lineage info from taxids. It expects a
        single-column file of taxids with no header and returns a table in the same order.
        Thanks go to taxonkit, don't forget to cite it if using:
        https://bioinf.shenwei.me/taxonkit/. For version info, run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-lineage-from-taxids -i taxids.txt -o lineages.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-file",
        metavar="<FILE>",
        help="Single-column input file of NCBI taxids, no header",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output table filename (default: "lineages.tsv")',
        default="lineages.tsv",
    )

    optional.add_argument(
        "-s",
        "--include-strain",
        help="Add strain info as an additional column if available",
        action="store_true",
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_files_are_found([args.input_file])

    get_lineage_from_taxids(args.input_file, args.output_file, args.include_strain)
