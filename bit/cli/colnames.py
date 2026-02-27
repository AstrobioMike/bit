import argparse
import sys
from bit.modules.general import check_files_are_found, colnames
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        Returns the column names (with numbers) from a delimited file. The
        delimiter is auto-detected (from tab, comma, pipe, semicolon, or space).
        For version info, run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-colnames input.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "input_file",
        metavar="<FILE>",
        help="Input delimited file",
    )

    add_help(optional)

    return parser


def main(args=None):

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args(args)

    check_files_are_found([args.input_file])

    colnames(args)


if __name__ == "__main__":
    main()
