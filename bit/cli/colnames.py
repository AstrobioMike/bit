import argparse
import sys
from bit.modules.general import check_files_are_found, colnames
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        Returns the column names (numbered) from a delimited file. The
        delimiter is auto-detected for tab, comma, pipe, semicolon, or space.
        For version info, run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-colnames input.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "input_file",
        metavar="<FILE>",
        nargs='?',
        default=sys.stdin,
        help="Input delimited file or stdin if none provided"
    )

    add_help(optional)

    return parser


def main(args=None):

    parser = build_parser()

    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args(args)

    if args.input_file is not sys.stdin:
        check_files_are_found([args.input_file])

    colnames(args)


if __name__ == "__main__":
    main()
