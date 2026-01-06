import sys
import argparse
from bit.modules.get_test_data import dl_test_data
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        This is a program for downloading test data files.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-get-test-data metagenomics",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "datatype",
        choices=["metagenomics"],
        help="The first positional argument should be what type of test data you'd like to download",
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_test_data(args)


