import sys
import os
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.get_workflow import dl_wf


def build_parser():

    desc = """
        This is a helper program for downloading bit workflows.
        Workflow version is included with the downloaded workflow.
        For bit version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-workflow metagenomics`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "workflow",
        choices=["metagenomics", "genome-summarize", "sra-download"],
        help="The first positional argument should be which one of these workflows you'd like to download",
    )

    optional.add_argument(
        "--list-available-versions",
        help="Provide this flag along with a specified workflow in order to get a printout of available versions",
        action="store_true"
    )

    optional.add_argument(
        "--wanted-version",
        metavar="VERSION",
        help="Specify the version you'd like to download (leaving out this argument will pull the latest by default)"
    )

    add_help(optional)

    return parser

def main():

    parser = build_parser()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_wf(args)
