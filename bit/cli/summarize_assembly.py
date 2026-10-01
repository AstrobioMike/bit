import sys
import argparse
from bit.modules.summarize_assembly import summarize_assemblies
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_version_arg)


def build_parser(parent_subparsers=None):

    desc = """
        This program outputs general summary stats for assemblies provided in fasta
        format. If an output file is specified, it writes the results there as a tsv
        rather than just printing to screen. "Ambiguous characters" reports
        total counts of any letter that is not "A", "T", "C", or "G".
    """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "summarize-assembly",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit summarize-assembly assembly.fasta`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )
    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "input_assemblies",
        metavar="<FILE(s)>",
        nargs="+",
        help="Input assembly file(s)"
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv file (if none provided, prints to screen)',
        default=False
    )
    optional.add_argument(
        "-t",
        "--transpose-output-tsv",
        help='Set this flag if we want to have the output table have genomes as rows rather than columns.',
        action="store_true"
    )
    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=10,
        help="Number of assemblies to summarize in parallel when multiple are provided (default: 10)"
    )

    add_help(optional)

    add_version_arg(optional)

    return parser


def main():
    parser = build_parser()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater")

    summarize_assemblies(
        input_assemblies = args.input_assemblies,
        output_tsv = args.output_tsv,
        transpose_output_tsv = args.transpose_output_tsv,
        jobs = args.jobs
    )
