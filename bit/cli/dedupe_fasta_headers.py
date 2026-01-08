import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.seqs import dedupe_fasta_headers


def build_parser():

    desc = """
        This script will append a number to headers if that exact ID has already
        appeared in the fasta file. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-dedupe-fasta-headers -i input.fasta -o output.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar="<FILE>",
        help="Starting fasta file",
        action="store",
        required=True
    )

    optional.add_argument(
        "-o",
        "--output-fasta",
        metavar="<FILE>",
        help='Output fasta file (default: "output.fasta").',
        default="output.fasta"
    )

    add_help(optional)

    return parser


def main():
    parser = build_parser()
    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dedupe_fasta_headers(args.input_fasta, args.output_fasta)
