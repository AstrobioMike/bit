import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.seqs import fasta_to_bed


def build_parser():

    desc = """
        This script takes a nucleotide multifasta as input and returns a tab-delimited bed file
        (see: https://bedtools.readthedocs.io/en/latest/content/general-usage.html).
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-fasta-to-bed -i input.fasta -o output.bed`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Name of output bed file (default: "output.bed")',
        default="output.bed"
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    bed_records = fasta_to_bed(args.input_fasta)

    with open(args.output_file, "w") as out:
        for name, start, end in bed_records:
            out.write(f"{name}\t{start}\t{end}\n")
