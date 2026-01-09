import sys
import argparse
from Bio import SeqIO
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.seqs import fasta_to_genbank


def build_parser():

    desc = """
        This script takes a nucleotide fasta file and converts it into minimal genbank format.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-fasta-to-genbank -i input.fasta -o output.gb`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        help="Input nucleotide fasta file",
        metavar="<FILE>",
        required=True
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output genbank file (default: "output.gb")',
        default="output.gb"
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    sequences = fasta_to_genbank(args.input_fasta)

    with open(args.output_file, "w") as out:
        SeqIO.write(sequences, out, "genbank")
