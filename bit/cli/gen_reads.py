import sys
import argparse
from bit.modules.gen_reads import generate_reads
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)


def main():

    desc = """
        This script generates perfect (no error model) paired-end reads in FASTQ format from one
        or multiple input FASTA files. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-gen-reads -i genome-1.fasta genome-2.fasta -p proportions.tsv -o perfect-reads",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )
    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fastas",
        metavar="<FILES>",
        nargs='+',
        help="List of input FASTA files",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help="Prefix for output FASTQ files (default: perfect-reads)",
        default="perfect-reads",
    )
    optional.add_argument(
        "-n",
        "--num-read-pairs",
        metavar="<INT>",
        help="Number of total read-pairs to generate (default: 1,000,000)",
        type=int,
        default=1000000,
    )
    optional.add_argument(
        "-r",
        "--read-length",
        metavar="<INT>",
        help="Length of each read (default: 150)",
        type=int,
        default=150,
    )
    optional.add_argument(
        "-f",
        "--fragment-size",
        metavar="<INT>",
        help="Size of the fragment from which paired reads are generated (default: 500)",
        type=int,
        default=500,
    )
    optional.add_argument(
        "-p",
        "--proportions-file",
        metavar="<FILE>",
        help="Proportions file (tab-delimited) specifying read proportions per input fasta \
              (1st column should be path to input fasta; 2nd column should be \
              wanted proportion with the total summing to one). If not provided, \
              equal proportions are assumed.",
    )
    optional.add_argument(
        "-s",
        "--seed",
        metavar="<INT>",
        help="Set the random seed if wanting control over the random number generator (default: None)",
        type=int,
        default=None,
    )
    add_help(optional)


    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    generate_reads(args)
