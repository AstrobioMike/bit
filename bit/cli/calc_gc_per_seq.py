import sys
import argparse
from bit.modules.seqs import calc_gc_per_seq
from bit.cli.common import CustomRichHelpFormatter

def main():

    desc = """
        This script takes a nucleotide multifasta and returns
        a tab-delimited file with 3 columns: header, sequence length,
        and GC. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-calc-gc-per-seq -i input.fasta -o output.tsv",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")
    required.add_argument(
        "-i",
        "--input-fasta",
        metavar="<FILE>",
        help='Input fasta file',
        required=True,
    )
    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv file (default: "GC-out.tsv")',
        default="GC-out.tsv",
    )
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit",
    )

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    gc_dict = calc_gc_per_seq(input_fasta=args.input_fasta)

    with open(args.output_tsv, "w") as out_tsv:
        out_tsv.write("header\tlength\tgc\n")
        for item in gc_dict:
            out_tsv.write(f"{item['header']}\t{item['length']}\t{item['gc']}\n")
