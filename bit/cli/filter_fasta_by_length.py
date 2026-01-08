import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.seqs import filter_fasta_by_length
from bit.modules.general import check_files_are_found


def build_parser():

    desc = """
        This script takes a multifasta as input and filters out sequences based on length.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-filter-fasta-by-length -i input.fasta -m 1000 -o filtered.fasta",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        help = "Input fasta file",
        metavar = "<FILE>",
        required = True
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar = "<FILE>",
        help = 'Name of output fasta file (default: "filtered.fasta")',
        default = "filtered.fasta"
    )
    optional.add_argument(
        "-m",
        "--min-length",
        metavar = "<INT>",
        help = "Minimum length retained",
        default = None
    )
    optional.add_argument(
        "-M",
        "--max-length",
        metavar = "<INT>",
        help = "Maximum length retained",
        default = "9223372036854775807",
    )
    add_help(optional)

    return parser


def main():
    parser = build_parser()
    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_files_are_found([args.input_fasta])

    in_fasta = args.input_fasta
    out_fasta = args.output_file
    min_length = int(args.min_length) if args.min_length else None
    max_length = int(args.max_length)

    (num_initial_seqs, num_seqs_retained,
     num_initial_bases, num_bases_retained) = filter_fasta_by_length(in_fasta, out_fasta,
                                                                     min_length, max_length)

    perc_seqs_retained = round(float(num_seqs_retained) / float(num_initial_seqs) * 100, 2)
    perc_bases_retained = round(float(num_bases_retained) / float(num_initial_bases) * 100, 2)

    print("\n    Retained " + f"{num_seqs_retained:,}" + " sequence(s) of the initial " + f"{num_initial_seqs:,}" + " (" + str(perc_seqs_retained) + "%).\n")
    print("    Retained " + f"{num_bases_retained:,}" + " bases of the initial " + f"{num_initial_bases:,}" + " (" + str(perc_bases_retained) + "%).\n")
