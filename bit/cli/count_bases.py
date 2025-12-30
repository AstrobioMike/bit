import sys
import argparse
from bit.modules.general import check_files_are_found
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.seqs import parse_fasta_lengths


def build_parser():

    desc = """
        This script takes a fasta as input and returns the total number of characters if the input
        holds a single sequence, otherwise it produces a tab-delimited file with two columns (header
        and number of characters for each sequence) and prints out some summary stats.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-count-bases -i input.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar="<FILE>",
        help="Input fasta file",
        action="store",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tab-delimited file (default: "num-chars.tsv")',
        action="store",
        default="num-chars.tsv",
    )

    add_help(optional)

    return parser


def main():
    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_files_are_found([args.input_fasta])

    stats = parse_fasta_lengths(args.input_fasta)

    # if only one sequence, just printing out full length
    if stats["stats"]["n_seqs"] == 1:
        length = next(iter(stats["lengths"].values()))
        print(str(length))
    else:
        write_lengths_to_tsv(stats["lengths"], args.output_tsv)
        print_summary(stats["stats"], args.output_tsv)


def write_lengths_to_tsv(seq_lengths, output_tsv):
    """Write sequence lengths to a TSV file."""
    with open(output_tsv, "w") as out_file:
        for seq_id, length in seq_lengths.items():
            out_file.write(f"{seq_id}\t{length}\n")


def print_summary(stats, output_tsv):
    """Print summary statistics."""
    print("\n    Number of seqs:  " + str(stats["n_seqs"]))
    print("    Min. length:     " + str(stats["min"]))
    print("    Max length:      " + str(stats["max"]))
    print("    Mean length:     " + str(stats["mean"]))
    print("    Median length:   " + str(stats["median"]) + "\n")
    print(f"  All seq lengths written to: '{output_tsv}'\n")
