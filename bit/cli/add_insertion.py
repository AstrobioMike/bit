import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_seed
from bit.modules.general import check_files_are_found, report_message, notify_premature_exit
from bit.modules.add_insertion import add_insertion


def build_parser():

    desc = """
        This script is for adding an insertion sequence to an input fasta. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-add-insertion -i input.fasta -I insertion-sequence.fasta -o output.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group('Required Parameters')
    optional = parser.add_argument_group('Optional Parameters')

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar="<FILE>",
        help="Input fasta file",
        required=True,
    )
    required.add_argument(
        "-I",
        "--insertion-fasta",
        metavar="<FILE>",
        help="Fasta file containing sequence to be inserted (currently only one sequence is supported)",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-fasta",
        metavar="<FILE>",
        help="Name of output fasta file (default: 'output.fasta')",
        default="output.fasta",
    )
    optional.add_argument(
        "-p",
        "--position",
        metavar="<CONTIG:POSITION>",
        help="Position to insert the sequence, e.g., 'contig_1:1000'. Otherwise, or after the first if more than one is requested, it will be random.",
        type=str,
    )
    optional.add_argument(
        "-n",
        "--num-insertions",
        metavar="<INT>",
        help="Number of insertions to add (default: 1).",
        default=1,
        type=int,
    )
    optional.add_argument(
        "--back-to-back",
        action="store_true",
        help="If adding multiple, stack them back-to-back instead of randomly distributing them (default: False)",
    )
    optional.add_argument(
        "-l",
        "--log",
        metavar="<FILE>",
        help="Path to a tsv log file if wanting to record what was done (default: None)",
        default=None,
    )

    add_seed(optional)

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    preflight_checks(args)

    result = add_insertion(args)

    report_results(result, output_fasta=args.output_fasta,
                   back_to_back=args.back_to_back, log_path=args.log)


def preflight_checks(args):

    check_files_are_found([args.input_fasta, args.insertion_fasta])

    if args.num_insertions < 1:
        report_message("`-n|--num-insertions` must be >= 1")
        notify_premature_exit()

    if args.back_to_back and args.num_insertions <= 1:
        report_message("`--back-to-back` only makes sense with `--num-insertions` > 1.")
        notify_premature_exit()


def report_results(result, output_fasta=None, back_to_back=False, log_path=None):

    print("")

    pattern = "back-to-back" if back_to_back else "randomly"

    ins_str = (f"    {result['num_insertions']} insertion(s) of length "
               f"{result['insertion_length']} bps added {pattern} "
               f"and written to:\n\n        {output_fasta}")

    print(ins_str)

    print("")

    if log_path is not None:
        with open(log_path, "w") as log:
            log.write("insertion_number\tcontig\tposition\n")
            for i, (contig, pos) in enumerate(result["positions"], 1):
                log.write(f"{i}\t{contig}\t{pos}\n")
