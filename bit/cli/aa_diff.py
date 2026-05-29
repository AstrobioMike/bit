import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from bit.modules.general import check_files_are_found, report_message, notify_premature_exit

def build_parser(parent_subparsers=None):

    desc = """
        This program compares an input query sequence (AA or nt) to a reference AA sequence
        and reports the differences relative to the reference. It expects them in single fasta format,
        no multifastas currently accepted.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "aa-diff",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit aa-diff -i query.fa -r ref.faa`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group('Required Parameters')
    optional = parser.add_argument_group('Optional Parameters')

    required.add_argument(
        "-i",
        "--input-query-fa",
        metavar="<FILE>",
        help="Input query fasta (AA or nt) to compare to reference",
        required=True,
    )

    required.add_argument(
        "-r",
        "--ref-faa",
        metavar="<FILE>",
        help="Reference AA fasta",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help="Prefix for output files (default: 'aa-diff')",
        default="aa-diff",
    )

    optional.add_argument(
        "-t",
        "--type",
        choices=["auto-detect", "aa", "nt"],
        help="Type of input query sequence (default: 'auto-detect')",
        default="auto-detect",
    )

    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    args = preflight_checks(args)

    from bit.modules.aa_diff import run_aa_diff
    run_aa_diff(args.input_query_fa, args.ref_faa, args.type, args.output_prefix)


def preflight_checks(args):

    check_files_are_found([args.ref_faa, args.input_query_fa])

    for fasta_path in [args.ref_faa, args.input_query_fa]:
        check_fasta_has_single_seq(fasta_path)

    from bit.modules.seqs import identify_seq_type
    if args.type == "auto-detect":
        args.type = identify_seq_type(args.input_query_fa)

    ref_type = identify_seq_type(args.ref_faa)
    if ref_type != "aa":
        report_message(f"Reference fasta '{args.ref_faa}' must be in amino-acid (AA) format.")
        notify_premature_exit()

    return args


def check_fasta_has_single_seq(fasta_path):

    from Bio import SeqIO # type: ignore

    seq_iter = SeqIO.parse(fasta_path, "fasta")
    first_seq = next(seq_iter, None)
    second_seq = next(seq_iter, None)

    if first_seq is None or second_seq is not None:
        report_message(f"Input fasta '{fasta_path}' must contain only a single sequence.")
        notify_premature_exit()

