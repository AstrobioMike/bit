import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import report_message, notify_premature_exit
from bit.modules.extract_seqs import (extract_seqs_by_coords,
                                      extract_seqs_by_primers,
                                      extract_seqs_by_headers)


def build_parser():

    desc = """
        This program extracts sequences from an input fasta through various methods, see subcommand-specific
        help menus for more info. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        # epilog="Use `bit-extract-seqs by-coords -h` or `bit-extract-seqs by-primers -h` to see subcommand-specific help.",
        formatter_class=CustomRichHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "-i",
            "--input-fasta",
            help = "Input fasta file",
            metavar = "<FILE>",
            required = True
        )

    def add_common_optional_arguments(group):
        group.add_argument(
            "-o",
            "--output-fasta",
            help = 'Output fasta file (default: "extracted-seqs.fasta")',
            metavar = "<FILE>",
            default = "extracted-seqs.fasta"
        )


    ### subcommand cli for extracting sequences by coordinates ###
    by_coords_desc = """
        This subcommand takes a fasta file and tab-delimited (bed) file specifying which contigs
        and coordinates are wanted, and it returns a multifasta of the chopped out sequences.
        """

    by_coords_parser = subparsers.add_parser(
        "by-coords",
        help="Extract sequences based on coordinates provided in a bed file",
        description=by_coords_desc,
        epilog="Ex. usage: `bit-extract-seqs by-coords -i input.fasta -b targets.bed`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    by_coords_required = by_coords_parser.add_argument_group("Required Parameters")
    by_coords_optional = by_coords_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(by_coords_required)
    by_coords_required.add_argument(
        "-b",
        "--bed-file",
        help="Tab-delimited bed file of desired contigs and coordinates (3 columns - contig, start, end - no header, 0-based counting)",
        metavar="<FILE>"
    )

    add_common_optional_arguments(by_coords_optional)

    add_help(by_coords_optional)

    by_coords_parser.set_defaults(func=extract_seqs_by_coords)


    ### subcommand cli for extracting sequences by header ###
    by_headers_desc = """
        This subcommand takes a fasta file and specified headers and extracts the sequences
        with those headers (or does the inverse if wanted).
        """

    by_headers_parser = subparsers.add_parser(
        "by-headers",
        help="Extract sequences based on specified headers",
        description=by_headers_desc,
        epilog="Ex. usage: `bit-extract-seqs by-headers -i input.fasta -h contig-1 contig-2`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    by_headers_required = by_headers_parser.add_argument_group("Required Parameters (choose one of `--headers` or `--file-with-headers)")
    by_headers_optional = by_headers_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(by_headers_required)

    by_headers_required.add_argument(
        "-H",
        "--headers",
        help="Headers of sequences (space-delimited if more than one)",
        metavar="<STR>",
        nargs="+"
    )

    by_headers_required.add_argument(
        "-f",
        "--file-with-headers",
        help="File with headers of sequences (one header per line)",
        metavar="<FILE>"
    )

    add_common_optional_arguments(by_headers_optional)

    by_headers_optional.add_argument(
        "--inverse",
        help="If specified, we will extract all sequences [bold]other[/bold] than the provided headers (default: False)",
        action="store_true"
    )

    add_help(by_headers_optional)

    by_headers_parser.set_defaults(func=extract_seqs_by_headers)

    ### subcommand cli for extracting sequences by primers ###
    by_primers_desc = """
        This subcommand takes a fasta file and forward and reverse primer sequences, and it
        returns a multifasta of the sequences including the specified primers. It currently doesn't
        allow for degenerate bases in the primers, but it does allow for up to 2 mismatches.
        """

    by_primers_parser = subparsers.add_parser(
        "by-primers",
        help="Extract sequences based on forward and reverse primer sequences",
        description=by_primers_desc,
        epilog="Ex. usage: `bit-extract-seqs by-primers -i input.fasta -f ForwardPrimerSeq -r ReversePrimerSeq`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    by_primers_required = by_primers_parser.add_argument_group("Required Parameters")
    by_primers_optional = by_primers_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(by_primers_required)
    by_primers_required.add_argument(
        "-f",
        "--forward-primer",
        help="Forward primer sequence",
        metavar="<STR>"
    )
    by_primers_required.add_argument(
        "-r",
        "--reverse-primer",
        help="Reverse primer sequence",
        metavar="<STR>"
    )

    add_common_optional_arguments(by_primers_optional)

    by_primers_optional.add_argument(
        "-m",
        "--max-mismatches",
        help="Maximum number of mismatches allowed between the primer and the target sequence (default: 0)",
        # metavar="<INT>",
        type=int,
        choices=[0, 1, 2],
        default=0
    )

    add_help(by_primers_optional)

    by_primers_parser.set_defaults(func=extract_seqs_by_primers)

    return parser

def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so approriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            parser.subparsers.choices[cmd].print_help(sys.stderr)
            sys.exit(0)

        print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if hasattr(args, "func"):
        if args.func == extract_seqs_by_headers:
            check_by_headers_required_inputs(args)

    args.func(args)

if __name__ == "__main__":
    main()


def check_by_headers_required_inputs(args):
    if not args.headers and not args.file_with_headers:
        report_message("You must provide either -H/--headers or -f/--file-with-headers.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()
    if args.headers and args.file_with_headers:
        report_message("You have provided both -H/--headers and -f/--file-with-headers parameters, please only provide one.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()
