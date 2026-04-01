import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.extract_seqs import extract_seqs_by_coords, extract_seqs_by_primers


def build_parser():

    desc = """
        This program extracts sequences from an input file based on either coordinates provided in a bed file
        or primer sequences provided to the command line. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Use `bit-extract-seqs by-coords -h` or `bit-extract-seqs by-primers -h` to see subcommand-specific help.",
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

    ### subcommand cli for extracting sequences by primers ###
    by_primers_desc = """
        This subcommand takes a fasta file and forward and reverse primer sequences, and it
        returns a multifasta of the sequences including the specified primers.
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

    args.func(args)

if __name__ == "__main__":
    main()
