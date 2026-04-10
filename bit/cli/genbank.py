import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        This program extracts different types of information and sequences from GenBank files.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "-i",
            "--input-genbank",
            help = "Input GenBank file",
            metavar = "<FILE>",
            required = True
        )

    def add_common_optional_arguments(group):
        group.add_argument(
            "-o",
            "--output-prefix",
            help = 'Output prefix (default: "extracted")',
            metavar = "<STR>",
            default = "extracted",
            type = str
        )

    ### subcommand cli for extracting full fasta ###
    to_fasta_desc = """
        This subcommand extracts the full fasta sequence.
        """

    to_fasta_parser = subparsers.add_parser(
        "to-fasta",
        help="Extract the full nucleotide fasta sequence",
        description=to_fasta_desc,
        epilog="Ex. usage: `bit-genbank to-fasta -i input.gbff`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    to_fasta_required = to_fasta_parser.add_argument_group("Required Parameters")
    to_fasta_optional = to_fasta_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(to_fasta_required)
    add_common_optional_arguments(to_fasta_optional)

    add_help(to_fasta_optional)

    to_fasta_parser.set_defaults(func="to-fasta")

    ### subcommand cli for extracting AA sequences ###
    to_AA_seqs_desc = """
        This subcommand extracts amino-acid sequences for complete coding sequences.
        It filters out gene calls that:
        1) are noted as "frameshifted", "internal stop", "incomplete", or "pseudo";
        2) have a "location" containing "join" (spanning multiple contigs) or "<" or ">" (running off a contig); or
        3) have an overlapping translation frame ("transl_except").
        """

    to_AA_seqs_parser = subparsers.add_parser(
        "to-AA-seqs",
        help="Extract amino-acid sequences for complete coding sequences",
        description=to_AA_seqs_desc,
        epilog="Ex. usage: `bit-genbank to-AA-seqs -i input.gbff`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    to_AA_seqs_required = to_AA_seqs_parser.add_argument_group("Required Parameters")
    to_AA_seqs_optional = to_AA_seqs_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(to_AA_seqs_required)
    add_common_optional_arguments(to_AA_seqs_optional)

    add_help(to_AA_seqs_optional)

    to_AA_seqs_parser.set_defaults(func="to-AA-seqs")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv) == 1: # pragma: no cover
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

    if args.subcommand == "to-fasta":
        from bit.modules.genbank import genbank_to_fasta
        output_fasta = f"{args.output_prefix}.fasta"
        genbank_to_fasta(args.input_genbank, output_fasta)

    elif args.subcommand == "to-AA-seqs":
        from bit.modules.genbank import genbank_to_AA_seqs
        output_fasta = f"{args.output_prefix}.faa"
        genbank_to_AA_seqs(args.input_genbank, output_fasta)
