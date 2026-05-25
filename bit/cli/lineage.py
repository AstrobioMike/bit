import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from bit.modules.general import check_files_are_found


def build_parser():

    desc = """
        This program has helpers for getting lineages from taxids and working with lineages. See subcommand-specific
        help menus for more info.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(parser)

    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers


    ### from-taxids subcommand ###
    from_taxids_desc = """
        This subcommand uses taxonkit to get NCBI lineage info from taxids. It expects a
        single-column file of taxids with no header and returns a table in the same order.
        Thanks go to taxonkit, don't forget to cite it if using:
        https://bioinf.shenwei.me/taxonkit/.
        """

    from_taxids_parser = subparsers.add_parser(
        "from-taxids",
        help="Get NCBI lineage info from taxids",
        description=from_taxids_desc,
        epilog="Ex. usage: `bit-lineage from-taxids -i taxids.txt -o lineages.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    from_taxids_required = from_taxids_parser.add_argument_group("Required Parameters")
    from_taxids_optional = from_taxids_parser.add_argument_group("Optional Parameters")

    from_taxids_required.add_argument(
        "-i",
        "--input-file",
        metavar="<FILE>",
        help="Single-column input file of NCBI taxids, no header",
        required=True,
    )

    from_taxids_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output table filename (default: "lineages.tsv")',
        default="lineages.tsv",
    )

    from_taxids_optional.add_argument(
        "-s",
        "--include-strain",
        help="Add strain info as an additional column if available",
        action="store_true",
    )

    add_help(from_taxids_optional)

    add_version_arg(from_taxids_optional)

    from_taxids_parser.set_defaults(func="from_taxids")


    ### to-tsv subcommand ###
    to_tsv_desc = """
        This subcommand converts condensed lineages (in, e.g., this format: "d__Bacteria;p__Campylobacterota;c__Campylobacteria",
        only including up to the standard 7 ranks d,p,c,o,f,g,s) into tsv format (e.g., the previous would be output
        tab-separated as "Bacteria\tCampylobacterota\tCampylobacteria\tNA\tNA\tNA\tNA"). It expects as input a 2-column
        tab-delimited file with column 1 holding an identifier and column 2 holding the lineage.
        """

    to_tsv_parser = subparsers.add_parser(
        "to-tsv",
        help="Convert condensed lineages to TSV format",
        description=to_tsv_desc,
        epilog="Ex. usage: `bit-lineage to-tsv -i input-lineages.tsv -o formatted-tax.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    to_tsv_required = to_tsv_parser.add_argument_group("Required Parameters")
    to_tsv_optional = to_tsv_parser.add_argument_group("Optional Parameters")

    to_tsv_required.add_argument(
        "-i",
        "--input-lineage",
        metavar="<FILE>",
        help="Input lineage file",
        required=True
    )

    to_tsv_optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output standard taxonomy tsv file (default: formatted-tax.tsv)',
        default="formatted-tax.tsv"
    )

    to_tsv_optional.add_argument(
        "--make-taxid",
        help="Provide this flag to make a unique taxid (string of all rank fields) for each lineage " \
             "(will be added as second column of output).",
        action="store_true"
    )

    add_help(to_tsv_optional)

    add_version_arg(to_tsv_optional)

    to_tsv_parser.set_defaults(func="to_tsv")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so appropriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in ("-v", "--version"):
            from bit.modules.general import report_version
            report_version()
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            if sys.stdin.isatty():
                parser.subparsers.choices[cmd].print_help(sys.stderr)
                sys.exit(0)
            # else: stdin is being piped, fall through to parse_args()
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    # check input files exist
    if hasattr(args, 'input_file'):
        check_files_are_found([args.input_file])
    elif hasattr(args, 'input_lineage'):
        check_files_are_found([args.input_lineage])

    # call appropriate function
    if args.func == "from_taxids":
        from bit.modules.ncbi.get_lineage_from_taxids import get_lineage_from_taxids
        get_lineage_from_taxids(args.input_file, args.output_file, args.include_strain)
    elif args.func == "to_tsv":
        from bit.modules.lineage_to_tsv import convert_lineage_to_tsv
        convert_lineage_to_tsv(args.input_lineage, args.output_tsv, args.make_taxid)
