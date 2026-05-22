import sys
import io
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        This program has utilities for working with tabular data. See subcommand-specific
        help menus for more info. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers


    ### colnames subcommand ###

    colnames_desc = """
        Returns the column names (numbered) from a delimited file. The delimiter is
        auto-detected for tab, comma, pipe, semicolon, or space.
        """

    colnames_parser = subparsers.add_parser(
        "colnames",
        help="List column names of a delimited file",
        description=colnames_desc,
        epilog="Ex. usage: `bit-table colnames input.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    colnames_required = colnames_parser.add_argument_group("Required Parameters")
    colnames_optional = colnames_parser.add_argument_group("Optional Parameters")

    colnames_required.add_argument(
        "input_file",
        metavar="<FILE>",
        nargs="?",
        default=sys.stdin,
        help="Input delimited file or stdin if none provided"
    )

    add_help(colnames_optional)
    colnames_parser.set_defaults(func="colnames")


    ### filter subcommand ###

    filter_desc = """
        Filters a table based on provided values (strings) in a specified column.
        """

    filter_parser = subparsers.add_parser(
        "filter",
        help="Filter a table by values (strings) in a column",
        description=filter_desc,
        epilog="Ex. usage: `bit-table filter -i input.tsv -w wanted-values.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    filter_required = filter_parser.add_argument_group("Required Parameters")
    filter_optional = filter_parser.add_argument_group("Optional Parameters")

    filter_required.add_argument(
        "-i",
        "--input-table",
        metavar="<FILE>",
        help="Input table",
        required=True,
    )

    filter_required.add_argument(
        "-w",
        "--wanted-values",
        metavar="<FILE>",
        help="Single-column file of wanted values (strings)",
        required=True,
    )

    filter_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output table filename (default: "filtered.tsv")',
        default="filtered.tsv",
    )

    filter_optional.add_argument(
        "-d",
        "--delimiter",
        metavar="<STR>",
        help='Delimiter (default: "\\t")',
        default="\t",
    )

    filter_optional.add_argument(
        "-c",
        "--column",
        metavar="<INT>",
        help="Index of column to filter on, 1-based (default: 1)",
        default=1,
        type=int,
    )

    filter_optional.add_argument(
        "--no-header",
        help="Add if there is no header",
        action="store_true",
    )

    filter_optional.add_argument(
        "--gz",
        help="Add if the input is gzipped (output will not be)",
        action="store_true",
    )

    add_help(filter_optional)
    filter_parser.set_defaults(func="filter")


    ### normalize subcommand ###

    normalize_desc = """
        Normalizes a tab-delimited table by either counts- or coverage-per-million (CPM) or
        with the median-ratio method as performed in DESeq2. Expects samples as columns and
        units (e.g. genes/KOs/OTUs/etc.) as rows.
        """

    normalize_parser = subparsers.add_parser(
        "normalize",
        help="Normalize a table by CPM or median-ratio",
        description=normalize_desc,
        epilog='Ex. usage: `bit-table normalize -i input-table.tsv -n CPM`',
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    normalize_required = normalize_parser.add_argument_group("Required Parameters")
    normalize_optional = normalize_parser.add_argument_group("Optional Parameters")

    normalize_required.add_argument(
        "-i",
        "--input-table",
        metavar="<FILE>",
        help="Input tab-delimited table",
        required=True,
    )

    normalize_optional.add_argument(
        "-n",
        "--normalization",
        help='Normalization method: "CPM" (counts/coverage-per-million) or "MR" (median-ratio) (default: "CPM")',
        choices=["CPM", "MR"],
        default="CPM",
    )

    normalize_optional.add_argument(
        "-o",
        "--output-table",
        metavar="<FILE>",
        help='Output filename (default: "normalized.tsv")',
        default="normalized.tsv",
    )

    add_help(normalize_optional)
    normalize_parser.set_defaults(func="normalize")


    ### summarize-column subcommand ###

    summarize_col_desc = """
        Outputs general summary stats for a numeric column. Can take stdin or a file as
        input. Runs on the first column by default.
        """

    summarize_col_parser = subparsers.add_parser(
        "summarize-column",
        help="Summarize stats of a numeric column",
        description=summarize_col_desc,
        epilog="Ex. usage: `bit-table summarize-column data.tsv -c 2`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    summarize_col_required = summarize_col_parser.add_argument_group("Required Parameters")
    summarize_col_optional = summarize_col_parser.add_argument_group("Optional Parameters")

    summarize_col_required.add_argument(
        "input_file",
        metavar="<FILE>",
        nargs="?",
        default=sys.stdin,
        help="Input file or stdin if none provided"
    )

    summarize_col_optional.add_argument(
        "-c",
        "--column",
        metavar="<STR>",
        help="Target column: 1-based index or column name (default: 1)",
        default=1,
    )

    summarize_col_optional.add_argument(
        "-d",
        "--delimiter",
        metavar="<STR>",
        help="Delimiter (default: '\\t')",
        default="\t",
    )

    add_help(summarize_col_optional)
    summarize_col_parser.set_defaults(func="summarize_column")


    return parser


################################################################################


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv) == 1 and sys.stdin.isatty(): # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so appropriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
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

    # stdin handling for subcommands that accept it (colnames, summarize-column)
    if hasattr(args, 'input_file') and args.input_file is sys.stdin:
        text = args.input_file.read()
        args.input_file = io.StringIO(text)
    elif hasattr(args, 'input_file'):
        from bit.modules.general import check_files_are_found
        check_files_are_found([args.input_file])

    func_map = {
        "colnames"        : _run_colnames,
        "filter"          : _run_filter,
        "normalize"       : _run_normalize,
        "summarize_column": _run_summarize_column,
    }

    func_map[args.func](args)


def _run_colnames(args):
    from bit.modules.general import colnames
    colnames(args)


def _run_filter(args):
    from bit.modules.general import check_files_are_found
    from bit.cli.filter_table import filter_table
    check_files_are_found([args.input_table, args.wanted_values])
    filter_table(
        input_table=args.input_table,
        wanted_file=args.wanted_values,
        output_file=args.output_file,
        delimiter=args.delimiter,
        column=args.column,
        no_header=args.no_header,
        gz=args.gz,
    )


def _run_normalize(args):
    from bit.modules.general import check_files_are_found
    from bit.cli.normalize_table import (remove_zero_columns, normalize_cpm,
                                         normalize_median_ratio, restore_zero_columns)
    import pandas as pd # type: ignore
    check_files_are_found([args.input_table])
    tab = pd.read_csv(args.input_table, sep="\t", index_col=0, low_memory=False)
    tab, zero_column_names, ordered_columns = remove_zero_columns(tab)
    norm_tab = normalize_cpm(tab) if args.normalization == "CPM" else normalize_median_ratio(tab)
    norm_tab = restore_zero_columns(norm_tab, zero_column_names, ordered_columns)
    norm_tab.to_csv(args.output_table, sep="\t")


def _run_summarize_column(args):
    from bit.cli.summarize_column import summarize_column
    input_file = args.input_file
    if isinstance(input_file, str):
        input_file = open(input_file, "r")
    summarize_column(input_file=input_file, column=args.column, delimiter=args.delimiter)
