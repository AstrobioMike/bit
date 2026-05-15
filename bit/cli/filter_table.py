import sys
import gzip
import argparse
from bit.modules.general import check_files_are_found, report_message
from bit.cli.common import CustomRichHelpFormatter, add_help


def build_parser():

    desc = """
        Filters a table based on values in a specified column. For version info,
        run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-filter-table -i input.tsv -w wanted-values.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-table",
        metavar="<FILE>",
        help="Input table",
        required=True,
    )

    required.add_argument(
        "-w",
        "--wanted-values",
        metavar="<FILE>",
        help="Single-column file of wanted values",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output table filename (default: "filtered.tsv")',
        default="filtered.tsv",
    )

    optional.add_argument(
        "-d",
        "--delimiter",
        metavar="<STR>",
        help='Delimiter (default: "\\t")',
        default="\t",
    )

    optional.add_argument(
        "-c",
        "--column",
        metavar="<INT>",
        help="Index of column to filter on, 1-based (default: 1)",
        default=1,
        type=int,
    )

    optional.add_argument(
        "--no-header",
        help="Add if there is no header",
        action="store_true",
    )

    optional.add_argument(
        "--gz",
        help="Add if the input is gzipped (output will not be)",
        action="store_true",
    )

    add_help(optional)

    return parser


def main(args=None):

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args(args)

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


def filter_table(input_table, wanted_file, output_file, delimiter, column, no_header, gz):

    targets = set(line.strip() for line in open(wanted_file))

    target_column = column - 1

    if gz:
        input_fh = gzip.open(input_table, "rt")
    else:
        input_fh = open(input_table, "r")

    with input_fh, open(output_file, "w") as output:

        num_targets_found = 0

        if not no_header:
            first_line = True

            for line in input_fh:
                if first_line:
                    output.write(line)
                    first_line = False
                    continue

                split_line = line.strip().split(delimiter)

                if split_line[target_column] in targets:
                    num_targets_found += 1
                    output.write(line)

        else:

            for line in input_fh:
                split_line = line.strip().split(delimiter)

                if split_line[target_column] in targets:
                    num_targets_found += 1
                    output.write(line)


    if num_targets_found == 0:
        import os
        report_message(f"None of the wanted values were found in column {column}.", trailing_newline=True)
        os.remove(output_file)
    else:
        report_message(f"Found {num_targets_found} row(s) with values in column {column} that were in the wanted-values file.")
        report_message(f"Output tsv written to {output_file}.", trailing_newline=True, color="none")


if __name__ == "__main__":
    main()
