import sys
import argparse
from bit.modules.general import report_failure, color_text
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)


def build_parser():

    desc = """
        This script outputs general summary stats for a numeric column. It can take stdin or
        a file as input. It will run on the first (or only column) if not specified. Otherwise
        you can indicate which column to summarize by column position or name. For version info,
        run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-summarize-column data.tsv -c 2",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )
    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "input_file",
        metavar="<FILE>",
        nargs='?',
        default=sys.stdin,
        help="Input file or stdin if none provided"
    )
    optional.add_argument(
        "-c",
        "--column",
        metavar="<STR>",
        help="Specify the target column to summarize. Can be a number specifying the column index (1-based, like unix cut/awk), \
              or can be a column name. (default: 1)",
        default=1,
    )
    optional.add_argument(
        "-d",
        "--delimiter",
        metavar="<STR>",
        help="Specify the delimiter (default is a tab: '\\t')",
        default="\t",
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.input_file is sys.stdin:
        import io
        text = args.input_file.read()
        args.input_file = io.StringIO(text)
    else:
        args.input_file = open(args.input_file, 'r')

    summarize_column(
        input_file=args.input_file,
        column=args.column,
        delimiter=args.delimiter
    )


def summarize_column(input_file, column, delimiter):
    target_series = load_series(input_file, column, delimiter)
    print_summary(target_series, column)


def load_series(input_file, column, delimiter):
    header = detect_header(input_file, column)
    import pandas as pd # type: ignore

    df = pd.read_csv(
        input_file,
        sep=delimiter,
        header=0 if header else None
    )
    try:
        col_num = int(column)
        series = df.iloc[:, col_num - 1]
    except ValueError:
        series = get_series_by_name(df, column)
    return pd.to_numeric(series, errors="coerce").dropna()


def detect_header(input_file, column):
    sample = input_file.read(4096)
    input_file.seek(0)
    import csv

    try:
        header = csv.Sniffer().has_header(sample)
    except csv.Error:
        # single-column / ambiguous case: if column isn't an int, treat first row as header
        try:
            int(column)
            header = False
        except ValueError:
            header = True

    if not header:
        try:
            column = int(column)
        except ValueError:
            message = "It seems you specified a column by name, but the input file doesn't appear to have a header."
            report_failure(message)

    return header


def get_series_by_name(df, column_name):
    if column_name not in df.columns:
        message = f"The specified column '{column_name}' does not exist in the input file."
        report_failure(message)
    return df[column_name]


def print_summary(series, column):
    target_percentiles = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    sum = series.sum()
    desc = series.describe(percentiles=target_percentiles).round(2)
    desc['std'] = series.std(ddof=0).round(2)

    def fmt(val):
        s = f"{val:.2f}"
        return s.rstrip("0").rstrip(".")

    print(color_text(f"\n  Column '{column}' summary\n", "yellow"))
    rows = [
        ("N:",     int(desc["count"])),
        ("Min:",   desc["min"]),
        ("Max:",   desc["max"]),
        ("Sum:",   sum),
        ("Mean:",  desc["mean"]),
        ("Median:",desc["50%"]),
        ("StDev:", desc["std"]),
    ]
    for name, val in rows:
        value = str(val) if isinstance(val, int) else fmt(val)
        print(f"    {name:<{16}}{value}")
    print(color_text("\n    Percentiles:\n", "yellow"))
    for percentile in target_percentiles:
        pct = int(percentile * 100)
        suffix = "st" if pct == 1 else "th"
        label = f"{pct}{suffix}:"
        key = f"{pct}%"
        value = desc[key]
        value = f"{value:.2f}"
        print(f"        {label:<{12}}{value}")

    print()
