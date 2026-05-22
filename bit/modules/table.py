import gzip


# ── filter ───────────────────────────────────────────────────────────────────

def filter_table(input_table, wanted_file, output_file, delimiter, column, no_header, gz):

    from bit.modules.general import report_message

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


# ── normalize ────────────────────────────────────────────────────────────────

def remove_zero_columns(tab):
    ordered_columns = tab.columns.tolist()
    column_sums = tab.sum()
    zero_column_names = column_sums[column_sums == 0].index.tolist()
    tab = tab.drop(zero_column_names, axis=1)
    return tab, zero_column_names, ordered_columns


def restore_zero_columns(tab, zero_column_names, ordered_columns):
    for col in zero_column_names:
        tab[col] = 0.0
    return tab[ordered_columns]


def normalize_cpm(tab):
    return tab / tab.sum() * 1000000


def normalize_median_ratio(tab):

    import numpy as np # type: ignore
    from scipy.stats.mstats import gmean # type: ignore

    # getting geometric means for each row
    with np.errstate(divide='ignore'):
        geomeans = gmean(tab, axis=1)

    # getting ratios of values to geometric means
    ratios_tab = (tab.T / geomeans).T

    # calculating size factors from rows with non-zero geometric means
    size_factors = ratios_tab[geomeans > 0].median().to_list()

    return tab / size_factors


# ── summarize-column ─────────────────────────────────────────────────────────

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
            from bit.modules.general import report_failure
            message = "It seems you specified a column by name, but the input file doesn't appear to have a header."
            report_failure(message)

    return header


def get_series_by_name(df, column_name):
    if column_name not in df.columns:
        from bit.modules.general import report_failure
        message = f"The specified column '{column_name}' does not exist in the input file."
        report_failure(message)
    return df[column_name]


def print_summary(series, column):
    from bit.modules.general import color_text

    target_percentiles = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    sum = series.sum()
    desc = series.describe(percentiles=target_percentiles).round(2)
    desc['std'] = series.std(ddof=0).round(2)

    def fmt(val):
        s = f"{val:.2f}"
        return s.rstrip("0").rstrip(".")

    print(color_text(f"\n  Column '{column}' summary\n", "yellow"))
    rows = [
        ("N:",      int(desc["count"])),
        ("Min:",    desc["min"]),
        ("Max:",    desc["max"]),
        ("Sum:",    sum),
        ("Mean:",   desc["mean"]),
        ("Median:", desc["50%"]),
        ("StDev:",  desc["std"]),
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
