import gzip
import io
import csv
import pytest # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from unittest.mock import patch
from bit.tests.utils import run_cli
from bit.modules.table import (
    filter_table,
    normalize_cpm,
    normalize_median_ratio,
    remove_zero_columns,
    restore_zero_columns,
    detect_header,
)


# ── filter fixtures ───────────────────────────────────────────────────────────

@pytest.fixture
def filter_table_tsv(tmp_path):
    text = (
        "id\tvalue\tcategory\n"
        "AAA\t10\talpha\n"
        "BBB\t20\tbeta\n"
        "CCC\t30\talpha\n"
        "DDD\t40\tgamma\n"
    )
    p = tmp_path / "table.tsv"
    p.write_text(text)
    return p


@pytest.fixture
def filter_table_no_header(tmp_path):
    text = (
        "AAA\t10\talpha\n"
        "BBB\t20\tbeta\n"
        "CCC\t30\talpha\n"
        "DDD\t40\tgamma\n"
    )
    p = tmp_path / "table_no_header.tsv"
    p.write_text(text)
    return p


@pytest.fixture
def filter_table_gz(tmp_path):
    text = (
        "id\tvalue\tcategory\n"
        "AAA\t10\talpha\n"
        "BBB\t20\tbeta\n"
        "CCC\t30\talpha\n"
        "DDD\t40\tgamma\n"
    )
    p = tmp_path / "table.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write(text)
    return p


@pytest.fixture
def wanted_col1(tmp_path):
    p = tmp_path / "wanted.txt"
    p.write_text("AAA\nCCC\n")
    return p


@pytest.fixture
def wanted_col3(tmp_path):
    p = tmp_path / "wanted_col3.txt"
    p.write_text("alpha\n")
    return p


# ── filter unit tests ─────────────────────────────────────────────────────────

def test_basic_filter_col1(filter_table_tsv, wanted_col1, tmp_path):
    out = tmp_path / "out.tsv"
    filter_table(str(filter_table_tsv), str(wanted_col1), str(out), "\t", 1, False, False)
    lines = out.read_text().splitlines()
    assert lines[0] == "id\tvalue\tcategory"   # header preserved
    assert len(lines) == 3                     # header + 2 matching rows
    assert [l.split("\t")[0] for l in lines[1:]] == ["AAA", "CCC"]


def test_filter_on_column_3(filter_table_tsv, wanted_col3, tmp_path):
    out = tmp_path / "out.tsv"
    filter_table(str(filter_table_tsv), str(wanted_col3), str(out), "\t", 3, False, False)
    lines = out.read_text().splitlines()
    assert lines[0] == "id\tvalue\tcategory"
    assert len(lines) == 3                     # header + AAA + CCC (both alpha)
    assert all(l.split("\t")[2] == "alpha" for l in lines[1:])


def test_filter_no_header(filter_table_no_header, wanted_col1, tmp_path):
    out = tmp_path / "out.tsv"
    filter_table(str(filter_table_no_header), str(wanted_col1), str(out), "\t", 1, True, False)
    lines = out.read_text().splitlines()
    assert len(lines) == 2
    assert [l.split("\t")[0] for l in lines] == ["AAA", "CCC"]


def test_filter_gz_input(filter_table_gz, wanted_col1, tmp_path):
    out = tmp_path / "out.tsv"
    filter_table(str(filter_table_gz), str(wanted_col1), str(out), "\t", 1, False, True)
    lines = out.read_text().splitlines()
    assert lines[0] == "id\tvalue\tcategory"
    assert len(lines) == 3
    assert [l.split("\t")[0] for l in lines[1:]] == ["AAA", "CCC"]


def test_no_matches_removes_output(filter_table_tsv, tmp_path):
    wanted = tmp_path / "nomatch.txt"
    wanted.write_text("ZZZ\n")
    out = tmp_path / "out.tsv"
    filter_table(str(filter_table_tsv), str(wanted), str(out), "\t", 1, False, False)
    assert not out.exists()


# ── filter cli tests ──────────────────────────────────────────────────────────

def test_cli_filter_basic(filter_table_tsv, wanted_col1, tmp_path):
    out = tmp_path / "out.tsv"
    run_cli(["bit", "table", "filter", "-i", str(filter_table_tsv), "-w", str(wanted_col1), "-o", str(out)])
    lines = out.read_text().splitlines()
    assert lines[0] == "id\tvalue\tcategory"
    assert len(lines) == 3
    assert [l.split("\t")[0] for l in lines[1:]] == ["AAA", "CCC"]


def test_cli_filter_column_flag(filter_table_tsv, wanted_col3, tmp_path):
    out = tmp_path / "out.tsv"
    run_cli(["bit", "table", "filter", "-i", str(filter_table_tsv), "-w", str(wanted_col3), "-o", str(out), "-c", "3"])
    lines = out.read_text().splitlines()
    assert len(lines) == 3
    assert all(l.split("\t")[2] == "alpha" for l in lines[1:])


# ── normalize fixtures ────────────────────────────────────────────────────────

@pytest.fixture
def norm_table(tmp_path):
    text = (
        "gene\ts1\ts2\ts3\n"
        "g1\t5\t50\t10\n"
        "g2\t10\t100\t20\n"
        "g3\t15\t150\t30\n"
    )
    p = tmp_path / "input.tsv"
    p.write_text(text)
    return p


@pytest.fixture
def basic_df():
    return pd.DataFrame(
        {"s1": [5, 10, 15], "s2": [50, 100, 150], "s3": [10, 20, 30]},
        index=["g1", "g2", "g3"],
    )


# ── normalize unit tests ──────────────────────────────────────────────────────

def test_cpm_columns_sum_to_million(basic_df):
    result = normalize_cpm(basic_df)
    for col in result.columns:
        assert result[col].sum() == pytest.approx(1_000_000)


def test_cpm_values(basic_df):
    result = normalize_cpm(basic_df)
    expected = {"g1": 166666.66666666666, "g2": 333333.3333333333, "g3": 500000.0}
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


def test_mr_values(basic_df):
    result = normalize_median_ratio(basic_df)
    expected = {"g1": 13.572088082974535, "g2": 27.14417616594907, "g3": 40.71626424892361}
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


def test_remove_zero_columns():
    df = pd.DataFrame({"s1": [1, 2], "s2": [0, 0], "s3": [3, 4]}, index=["g1", "g2"])
    cleaned, zeros, ordered = remove_zero_columns(df)
    assert "s2" not in cleaned.columns
    assert zeros == ["s2"]
    assert ordered == ["s1", "s2", "s3"]


def test_restore_zero_columns():
    df = pd.DataFrame({"s1": [1.0, 2.0], "s3": [3.0, 4.0]}, index=["g1", "g2"])
    restored = restore_zero_columns(df, ["s2"], ["s1", "s2", "s3"])
    assert list(restored.columns) == ["s1", "s2", "s3"]
    assert (restored["s2"] == 0.0).all()


# ── normalize cli tests ───────────────────────────────────────────────────────

def test_cli_cpm(norm_table, tmp_path):
    out = tmp_path / "out_cpm.tsv"
    run_cli(["bit", "table", "normalize", "-i", str(norm_table), "-n", "CPM", "-o", str(out)])
    result = pd.read_csv(out, sep="\t", index_col=0)
    expected = {"g1": 166666.66666666666, "g2": 333333.3333333333, "g3": 500000.0}
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


def test_cli_mr(norm_table, tmp_path):
    out = tmp_path / "out_mr.tsv"
    run_cli(["bit", "table", "normalize", "-i", str(norm_table), "-n", "MR", "-o", str(out)])
    result = pd.read_csv(out, sep="\t", index_col=0)
    expected = {"g1": 13.572088082974535, "g2": 27.14417616594907, "g3": 40.71626424892361}
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


# ── summarize-column fixtures ─────────────────────────────────────────────────

@pytest.fixture
def summarize_table(tmp_path):
    text = "\t".join(["mike", "lee"]) + "\n" + \
           "\n".join("\t".join(map(str, row)) for row in [(2, 20), (3, 30), (5, 40)]) + "\n"
    p = tmp_path / "data.tsv"
    p.write_text(text)
    return p


def parse_summary_lines(output):
    return [line.rstrip() for line in output.splitlines() if line.strip()]


# ── summarize-column unit tests ───────────────────────────────────────────────

def test_detect_header_sniffer_error_numeric_column():
    fake_file = io.StringIO("1\n2\n3\n")
    with patch("csv.Sniffer.has_header", side_effect=csv.Error):
        header = detect_header(fake_file, column="1")
    assert header is False


def test_detect_header_sniffer_error_named_column():
    fake_file = io.StringIO("value\n1\n2\n3\n")
    with patch("csv.Sniffer.has_header", side_effect=csv.Error):
        header = detect_header(fake_file, column="value")
    assert header is True


def test_detect_header_no_header_named_column_fails():
    fake_file = io.StringIO("1\n2\n3\n")
    with patch("csv.Sniffer.has_header", return_value=False):
        with pytest.raises(SystemExit) as e:
            detect_header(fake_file, column="col1")
    assert e.value.code == 1


# ── summarize-column cli tests ────────────────────────────────────────────────

def test_summarize_first_column_from_stdin(summarize_table):
    result = run_cli(["bit", "table", "summarize-column", "-c", "1"], input=summarize_table.read_text())
    lines = parse_summary_lines(result.stdout)
    assert lines[0] == "  Column '1' summary"
    stats = {"N:": "3", "Min:": "2", "Max:": "5", "Sum:": "10", "Mean:": "3", "Median:": "3", "StDev:": "1.25"}
    for label, expect in stats.items():
        found = [L for L in lines if L.startswith(f"    {label}")]
        assert found, f"Did not find summary line for {label}"
        assert found[0].endswith(expect), f"{label!r} line was {found[0]!r}, expected to end with {expect!r}"


def test_summarize_second_column_by_name(summarize_table):
    result = run_cli(["bit", "table", "summarize-column", str(summarize_table), "-c", "lee"])
    lines = parse_summary_lines(result.stdout)
    assert lines[0] == "  Column 'lee' summary"
    stats = {"N:": "3", "Min:": "20", "Max:": "40", "Sum:": "90", "Mean:": "30", "Median:": "30", "StDev:": "8.16"}
    for label, expect in stats.items():
        found = [L for L in lines if L.startswith(f"    {label}")]
        assert found, f"Did not find summary line for {label}"
        assert found[0].endswith(expect), f"{label!r} line was {found[0]!r}, expected to end with {expect!r}"
