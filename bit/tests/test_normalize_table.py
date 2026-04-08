import pytest # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
from bit.tests.utils import run_cli
from bit.cli.normalize_table import (
    normalize_cpm,
    normalize_median_ratio,
    remove_zero_columns,
    restore_zero_columns,
)


@pytest.fixture
def sample_table(tmp_path):
    """A small tab-delimited table with 3 genes x 3 samples."""
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
def table_with_zero_column(tmp_path):
    """Table with one column that sums to zero."""
    text = (
        "gene\tsampleA\tsampleB\tsampleC\n"
        "g1\t10\t0\t30\n"
        "g2\t40\t0\t60\n"
        "g3\t50\t0\t10\n"
    )
    p = tmp_path / "input_zero.tsv"
    p.write_text(text)
    return p


@pytest.fixture
def basic_df():
    return pd.DataFrame(
        {"s1": [5, 10, 15], "s2": [50, 100, 150], "s3": [10, 20, 30]},
        index=["g1", "g2", "g3"],
    )


# CPM tests
def test_cpm_columns_sum_to_million(basic_df):
    result = normalize_cpm(basic_df)
    for col in result.columns:
        assert result[col].sum() == pytest.approx(1_000_000)


def test_cpm_values(basic_df):
    result = normalize_cpm(basic_df)
    expected = {
        "g1": 166666.66666666666,
        "g2": 333333.3333333333,
        "g3": 500000.0,
    }
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


# MR test
def test_mr_values(basic_df):
    result = normalize_median_ratio(basic_df)
    expected = {
        "g1": 13.572088082974535,
        "g2": 27.14417616594907,
        "g3": 40.71626424892361,
    }
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


# removing/restoring zero columns tests
def test_remove_zero_columns():
    df = pd.DataFrame(
        {"s1": [1, 2], "s2": [0, 0], "s3": [3, 4]},
        index=["g1", "g2"],
    )
    cleaned, zeros, ordered = remove_zero_columns(df)
    assert "s2" not in cleaned.columns
    assert zeros == ["s2"]
    assert ordered == ["s1", "s2", "s3"]


def test_restore_zero_columns():
    df = pd.DataFrame(
        {"s1": [1.0, 2.0], "s3": [3.0, 4.0]},
        index=["g1", "g2"],
    )
    restored = restore_zero_columns(df, ["s2"], ["s1", "s2", "s3"])
    assert list(restored.columns) == ["s1", "s2", "s3"]
    assert (restored["s2"] == 0.0).all()


# cli tests
def test_cli_cpm(sample_table, tmp_path):
    out = tmp_path / "out_cpm.tsv"
    cmd = [
        "bit-normalize-table",
        "-i", str(sample_table),
        "-n", "CPM",
        "-o", str(out),
    ]
    run_cli(cmd)
    result = pd.read_csv(out, sep="\t", index_col=0)
    expected = {
        "g1": 166666.66666666666,
        "g2": 333333.3333333333,
        "g3": 500000.0,
    }
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)


def test_cli_mr(sample_table, tmp_path):
    out = tmp_path / "out_mr.tsv"
    cmd = [
        "bit-normalize-table",
        "-i", str(sample_table),
        "-n", "MR",
        "-o", str(out),
    ]
    run_cli(cmd)
    result = pd.read_csv(out, sep="\t", index_col=0)
    expected = {
        "g1": 13.572088082974535,
        "g2": 27.14417616594907,
        "g3": 40.71626424892361,
    }
    for gene, val in expected.items():
        for col in result.columns:
            assert result.loc[gene, col] == pytest.approx(val)
