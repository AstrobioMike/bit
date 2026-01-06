import io
import pytest
import csv
from unittest.mock import patch
from bit.tests.utils import run_cli
from bit.cli.summarize_column import detect_header


@pytest.fixture
def sample_table(tmp_path):
    text = "\t".join(["mike","lee"]) + "\n" + \
           "\n".join(
               "\t".join(map(str,row))
               for row in [(2,20),(3,30),(5,40)]
           ) + "\n"
    p = tmp_path/"data.tsv"
    p.write_text(text)
    return p


def parse_lines(output):
    """
    Split into non-empty lines, strip right whitespace
    """
    return [line.rstrip() for line in output.splitlines() if line.strip()]


def test_summarize_first_column_from_stdin(sample_table):
    cmd = [
        "bit-summarize-column",
        "-c", "1",
    ]
    data = sample_table.read_text()
    result = run_cli(cmd, input=data)

    lines = parse_lines(result.stdout)

    assert lines[0] == "  Column '1' summary"

    stats = {
        "N:": "3",
        "Min:": "2",
        "Max:": "5",
        "Sum:": "10",
        "Mean:": "3",
        "Median:": "3",
        "StDev:": "1.25"
    }
    for label, expect in stats.items():
        found = [L for L in lines if L.startswith(f"    {label}")]
        assert found, f"Did not find summary line for {label}"
        assert found[0].endswith(expect), f"{label!r} line was {found[0]!r}, expected to end with {expect!r}"


def test_summarize_second_column_by_name(sample_table):
    cmd = [
        "bit-summarize-column",
        "-i", str(sample_table),
        "-c", "lee",
    ]
    result = run_cli(cmd)

    lines = parse_lines(result.stdout)

    assert lines[0] == "  Column 'lee' summary"

    stats = {
        "N:": "3",
        "Min:": "20",
        "Max:": "40",
        "Sum:": "90",
        "Mean:": "30",
        "Median:": "30",
        "StDev:": "8.16"
    }
    for label, expect in stats.items():
        found = [L for L in lines if L.startswith(f"    {label}")]
        assert found, f"Did not find summary line for {label}"
        assert found[0].endswith(expect), f"{label!r} line was {found[0]!r}, expected to end with {expect!r}"


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
