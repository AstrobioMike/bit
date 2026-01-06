import shutil
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
import bit.cli.summarize_assembly

test_assembly = get_package_path("tests/data/ez-screen-targets.fasta")

def test_summarize_assembly(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    cmd = [
        "bit-summarize-assembly",
        str(test_assembly),
        "-o", f"{out_dir}/test.tsv"
    ]

    run_cli(cmd)

    expected = {
        "Assembly":                 "ez-screen-targets",
        "Total contigs":            "2",
        "Total length":             "1209",
        "Ambiguous characters":     "0",
        "GC content":               "43.09",
        "Maximum contig length":    "660",
        "Minimum contig length":    "549",
        "N50":                      "660",
        "N75":                      "549",
        "N90":                      "549",
        "L50":                      "1",
        "L75":                      "2",
        "L90":                      "2",
        "Num. contigs >= 100":      "2",
        "Num. contigs >= 500":      "2",
        "Num. contigs >= 1000":     "0",
        "Num. contigs >= 5000":     "0",
        "Num. contigs >= 10000":    "0",
        "Num. contigs >= 50000":    "0",
        "Num. contigs >= 100000":   "0",
    }

    summary_tsv = out_dir / "test.tsv"
    assert summary_tsv.exists(), f"Summary TSV not found at {summary_tsv}"
    text = summary_tsv.read_text().splitlines()

    for line in text[0:]:
        key, value = line.split("\t", 1)
        assert key in expected, f"Unexpected key in summary TSV: {key}"
        assert value == expected[key], f"Unexpected value for {key} in summary TSV: {value} (expected {expected[key]})"
