import pandas as pd
import numpy as np
import pytest
from bit.modules import kraken2_to_taxon_summaries as k


@pytest.fixture
def report_min(tmp_path):
    """
    Tiny kraken2-like report covering:
      - unclassified (U)
      - root (R)
      - domain given as R1 with Bacteria (to test normalize_rank_code)
      - genus lines with spaces around names (to test stripping)
      - one non-standard rank '-' which should produce a row but not change lineage
    """
    text = "\n".join([
        # percent clade_reads taxon_reads rank taxid name(with indent)
        "10.00 100 100 U 0 unclassified",
        "90.00 900 50 R 1 root",
        "80.00 800 0 R1 2 Bacteria",
        "25.00 250 250 G 123    GenX",     # leading spaces
        "15.00 150 150 G 124 GenY   ",     # trailing spaces
        " 0.50   5   5 -  9999    12345_like_strain",  # non-standard rank
        "",  # blank line should be ignored
    ]) + "\n"
    p = tmp_path / "k.report"
    p.write_text(text)
    return p


def test_parse_report_line_basic_space_split():
    line = "12.3 120 12 G 1234    Some Genus"
    rec = k.parse_report_line(line)
    assert rec == {
        "clade_reads": 120,
        "taxon_reads": 12,
        "rank": "G",
        "taxid": 1234,
        "name": "Some Genus",
    }

def test_parse_report_line_tab_fallback():
    line = "12.3\t120\t12\tG\t1234\t   Some\tGenus  "
    rec = k.parse_report_line(line)
    assert rec["taxid"] == 1234
    assert rec["rank"] == "G"
    assert rec["name"] == "Some\tGenus"

def test_normalize_rank_code_domain_R1():
    assert k.normalize_rank_code("R1", "Bacteria") == "D"
    assert k.normalize_rank_code("R1", "Viruses") == "D"
    assert k.normalize_rank_code("R1", "WeirdDomain") == "R1"
    assert k.normalize_rank_code("G", "Genus") == "G"


def test_parse_report_builds_lineages(report_min):
    df = k.parse_report(str(report_min))

    un = df[df["taxid"] == 0].iloc[0]
    assert all(un[r] == "Unclassified" for r in k.STD_RANKS)
    assert un["read_counts"] == 100

    gx = df[df["taxid"] == 123].iloc[0]
    gy = df[df["taxid"] == 124].iloc[0]
    assert gx["genus"] == "GenX"
    assert gy["genus"] == "GenY"
    assert gx["domain"] == "Bacteria"
    assert gy["domain"] == "Bacteria"

    row_dash = df[df["taxid"] == 9999].iloc[0]
    assert row_dash["domain"] == "Bacteria"
    assert row_dash["domain"] != "Unclassified"


def test_refine_df_aggregates_duplicates_and_computes_percents():
    rows = [
        {"taxid": 5, "domain":"Bacteria","phylum":"Firmicutes","class":"Bacilli","order":"NA","family":"NA","genus":"GenZ","species":"NA","read_counts":30},
        {"taxid": 5, "domain":"Bacteria","phylum":"Firmicutes","class":"Bacilli","order":"NA","family":"NA","genus":"GenZ","species":"NA","read_counts":70},
        {"taxid": 0, **{r:"Unclassified" for r in k.STD_RANKS}, "read_counts":100},
    ]
    df = pd.DataFrame(rows)
    out = k.refine_df(df.copy())

    agg = out[out["taxid"] == 5].iloc[0]
    assert agg["read_counts"] == 100

    z = out.set_index("taxid")["percent_of_reads"].to_dict()
    assert np.isclose(z[5], 50.0, atol=1e-6)
    assert np.isclose(z[0], 50.0, atol=1e-6)

def test_refine_df_zero_total_reads():
    rows = [
        {"taxid": 1, **{r:"NA" for r in k.STD_RANKS}, "read_counts":0},
        {"taxid": 2, **{r:"NA" for r in k.STD_RANKS}, "read_counts":0},
    ]
    df = pd.DataFrame(rows)
    out = k.refine_df(df.copy())
    # zero-read rows are dropped
    assert out.empty


def test_sort_df_custom_groups_and_stability():
    # building rows to exercise sort groups:
    #  - Unclassified (taxid=0) should come first
    #  - All-NA lineage (group 1) second
    #  - The rest (group 2) sorted by lineage then taxid (mergesort: stable)
    rows = [
        {"taxid": 3, "domain":"Bacteria","phylum":"Actino","class":"C1","order":"O1","family":"F1","genus":"G1","species":"S1","read_counts":10},
        {"taxid": 0, **{r:"Unclassified" for r in k.STD_RANKS}, "read_counts":10},
        {"taxid": 2, **{r:"NA" for r in k.STD_RANKS}, "read_counts":10},
        {"taxid": 4, "domain":"Archaea","phylum":"Eury","class":"C1","order":"O1","family":"F1","genus":"G1","species":"S1","read_counts":10},
        {"taxid": 5, "domain":"Archaea","phylum":"Eury","class":"C1","order":"O1","family":"F1","genus":"G1","species":"S1","read_counts":10},
    ]
    df = pd.DataFrame(rows)
    out = k.sort_df(df.copy())

    # expected order by sort_group then lineage then taxid:
    # 0 first, then 2, then archaea 4,5 (taxid ascending), then bacteria 3
    assert list(out["taxid"]) == [0, 2, 4, 5, 3]


def test_kraken2_to_taxon_summaries_writes_tsv(tmp_path, report_min, monkeypatch):
    # avoid touching the filesystem checker in unit tests
    monkeypatch.setattr(k, "check_files_are_found", lambda paths: None)

    out = tmp_path / "summary.tsv"
    k.kraken2_to_taxon_summaries(str(report_min), str(out))

    assert out.exists()

    df = pd.read_csv(out, sep="\t")
    expected_cols = ["taxid"] + k.STD_RANKS + ["read_counts", "percent_of_reads"]
    assert list(df.columns) == expected_cols
    assert (df["taxid"] == 0).any()
    assert df["percent_of_reads"].map(lambda x: isinstance(x, float)).all()


def test_parse_report_ignores_blank_and_handles_U_and_R1(tmp_path):
    text = "\n".join([
        "  ",  # blank
        "5.0 50 50 U 0 unclassified",
        "95.0 950 0 R1 2 Bacteria",
        "20.0 200 200 G 9  GenA",
    ])
    p = tmp_path / "mini.report"
    p.write_text(text)

    df = k.parse_report(str(p))
    assert (df["taxid"] == 0).any()
    gen = df[df["taxid"] == 9].iloc[0]
    assert gen["domain"] == "Bacteria"


def test_preflight_checks_calls_validator(monkeypatch):
    called = {}
    monkeypatch.setattr(k, "check_files_are_found", lambda paths: called.setdefault("ok", paths))
    k.preflight_checks("abc.txt")
    assert "ok" in called and called["ok"] == ["abc.txt"]


def test_parse_report_line_bad_line_uses_report_failure(monkeypatch):
    msgs = {}
    def fake_report_failure(msg):
        raise ValueError(msg)
    monkeypatch.setattr(k, "report_failure", fake_report_failure)

    with pytest.raises(ValueError):
        k.parse_report_line("not enough fields")
