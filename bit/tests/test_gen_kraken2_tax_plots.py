import io
from pathlib import Path
import pandas as pd
import pytest
from bit.modules import gen_kraken2_tax_plots as g


@pytest.fixture(autouse=True)
def _use_non_gui_matplotlib(monkeypatch):
    import matplotlib
    matplotlib.use("Agg", force=True)
    calls = []
    def fake_savefig(path, *a, **k):
        calls.append(Path(path))
        # creating a tiny file so the path exists for assertions
        Path(path).write_bytes(b"ok")
    monkeypatch.setattr(g.plt, "savefig", fake_savefig)
    return calls


@pytest.fixture
def small_report(tmp_path):
    """
    Make a tiny Kraken-style report with:
      - U (unclassified)
      - R (root)
      - D (Bacteria)
      - G (two genera)
      - S (one species)
    Note: taxon_name is intentionally indented to mimic real reports; the parser must strip it.
    """
    text = "\n".join([
        # percent  clade  this_rank  taxid  taxon
        "10.0\t100\t100\tU\t0\tunclassified",
        "90.0\t900\t50\tR\t1\troot",
        "80.0\t800\t0\tD\t2\tBacteria",
        "30.0\t300\t300\tG\t123\t  GenX",   # leading spaces should be stripped
        "20.0\t200\t200\tG\t124\tGenY  ",   # trailing spaces should be stripped
        "15.0\t150\t150\tS\t200\tGenX speciesA",
    ]) + "\n"
    p = tmp_path / "kraken.report"
    p.write_text(text)
    return p


def test_read_in_kraken_report_strips_whitespace(small_report):
    df = g.read_in_kraken_report(str(small_report))
    assert "taxon_name" in df.columns
    assert "GenX" in df.loc[df["rank_code"] == "G", "taxon_name"].tolist()
    assert "GenY" in df.loc[df["rank_code"] == "G", "taxon_name"].tolist()


def test_get_stats_values(small_report):
    df = g.read_in_kraken_report(str(small_report))
    total_reads, unclassified_pct, stuck_pct = g.get_stats(df)
    # From the fixture:
    # total_reads = clade(U) + clade(root) = 100 + 900 = 1000
    assert total_reads == 1000
    # unclassified percent = 10.0
    assert pytest.approx(unclassified_pct, rel=1e-6) == 10.0
    # stuck_at_root = this_level(root) / total_reads * 100 = 50/1000*100 = 5.0
    # function returns as a trimmed string ("5" not "5.0")
    assert stuck_pct == "5"


def test_reduce_kraken_report_filters_columns_and_ranks(small_report):
    df = g.read_in_kraken_report(str(small_report))
    reduced = g.reduce_kraken_report(df)
    # only these columns remain:
    assert list(reduced.columns) == ["rank_code", "taxon_name", "percent", "clade_reads"]
    # only target ranks (U, R, R1, and the rank letters in RANK_DICT)
    allowed = set(g.RANK_DICT.keys()) | {"U", "R", "R1"}
    assert set(reduced["rank_code"].unique()).issubset(allowed)


def test_add_domain_letter_prefixes_after_domain_rows(small_report):
    df = g.read_in_kraken_report(str(small_report))
    reduced = g.reduce_kraken_report(df)
    out = g.add_domain_letter_to_taxon_name(reduced.copy())

    # domain row itself should NOT be prefixed
    row_D = out[out["rank_code"] == "D"].iloc[0]
    assert row_D["taxon_name"] == "Bacteria"

    # rows after domain should be prefixed with (B) in this case
    gen_rows = out[out["rank_code"] == "G"]["taxon_name"].tolist()
    assert any(name.startswith("(B) ") for name in gen_rows)
    assert all("(B) " in name for name in gen_rows)


def test_plot_barplot_and_save_other_is_last(tmp_path, _use_non_gui_matplotlib, monkeypatch):
    df = pd.DataFrame({
        "taxon_name": ["Alpha", "Other", "Beta"],
        "percent": [60.0, 10.0, 30.0],
    })

    labels_rendered = {}
    def fake_set_xticklabels(labels, *a, **k):
        labels_rendered["labels"] = labels

    class _FakeAxes:
        def bar(self, *a, **k): pass
        def set_xticks(self, *a, **k): pass
        set_xticklabels = staticmethod(fake_set_xticklabels)
        def set_title(self, *a, **k): pass
        def set_ylabel(self, *a, **k): pass
        def text(self, *a, **k): pass

    def fake_subplots(*a, **k):
        class _FakeFig:
            def tight_layout(self):
                pass
        return _FakeFig(), _FakeAxes()

    monkeypatch.setattr(g.plt, "subplots", fake_subplots)

    out_prefix = str(tmp_path / "out")
    g.plot_barplot_and_save(
        df_to_plot=df,
        rank_name="Genus",
        out_prefix=out_prefix,
        total_reads=1000,
        unclassified_percent=10.0,
        percent_reads_stuck_at_root="5",
        no_annots=True,
    )

    assert "labels" in labels_rendered
    joined = "\t".join(labels_rendered["labels"])
    assert joined.strip().endswith("Other")


def test_gen_kraken2_tax_plots_topN_creates_pngs(tmp_path, small_report, monkeypatch):
    monkeypatch.setattr(g, "check_files_are_found", lambda paths: None)

    out_prefix = str(tmp_path / "sample")
    g.gen_kraken2_tax_plots(
        input_kraken=str(small_report),
        output_prefix=out_prefix,
        max_taxa=5,
        min_percent=0.0,
        no_annots=True,
    )

    # expecting at least domain and genus plots to be produced
    expected_some = [
        Path(f"{out_prefix}-domain-barplot.png"),
        Path(f"{out_prefix}-genus-barplot.png"),
    ]
    assert any(p.exists() for p in expected_some), "Expected at least one of the top-N PNGs to exist"


def test_gen_kraken2_tax_plots_threshold_creates_other_and_skips_when_empty(tmp_path, small_report, monkeypatch, capsys):
    monkeypatch.setattr(g, "check_files_are_found", lambda paths: None)

    out_prefix = str(tmp_path / "sample2")
    # setting min_percent threshold to 25% so at Genus rank: GenX=30 meets threshold, GenY=20 falls into 'Other'
    g.gen_kraken2_tax_plots(
        input_kraken=str(small_report),
        output_prefix=out_prefix,
        max_taxa=3,
        min_percent=25.0,
        no_annots=True,
    )

    # should produce at least a genus plot
    genus_png = Path(f"{out_prefix}-genus-barplot.png")
    assert genus_png.exists()
