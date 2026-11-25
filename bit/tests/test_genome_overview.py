import pytest
from bit.modules import genome_overview as m

STATS_HEADER = ["Assembly", "Total contigs", "Total length", "Ambiguous characters", "GC content", "Maximum contig length",
                "Minimum contig length", "N50", "N75", "N90", "L50", "L75", "L90", "Num. contigs >= 100",
                "Num. contigs >= 500", "Num. contigs >= 1000", "Num. contigs >= 5000", "Num. contigs >= 10000",
                "Num. contigs >= 50000", "Num. contigs >= 100000"]


def write_tsv(path, rows):
    path.write_text("".join("\t".join(str(c) for c in row) + "\n" for row in rows))
    return path


@pytest.fixture
def stats(tmp_path):
    rows = [STATS_HEADER,
            ["g1", 46, 2095693, 0, 41.87, 134644, 1212, 83489, 42329, 23490, 10, 19, 28, 46, 46, 46, 39, 39, 16, 7],
            ["g2", 109, 4199784, 291, 63.45, 402253, 212, 105665, 75212, 24362, 10, 22, 37, 109, 88, 85, 59, 52, 26, 11],
            ["g3", 2, 1000, 0, 50.0, 600, 400, 600, 400, 400, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0]]
    return write_tsv(tmp_path / "stats.tsv", rows)


@pytest.fixture
def checkm2(tmp_path):
    return write_tsv(tmp_path / "checkm2.tsv", [["Name", "Completeness", "Contamination", "Completeness_Model_Used"],
                                                ["g2", 98.88, 2.24, "x"], ["g1", 96.06, 0.06, "x"]])


@pytest.mark.parametrize("lineage, expected", [
    ("d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus",
     ["Bacteria", "Bacillota", "Bacilli", "Staphylococcales", "Staphylococcaceae", "Staphylococcus", "Staphylococcus aureus"]),
    ("d__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__", ["Bacteria", "Bacillota", "Bacilli", "NA", "NA", "NA", "NA"]),
    ("Unclassified Bacteria", ["Bacteria", "NA", "NA", "NA", "NA", "NA", "NA"]),
    ("Unclassified", ["NA"] * 7),
])
def test_split_gtdb_lineage(lineage, expected):
    assert m.split_gtdb_lineage(lineage) == expected


def test_overview_without_taxonomy_keeps_the_stats_order(tmp_path, stats, checkm2):
    out = tmp_path / "out.tsv"
    m.write_overview(stats, out, checkm2_path=checkm2)
    lines = [l.split("\t") for l in out.read_text().splitlines()]
    assert lines[0] == ["Assembly"] + m.STATS_COLUMNS + m.QUALITY_COLUMNS
    assert [l[0] for l in lines[1:]] == ["g1", "g2", "g3"]
    assert lines[1] == ["g1", "46", "2095693", "41.87", "134644", "83489", "10", "39", "16", "7", "96.06", "0.06"]
    # a genome checkm2 gave nothing for
    assert lines[3][-2:] == ["NA", "NA"]


def test_overview_with_gtdb_taxonomy_is_ordered_by_it(tmp_path, stats, checkm2):
    bac = write_tsv(tmp_path / "gtdbtk.bac120.summary.tsv", [["user_genome", "classification", "other"],
                    ["g1", "d__Bacteria;p__Bacillota;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella_A;s__", "x"]])
    ar = write_tsv(tmp_path / "gtdbtk.ar53.summary.tsv", [["user_genome", "classification", "other"],
                   ["g2", "d__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;f__Halococcaceae;g__Halococcus;s__Halococcus salifodinae", "x"]])
    out = tmp_path / "out.tsv"
    m.write_overview(stats, out, checkm2_path=checkm2, gtdbtk_paths=[bac, ar])
    lines = [l.split("\t") for l in out.read_text().splitlines()]
    assert lines[0][-7:] == m.RANKS
    assert [l[0] for l in lines[1:]] == ["g2", "g1", "g3"]
    assert lines[2][-7:] == ["Bacteria", "Bacillota", "Negativicutes", "Veillonellales", "Veillonellaceae", "Veillonella_A", "NA"]
    assert lines[3][-7:] == ["NA"] * 7
    assert len(set(len(l) for l in lines)) == 1


def test_overview_from_eukcc_and_a_rank_table(tmp_path, stats):
    eukcc = write_tsv(tmp_path / "eukcc.tsv", [["Assembly", "Est. Comp.", "Est. Redund."], ["g1", 90.1, 1.2]])
    ranks = write_tsv(tmp_path / "tax.tsv", [["Assembly"] + m.RANKS,
                      ["g1", "Eukaryota", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "Saccharomycetaceae", "Saccharomyces", "NA"]])
    out = tmp_path / "out.tsv"
    m.write_overview(stats, out, eukcc_path=eukcc, rank_table_path=ranks)
    lines = [l.split("\t") for l in out.read_text().splitlines()]
    g1 = next(l for l in lines if l[0] == "g1")
    assert g1[-9:-7] == ["90.1", "1.2"] and g1[-7] == "Eukaryota" and g1[-1] == "NA"


def test_a_note_instead_of_stats_is_passed_through(tmp_path):
    note = tmp_path / "stats.tsv"
    note.write_text("There were no bins recovered.\n")
    out = tmp_path / "out.tsv"
    m.write_overview(note, out, checkm2_path=note)
    assert out.read_text() == "There were no bins recovered.\n"


def test_missing_stats_column_is_an_error(tmp_path):
    bad = write_tsv(tmp_path / "stats.tsv", [["Assembly", "Total contigs"], ["g1", 1]])
    with pytest.raises(ValueError):
        m.write_overview(bad, tmp_path / "out.tsv")
