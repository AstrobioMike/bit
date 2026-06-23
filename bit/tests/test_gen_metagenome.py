"""
Tests for the non-trivial *logic* that lives in the gen_metagenome orchestrator
(accession normalization, downloaded-fasta matching, the suppression-screen
backfill loop, the GTDB info summary, and small parsing helpers). The phase
functions themselves are thin wiring over already-tested modules and are not
unit-tested here.
"""
import os
import pandas as pd
import pytest

from bit.modules import gen_metagenome as G


# ─── accession normalization helpers ───────────────────────────────────────

@pytest.mark.parametrize("acc,core", [
    ("GCF_000005845.2", "000005845"),
    ("GCA_000005845.1", "000005845"),     # GCA/GCF twin -> same core
    ("RS_GCF_000005845.2", "000005845"),  # GTDB RefSeq prefix
    ("GB_GCA_000005845.1", "000005845"),  # GTDB GenBank prefix
    ("GCF_000005845", "000005845"),       # no version
])
def test_acc_digits_core(acc, core):
    assert G._acc_digits_core(acc) == core


@pytest.mark.parametrize("path,acc", [
    ("/g/GCF_000005845.2.fasta.gz", "GCF_000005845.2"),
    ("/g/GCA_000009.3.fasta", "GCA_000009.3"),
    ("/g/GCF_111.1.fna.gz", "GCF_111.1"),
    (None, None),
])
def test_accession_from_fasta_name(path, acc):
    assert G._accession_from_fasta_name(path) == acc


# ─── find_downloaded_fasta: version/prefix tolerance (the 500->499 bug) ─────

@pytest.fixture
def genomes_dir(tmp_path):
    for fn in ["GCF_000005845.2.fasta.gz",
               "GCA_000009999.3.fasta.gz",
               "GCF_111111111.1.fasta.gz"]:
        (tmp_path / fn).write_text("")
    return str(tmp_path)


def test_find_downloaded_fasta_exact(genomes_dir):
    p = G.find_downloaded_fasta(genomes_dir, "GCF_000005845.2")
    assert os.path.basename(p) == "GCF_000005845.2.fasta.gz"


def test_find_downloaded_fasta_version_drift(genomes_dir):
    p = G.find_downloaded_fasta(genomes_dir, "GCA_000009999.1")
    assert os.path.basename(p) == "GCA_000009999.3.fasta.gz"


def test_find_downloaded_fasta_gca_gcf_twin(genomes_dir):
    p = G.find_downloaded_fasta(genomes_dir, "GCA_000005845.1")
    assert os.path.basename(p) == "GCF_000005845.2.fasta.gz"


def test_find_downloaded_fasta_genuinely_absent(genomes_dir):
    assert G.find_downloaded_fasta(genomes_dir, "GCA_222222222.1") is None


# ─── suppression-screen backfill loop ──────────────────────────────────────

def _gtdb_row(acc, gen, sp, rep="t"):
    return {"ncbi_genbank_assembly_accession": acc, "accession": "RS_" + acc,
            "domain": "Bacteria", "phylum": "P", "class": "C", "order": "O",
            "family": "F", "genus": gen, "species": sp,
            "gtdb_representative": rep, "ncbi_refseq_category": "",
            "checkm2_completeness": 98.0, "checkm2_contamination": 1.0}


@pytest.fixture
def screen_pool():
    # genus A has 3 species-reps, B has 2, C has 1
    return pd.DataFrame([
        _gtdb_row("GCA_000000001.1", "A", "a1"), _gtdb_row("GCA_000000002.1", "A", "a2"),
        _gtdb_row("GCA_000000003.1", "A", "a3"),
        _gtdb_row("GCA_000000004.1", "B", "b1"), _gtdb_row("GCA_000000005.1", "B", "b2"),
        _gtdb_row("GCA_000000006.1", "C", "c1"),
    ])


def _ai_file(tmp_path, live_accs):
    ai = tmp_path / "ncbi-assembly-info.tsv"
    with open(ai, "w") as fh:
        fh.write("#header\n")
        for a in live_accs:
            fh.write(a + "\t" + "\t".join(["x"] * 22) + "\n")
    return str(ai)


def test_screen_all_live_no_screening(screen_pool, tmp_path):
    live = [f"GCA_00000000{i}.1" for i in range(1, 7)]
    ai = _ai_file(tmp_path, live)
    sel, w = G._select_with_suppression_screen(
        screen_pool, derep_rank="genus", domains=None, num_genomes=3, seed=1,
        user_excluded_groups=set(), assembly_info_path=ai)
    assert len(sel) == 3
    assert set(sel["genus"]) == {"A", "B", "C"}
    assert w == []


def test_screen_backfills_same_group_sibling(screen_pool, tmp_path):
    # genus A's seed-1 pick is suppressed; backfill keeps genus A via a sibling
    live = ["GCA_000000001.1", "GCA_000000002.1",
            "GCA_000000004.1", "GCA_000000005.1", "GCA_000000006.1"]  # 003 suppressed
    ai = _ai_file(tmp_path, live)
    sel, w = G._select_with_suppression_screen(
        screen_pool, derep_rank="genus", domains=None, num_genomes=3, seed=1,
        user_excluded_groups=set(), assembly_info_path=ai)
    assert "GCA_000000003.1" not in set(sel["ncbi_genbank_assembly_accession"])
    assert set(sel["genus"]) == {"A", "B", "C"}
    assert len(sel) == 3
    assert w == []


def test_screen_recovers_when_first_round_all_suppressed(screen_pool, tmp_path):
    # the seed-1 head(1) picks (003/005/006) are ALL suppressed, but live siblings
    # exist for A and B. The loop must not give up on the first fruitless round:
    # it should exclude the dead picks, re-select, and recover 001 + 004.
    live = ["GCA_000000001.1", "GCA_000000004.1"]   # only these two live
    ai = _ai_file(tmp_path, live)
    sel, w = G._select_with_suppression_screen(
        screen_pool, derep_rank="genus", domains=None, num_genomes=3, seed=1,
        user_excluded_groups=set(), assembly_info_path=ai)
    assert set(sel["ncbi_genbank_assembly_accession"]) == {"GCA_000000001.1", "GCA_000000004.1"}
    assert set(sel["genus"]) == {"A", "B"}          # genus C had no live genome
    assert len(sel) == 2                            # shortfall of 1
    assert len(w) == 1 and "returning 2" in w[0].lower()


def test_screen_no_assembly_info_assumes_all_live(screen_pool):
    sel, w = G._select_with_suppression_screen(
        screen_pool, derep_rank="genus", domains=None, num_genomes=3, seed=1,
        user_excluded_groups=set(), assembly_info_path=None)
    assert len(sel) == 3
    assert w == []


def test_screen_terminates_when_nothing_live(screen_pool, tmp_path):
    # no live genomes at all -> returns empty, warns, and does NOT hang
    ai = _ai_file(tmp_path, [])
    sel, w = G._select_with_suppression_screen(
        screen_pool, derep_rank="genus", domains=None, num_genomes=3, seed=1,
        user_excluded_groups=set(), assembly_info_path=ai)
    assert len(sel) == 0
    assert len(w) == 1


# ─── write_gtdb_summary ────────────────────────────────────────────────────

def test_write_gtdb_summary_populates_and_nas(tmp_path):
    gtdb = pd.DataFrame([{
        "ncbi_genbank_assembly_accession": "GCA_000005845.2",
        "domain": "Bacteria", "phylum": "Pseudomonadota", "class": "G",
        "order": "E", "family": "Ent", "genus": "Escherichia", "species": "E. coli",
        "ambiguous_bases": 0, "checkm2_completeness": 99.9, "checkm2_contamination": 0.1,
        "coding_bases": 4_000_000, "coding_density": 88.5, "contig_count": 1,
        "gc_count": 2_300_000, "gc_percentage": 50.8, "genome_size": 4_600_000}])
    merged = pd.DataFrame({"accession": ["GCA_000005845.2", "GCA_999999999.1"]})
    gdir = tmp_path / "genomes"; gdir.mkdir()
    run = type("R", (), {"merged": merged, "gtdb_tab": gtdb, "genomes_dir": str(gdir)})()
    G.write_gtdb_summary(run)

    df = pd.read_csv(run.gtdb_summary_path, sep="\t", keep_default_na=False)
    assert list(df.columns)[0] == "accession"
    ecoli = df[df["accession"] == "GCA_000005845.2"].iloc[0]
    absent = df[df["accession"] == "GCA_999999999.1"].iloc[0]
    assert ecoli["genus"] == "Escherichia"
    assert float(ecoli["checkm2_completeness"]) == 99.9
    assert absent["genus"] in ("", "NA")
    assert absent["checkm2_completeness"] in ("", "NA")


def test_write_gtdb_summary_matches_across_version_and_prefix(tmp_path):
    gtdb = pd.DataFrame([{
        "ncbi_genbank_assembly_accession": "GCA_000005845.2",
        "domain": "Bacteria", "phylum": "P", "class": "C", "order": "O",
        "family": "F", "genus": "Escherichia", "species": "E. coli",
        "genome_size": 4_600_000}])
    merged = pd.DataFrame({"accession": ["GCF_000005845.1"]})
    gdir = tmp_path / "genomes"; gdir.mkdir()
    run = type("R", (), {"merged": merged, "gtdb_tab": gtdb, "genomes_dir": str(gdir)})()
    G.write_gtdb_summary(run)
    df = pd.read_csv(run.gtdb_summary_path, sep="\t", keep_default_na=False)
    assert df.iloc[0]["genus"] == "Escherichia"


# ─── small parsing helpers ─────────────────────────────────────────────────

def test_read_retrieved_date(tmp_path):
    (tmp_path / "date-retrieved.txt").write_text("2026,06,20\n")
    assert G._read_retrieved_date(str(tmp_path)) == "Jun 20, 2026"


def test_read_retrieved_date_missing(tmp_path):
    assert G._read_retrieved_date(str(tmp_path)) == "(date unknown)"
    assert G._read_retrieved_date(None) == "(date unknown)"


def test_read_retrieved_date_malformed(tmp_path):
    (tmp_path / "date-retrieved.txt").write_text("not-a-date\n")
    assert G._read_retrieved_date(str(tmp_path)) == "not-a-date"


@pytest.mark.parametrize("v", [None, float("nan"), 0])
def test_size_or_na_blanks(v):
    assert pd.isna(G._size_or_na(v))


def test_size_or_na_keeps_real():
    assert G._size_or_na(4_600_000) == 4_600_000


@pytest.mark.parametrize("v,expected", [("4.5", 4.5), ("10", 10.0)])
def test_num_or_na_numbers(v, expected):
    assert G._num_or_na(v) == expected


@pytest.mark.parametrize("v", [None, "", "abc"])
def test_num_or_na_blanks(v):
    assert pd.isna(G._num_or_na(v))


# ─── load_user_accessions ──────────────────────────────────────────────────

def test_load_user_accessions_bare_list(tmp_path):
    f = tmp_path / "accs.txt"
    f.write_text("GCF_000005845.2\nGCA_000009999.1\n")
    args = type("A", (), {"accessions": str(f)})()
    udf = G.load_user_accessions(args)
    assert list(udf["accession"]) == ["GCF_000005845.2", "GCA_000009999.1"]
    assert udf["user_supplied"].all()


def test_load_user_accessions_with_pins(tmp_path):
    f = tmp_path / "accs.tsv"
    f.write_text("accession\trel_abundance\tmutation_rate\n"
                 "GCF_000005845.2\t0.7\t0.02\n"
                 "GCA_000009999.1\t0.3\t\n")
    args = type("A", (), {"accessions": str(f)})()
    udf = G.load_user_accessions(args).set_index("accession")
    assert udf.loc["GCF_000005845.2", "pinned_rel_abundance"] == 0.7
    assert udf.loc["GCF_000005845.2", "pinned_mutation_rate"] == 0.02
    assert pd.isna(udf.loc["GCA_000009999.1", "pinned_mutation_rate"])
