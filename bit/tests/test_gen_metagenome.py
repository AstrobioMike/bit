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
from types import SimpleNamespace

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


# ─── per-read truth build (both taxonomies, in-process) ────────────────────

def test_build_truth_for_taxonomy_per_read_both(tmp_path):
    """ each taxonomy's per-read table is built correctly in-process (the polars
    streaming path that replaced the old multiprocessing parallelization). """
    from bit.modules.gen_mg import truth as TRU
    src = tmp_path / "mg-read-sources.tsv.gz"
    pd.DataFrame({
        "read_id": [f"r{i}" for i in range(1500)],
        "source_fasta": [f"/g/GCF_00000000{i%3}.1.fasta.gz" for i in range(1500)],
        "contig": "c0", "start": range(1500), "end": range(1500),
        "strand": "+", "wrapped": "false",
    }).to_csv(src, sep="\t", index=False, compression="gzip")
    accs = ["GCF_000000000.1", "GCF_000000001.1", "GCF_000000002.1"]
    merged = pd.DataFrame({
        "accession": accs, "assigned_reads": [500, 500, 500],
        "assigned_coverage": [1.0, 1.0, 1.0],
        "assigned_rel_abundance": [0.34, 0.33, 0.33],
        "used_genome_size": [1000, 1000, 1000],
        "downloaded_genome_size": [1000, 1000, 1000]})
    for r in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        merged["gtdb_" + r] = f"g_{r}"
        merged["ncbi_" + r] = f"n_{r}"
    mut = {a: {"rate": 0.0} for a in accs}
    f2a = {f"/g/{a}.fasta.gz": a for a in accs}

    gt_root = tmp_path / "ground-truth"; gt_root.mkdir()
    written = {t: TRU.build_truth_for_taxonomy(
        t, merged, mut, str(gt_root), read_sources_tsv=str(src),
        fasta_to_accession=f2a, per_read=True) for t in TRU.TAXONOMIES}

    assert set(written) == set(TRU.TAXONOMIES)
    for t in TRU.TAXONOMIES:
        df = pd.read_csv(written[t]["per_read"], sep="\t", dtype=str, keep_default_na=False)
        assert len(df) == 1500
        assert (df["accession"] != "NA").all()
        assert df["genus"].iloc[0] == ("g_genus" if t == "gtdb" else "n_genus")


# ─── _load_gtdb_table reads only used columns (usecols intersect) ───────────

def test_load_gtdb_table_usecols_intersect(tmp_path, monkeypatch):
    """ only the used columns are read; extras are dropped and columns absent in
    older GTDB releases (coding_density, checkm2_*) don't raise. """
    gd = tmp_path / "gtdb"; gd.mkdir()
    cols = ["accession", "ncbi_genbank_assembly_accession", "ncbi_taxid",
            "gtdb_representative", "ncbi_refseq_category",
            "domain", "phylum", "class", "order", "family", "genus", "species",
            "checkm_completeness", "checkm_contamination",  # old style, no checkm2_*
            "genome_size", "contig_count", "gc_count", "gc_percentage",
            "ambiguous_bases", "coding_bases",
            "filler1", "filler2"]  # extras to drop; coding_density absent
    row = ["RS_GCF_000000001.1", "GCA_000000001.1", "123", "t", "reference genome",
           "Bacteria", "p", "c", "o", "f", "g", "s", "98.0", "1.0",
           "3000000", "1", "1500000", "50.0", "0", "2700000", "junk", "junk"]
    (gd / "GTDB-arc-and-bac-metadata.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(row) + "\n")
    (gd / "GTDB-version-info.txt").write_text("r220\n2024-04-24\n")

    monkeypatch.setattr(G, "get_gtdb_data", lambda quiet=True: str(gd))
    run = SimpleNamespace(log_file=None)
    tab = G._load_gtdb_table(SimpleNamespace(), run)

    assert not any(c.startswith("filler") for c in tab.columns)   # extras dropped
    assert "coding_density" not in tab.columns                    # absent, no error
    assert "checkm_completeness" in tab.columns                   # old-style kept
    assert "checkm2_completeness" not in tab.columns
    # every present used-column made it in
    for c in ["accession", "ncbi_genbank_assembly_accession", "domain", "species",
              "genome_size", "ncbi_taxid"]:
        assert c in tab.columns


# ─── parallel per-read truth build (processes) ─────────────────────────────

def _parallel_truth_inputs(tmp_path):
    import gzip
    src = tmp_path / "rs.tsv.gz"
    n = 1500
    with gzip.open(src, "wt") as fh:
        fh.write("read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n")
        for i in range(n):
            fh.write(f"r{i}\t/g/GCF_{i%3:09d}.fasta.gz\tc0\t0\t150\t+\tfalse\n")
    accs = [f"GCF_{i:09d}" for i in range(3)]
    merged = pd.DataFrame({
        "accession": accs, "assigned_reads": [n // 3] * 3,
        "assigned_coverage": [1.0] * 3, "assigned_rel_abundance": [1/3] * 3,
        "used_genome_size": [1000] * 3, "downloaded_genome_size": [1000] * 3})
    for r in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        merged["gtdb_" + r] = f"g_{r}"
        merged["ncbi_" + r] = f"n_{r}"
    run = SimpleNamespace(merged=merged,
                          mutation_results={a: {"rate": 0.0} for a in accs},
                          read_sources_tsv=str(src))
    f2a = {f"/g/{a}.fasta.gz": a for a in accs}
    return run, f2a, n


class _FakeBar:
    """ minimal stand-in for a tqdm bar (records updates, no terminal I/O). """
    def __init__(self): self.n = 0
    def update(self, k=1): self.n += k


def test_build_truth_processes(tmp_path):
    from bit.modules.gen_mg import truth as TRU
    run, f2a, n = _parallel_truth_inputs(tmp_path)
    gt = tmp_path / "gt"; gt.mkdir()
    written = G._build_truth_processes(run, str(gt), f2a, 500, _FakeBar())
    assert set(written) == set(TRU.TAXONOMIES)
    for t in TRU.TAXONOMIES:
        df = pd.read_csv(written[t]["per_read"], sep="\t", dtype=str, keep_default_na=False)
        assert len(df) == n
        assert "wrapped" not in df.columns
        assert df["genus"].iloc[0] == ("g_genus" if t == "gtdb" else "n_genus")


def test_processes_match_serial(tmp_path):
    from bit.modules.gen_mg import truth as TRU
    run, f2a, n = _parallel_truth_inputs(tmp_path)
    s_root = tmp_path / "s"; s_root.mkdir()
    ser = {t: TRU.build_truth_for_taxonomy(
        t, run.merged, run.mutation_results, str(s_root),
        read_sources_tsv=run.read_sources_tsv, fasta_to_accession=f2a,
        per_read=True, chunksize=500) for t in TRU.TAXONOMIES}
    p_root = tmp_path / "p"; p_root.mkdir()
    par = G._build_truth_processes(run, str(p_root), f2a, 500, _FakeBar())
    for t in TRU.TAXONOMIES:
        a = pd.read_csv(ser[t]["per_read"], sep="\t", dtype=str, keep_default_na=False).sort_values("read_id").reset_index(drop=True)
        b = pd.read_csv(par[t]["per_read"], sep="\t", dtype=str, keep_default_na=False).sort_values("read_id").reset_index(drop=True)
        assert a.equals(b)


# ─── input-TSV-driven mode resolution ──────────────────────────────────────

import bit.cli.gen_metagenome as CLI


def _acc_file(tmp_path, content):
    p = tmp_path / "accs.tsv"
    p.write_text(content)
    return str(p)


def _mode_args(tmp_path, content, **over):
    base = dict(accessions=_acc_file(tmp_path, content), abundance_mode=None,
                mutation_mode=None, num_genomes=None, mutation_rate=None)
    base.update(over)
    return SimpleNamespace(**base)


def test_resolve_coverage_column_autoswitch(tmp_path):
    a = _mode_args(tmp_path, "accession\tcoverage\nGCF_1\t50\nGCF_2\t30\n")
    CLI.resolve_input_driven_modes(a)
    assert a.abundance_mode == "coverage"


def test_resolve_coverage_column_conflicts_with_relative(tmp_path):
    a = _mode_args(tmp_path, "accession\tcoverage\nGCF_1\t50\n", abundance_mode="relative")
    with pytest.raises(SystemExit):
        CLI.resolve_input_driven_modes(a)


def test_resolve_rel_abundance_column_relative(tmp_path):
    a = _mode_args(tmp_path, "accession\trel_abundance\nGCF_1\t0.3\nGCF_2\t0.2\n")
    CLI.resolve_input_driven_modes(a)
    assert a.abundance_mode == "relative"


def test_resolve_both_columns_default_relative(tmp_path):
    a = _mode_args(tmp_path, "accession\trel_abundance\tcoverage\nGCF_1\t0.3\t50\n")
    CLI.resolve_input_driven_modes(a)
    assert a.abundance_mode == "relative"


def test_resolve_both_columns_explicit_coverage_ok(tmp_path):
    a = _mode_args(tmp_path, "accession\trel_abundance\tcoverage\nGCF_1\t0.3\t50\n",
                   abundance_mode="coverage")
    CLI.resolve_input_driven_modes(a)
    assert a.abundance_mode == "coverage"


def test_resolve_pins_sum_ge_one_with_num_genomes_exits(tmp_path):
    a = _mode_args(tmp_path, "accession\trel_abundance\nGCF_1\t0.6\nGCF_2\t0.5\n",
                   num_genomes=10)
    with pytest.raises(SystemExit):
        CLI.resolve_input_driven_modes(a)


def test_resolve_mutation_uniform_column(tmp_path):
    a = _mode_args(tmp_path, "accession\tmutation_rate\nGCF_1\t0.02\nGCF_2\t0.02\n")
    CLI.resolve_input_driven_modes(a)
    assert a.mutation_mode == "uniform"
    assert a.mutation_rate == 0.02


def test_resolve_mutation_varied_column_autoswitch(tmp_path):
    a = _mode_args(tmp_path, "accession\tmutation_rate\nGCF_1\t0.02\nGCF_2\t0.05\n")
    CLI.resolve_input_driven_modes(a)
    assert a.mutation_mode == "varied"


def test_resolve_mutation_all_empty_column_is_off(tmp_path):
    # a mutation_rate column with no values is treated as absent (no mutation)
    a = _mode_args(tmp_path, "accession\tmutation_rate\nGCF_1\t\nGCF_2\t\n")
    CLI.resolve_input_driven_modes(a)
    assert a.mutation_mode == "off"


def test_resolve_mutation_column_conflicts_with_off(tmp_path):
    a = _mode_args(tmp_path, "accession\tmutation_rate\nGCF_1\t0.02\n", mutation_mode="off")
    with pytest.raises(SystemExit):
        CLI.resolve_input_driven_modes(a)


def test_resolve_bare_list_defaults(tmp_path):
    a = _mode_args(tmp_path, "GCF_1\nGCF_2\n")
    CLI.resolve_input_driven_modes(a)
    assert a.abundance_mode == "relative"
    assert a.mutation_mode == "off"


# ─── reproducible-input.tsv ─────────────────────────────────────────────────

def test_write_reproducible_input_relative_with_mutation(tmp_path):
    merged = pd.DataFrame({"accession": ["GCF_1", "GCF_2"],
                           "assigned_rel_abundance": [0.6, 0.4],
                           "assigned_coverage": [50.0, 30.0]})
    run = SimpleNamespace(merged=merged, out_dir=str(tmp_path),
                          mutation_results={"GCF_1": {"rate": 0.02}, "GCF_2": {"rate": 0.02}})
    G.write_reproducible_input(
        SimpleNamespace(abundance_mode="relative", mutation_mode="uniform"), run)
    df = pd.read_csv(tmp_path / "reproducible-input.tsv", sep="\t")
    # both abundance columns are emitted, plus mutation_rate
    assert list(df.columns) == ["accession", "rel_abundance", "coverage", "mutation_rate"]


def test_write_reproducible_input_no_mutation(tmp_path):
    merged = pd.DataFrame({"accession": ["GCF_1", "GCF_2"],
                           "assigned_rel_abundance": [0.6, 0.4],
                           "assigned_coverage": [50.0, 30.0]})
    run = SimpleNamespace(merged=merged, out_dir=str(tmp_path), mutation_results={})
    G.write_reproducible_input(
        SimpleNamespace(abundance_mode="coverage", mutation_mode="off"), run)
    df = pd.read_csv(tmp_path / "reproducible-input.tsv", sep="\t")
    # both abundance columns, no mutation_rate when mutation was off
    assert list(df.columns) == ["accession", "rel_abundance", "coverage"]
    assert "mutation_rate" not in df.columns


def test_load_user_accessions_skips_bare_header(tmp_path):
    p = tmp_path / "accs.txt"
    p.write_text("GCF_000005845.2\nGCF_000009999.1\n")
    udf = G.load_user_accessions(SimpleNamespace(accessions=str(p)))
    assert list(udf["accession"]) == ["GCF_000005845.2", "GCF_000009999.1"]
