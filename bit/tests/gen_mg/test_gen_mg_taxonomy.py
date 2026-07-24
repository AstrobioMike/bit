import pandas as pd # type: ignore
import pytest # type: ignore

from bit.modules.gen_mg import selection as SEL
from bit.modules.gen_mg import taxonomy as TAX
from bit.modules.gen_mg.selection import RANKS


# helpers ────────────────────────────────────────────────────────────────────

def gtdb_meta_row(col0_acc, gca, taxid, genus, species, domain="Bacteria"):
    """ a GTDB metadata row: col-0 'accession' (prefixed RS_/GB_) + GCA twin. """
    return {
        "accession": col0_acc,
        "ncbi_genbank_assembly_accession": gca,
        "ncbi_taxid": taxid,
        "domain": domain, "phylum": "P", "class": "C", "order": "O",
        "family": "F", "genus": genus, "species": species,
    }


@pytest.fixture
def gtdb_tab():
    return pd.DataFrame([
        # RefSeq entry: col-0 RS_GCF, GCA twin -> keyed under both forms
        gtdb_meta_row("RS_GCF_000005845.2", "GCA_000005845.2", "511145",
                      "Escherichia", "Escherichia coli"),
        # GenBank-only entry
        gtdb_meta_row("GB_GCA_000009.1", "GCA_000009.1", "1280",
                      "Staphylococcus", "Staphylococcus aureus"),
    ])


def make_resolver(taxid_to_lineage):
    """ a fake _resolver: taxids -> {taxid: 'd;p;c;o;f;g;s'} (the seam
    fill_ncbi_taxonomy calls). Mirrors what parquet_lineage_resolver returns, but
    from an in-memory map so tests don't need a Parquet. """
    def _resolve(taxids):
        want = {str(t).strip() for t in taxids if str(t).strip()}
        return {t: taxid_to_lineage[t] for t in want if t in taxid_to_lineage}
    return _resolve


# ─── GTDB maps: dual GCA/GCF keying ────────────────────────────────────────

def test_gtdb_taxid_map_keys_both_accession_forms(gtdb_tab):
    m = TAX.gtdb_taxid_map(gtdb_tab)
    # the RefSeq row is reachable by GCF and by its GCA twin
    assert m["GCF_000005845"] == "511145"
    assert m["GCA_000005845"] == "511145"
    assert m["GCA_000009"] == "1280"


def test_gtdb_lineage_map_keys_both_forms(gtdb_tab):
    m = TAX.gtdb_lineage_map(gtdb_tab)
    assert m["GCF_000005845"]["genus"] == "Escherichia"
    assert m["GCA_000005845"]["genus"] == "Escherichia"


def test_gtdb_maps_missing_columns_return_empty():
    assert TAX.gtdb_taxid_map(pd.DataFrame({"foo": [1]})) == {}
    assert TAX.gtdb_lineage_map(pd.DataFrame({"foo": [1]})) == {}


# ─── fill_gtdb_taxonomy ────────────────────────────────────────────────────

def test_fill_gtdb_taxonomy_resolves_user_gcf(gtdb_tab):
    """ a user-supplied GCF accession matches the GCA-keyed GTDB row. """
    u = SEL._normalize_ncbi_rows(pd.DataFrame({"accession": ["GCF_000005845.2"]}),
                                 lineage_lookup=None, user_supplied=True)
    out = TAX.fill_gtdb_taxonomy(u, gtdb_tab)
    assert out["gtdb_genus"].iloc[0] == "Escherichia"


def test_fill_gtdb_taxonomy_leaves_non_gtdb_na(gtdb_tab):
    u = SEL._normalize_ncbi_rows(pd.DataFrame({"accession": ["GCA_999999.1"]}),
                                 lineage_lookup=None, user_supplied=True)
    out = TAX.fill_gtdb_taxonomy(u, gtdb_tab)
    assert out["gtdb_genus"].iloc[0] == "NA"


# ─── assembly-info taxid map (Parquet) ─────────────────────────────────────

# The assembly-info table is now a Parquet file with named columns, so the old
# format concerns (legacy positional vs. slim header, header-row-as-data) no longer
# exist -- those tests are replaced by the behavioural checks below. Columns are read
# by name via pushdown.

import pyarrow as _pa  # type: ignore
import pyarrow.parquet as _pq  # type: ignore

_AI_COLUMNS = [
    "assembly_accession", "taxid", "organism_name", "infraspecific_name",
    "version_status", "assembly_level", "asm_name", "ftp_path",
    "domain", "phylum", "class", "order", "family", "genus", "species",
]


def _write_ai_parquet(path, rows):
    """rows: (accession, taxid) or (accession, taxid, lineage_str). Lineage is a
    ';'-joined 7-rank string; when omitted the rank columns are 'na'."""
    data = {c: [] for c in _AI_COLUMNS}
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for row in rows:
        acc, taxid = row[0], row[1]
        lineage = row[2].split(";") if len(row) > 2 else ["na"] * 7
        data["assembly_accession"].append(acc)
        data["taxid"].append(taxid)
        for c in ("organism_name", "infraspecific_name", "version_status",
                  "assembly_level", "asm_name", "ftp_path"):
            data[c].append("na")
        for r, v in zip(ranks, lineage):
            data[r].append(v)
    _pq.write_table(_pa.table({c: _pa.array(data[c]) for c in _AI_COLUMNS}), str(path))


def test_assembly_info_taxid_map_reads_taxid_by_name(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "4932")])
    m = TAX.assembly_info_taxid_map(["GCA_900.1"], str(ai))
    assert m["GCA_900"] == "4932"


def test_assembly_info_taxid_map_only_keeps_wanted_rows(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "4932"), ("GCA_111.1", "999")])
    m = TAX.assembly_info_taxid_map(["GCA_900.1"], str(ai))
    assert m == {"GCA_900": "4932"}


def test_assembly_info_taxid_map_skips_blank_taxid(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "")])
    assert TAX.assembly_info_taxid_map(["GCA_900.1"], str(ai)) == {}


def test_dead_accession_cores_reads_parquet(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_000005845.2", "562")])
    # a GCF pick whose GCA twin (same digits) is present should count as live
    assert TAX.dead_accession_cores(["GCF_000005845.1"], str(ai)) == set()


def test_assembly_info_taxid_map_missing_file():
    assert TAX.assembly_info_taxid_map(["GCA_1.1"], None) == {}
    assert TAX.assembly_info_taxid_map(["GCA_1.1"], "/no/such/file.parquet") == {}


# ─── tiered gather_taxids ──────────────────────────────────────────────────

def test_gather_taxids_prefers_gtdb_then_assembly_info(gtdb_tab, tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "4932")])     # a euk, not in GTDB

    taxids = TAX.gather_taxids(
        ["GCF_000005845.2", "GCA_900.1"], gtdb_tab, str(ai))
    assert taxids["GCF_000005845"] == "511145"   # from GTDB table
    assert taxids["GCA_900"] == "4932"            # from assembly-info


def test_gather_taxids_no_sources_returns_empty():
    taxids = TAX.gather_taxids(["GCA_xyz.1"], None, None)
    assert taxids == {}


# ─── parquet_lineage_resolver (taxid -> lineage from the NCBI Parquet) ──────

def test_parquet_lineage_resolver_reads_lineage_by_taxid(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [
        ("GCF_1.1", "511145",
         "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
         "Enterobacteriaceae;Escherichia;Escherichia coli"),
    ])
    resolve = TAX.parquet_lineage_resolver(str(ai))
    out = resolve(["511145"])
    assert out["511145"].split(";")[5] == "Escherichia"


def test_parquet_lineage_resolver_empty_and_missing(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCF_1.1", "511145", "Bacteria;P;C;O;F;Escherichia;E coli")])
    resolve = TAX.parquet_lineage_resolver(str(ai))
    assert resolve([]) == {}                 # no taxids
    assert resolve(["999999"]) == {}         # taxid not in table
    # no file -> empty, no crash
    assert TAX.parquet_lineage_resolver("/no/such.parquet")(["511145"]) == {}


def test_parquet_lineage_resolver_preserves_na_ranks(tmp_path):
    """where the Parquet has 'NA' at a rank (NCBI itself has no value), it comes
    through as NA -- same as taxonkit's reformat2 -r NA behavior."""
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_9.1", "4444", "Bacteria;P;C;O;F;NA;NA")])
    resolve = TAX.parquet_lineage_resolver(str(ai))
    parts = resolve(["4444"])["4444"].split(";")
    assert parts[5] == "NA" and parts[6] == "NA"


# ─── fill_ncbi_taxonomy + resolve_all (full, injected resolver) ────────────

@pytest.fixture
def merged_mixed(gtdb_tab):
    """ a GTDB-selected genome + a user euk accession (not in GTDB). """
    g = SEL._normalize_gtdb_rows(pd.DataFrame([{
        "ncbi_genbank_assembly_accession": "GCA_000005845.2",
        "domain": "Bacteria", "phylum": "P", "class": "C", "order": "O",
        "family": "F", "genus": "Escherichia", "species": "Escherichia coli"}]))
    u = SEL._normalize_ncbi_rows(pd.DataFrame({"accession": ["GCA_900.1"]}),
                                 lineage_lookup=None, user_supplied=True)
    merged, _ = SEL.merge_sources(g, u)
    return merged


def test_resolve_all_fills_both_taxonomies(gtdb_tab, merged_mixed, tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "4932")])

    resolver = make_resolver({
        "511145": "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                  "Enterobacteriaceae;Escherichia;Escherichia coli",
        "4932": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"})

    out = TAX.resolve_all(merged_mixed, gtdb_tab, str(ai), _resolver=resolver)
    r = out.set_index("accession")

    # GTDB genome: gtdb already there, ncbi resolved via GTDB-table taxid
    assert r.loc["GCA_000005845.2", "gtdb_genus"] == "Escherichia"
    assert r.loc["GCA_000005845.2", "ncbi_genus"] == "Escherichia"
    # euk user genome: gtdb stays NA, ncbi resolved via assembly-info taxid
    assert r.loc["GCA_900.1", "gtdb_genus"] == "NA"
    assert r.loc["GCA_900.1", "ncbi_genus"] == "Saccharomyces"


def test_resolve_all_uses_parquet_as_terminal_taxid_source(gtdb_tab, merged_mixed, tmp_path):
    """
    The euk user accession isn't in GTDB, so its taxid must come from the NCBI
    Parquet (the terminal source now that the datasets CLI fallback is gone). If the
    Parquet lookup works, ncbi ranks fill; there is no other fallback.
    """
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [("GCA_900.1", "4932")])

    resolver = make_resolver({
        "511145": "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                  "Enterobacteriaceae;Escherichia;Escherichia coli",
        "4932": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"})
    out = TAX.resolve_all(merged_mixed, gtdb_tab, str(ai), _resolver=resolver)
    r = out.set_index("accession")
    # euk genome resolved purely via the Parquet taxid -> lineage
    assert r.loc["GCA_900.1", "ncbi_genus"] == "Saccharomyces"


def test_datasets_fallback_is_gone(gtdb_tab, merged_mixed, tmp_path):
    """Regression guard: the datasets CLI fallback has been removed entirely, so
    neither the function nor the plumbing kwarg should exist."""
    assert not hasattr(TAX, "datasets_taxid_map")
    import inspect
    assert "use_datasets_fallback" not in inspect.signature(TAX.gather_taxids).parameters
    assert "use_datasets_fallback" not in inspect.signature(TAX.resolve_all).parameters


def test_taxonkit_is_gone_from_gen_mg(gtdb_tab, merged_mixed, tmp_path):
    """Regression guard: gen-mg no longer resolves lineages via taxonkit; the
    default resolver reads them from the NCBI Parquet."""
    assert not hasattr(TAX, "resolve_lineages")
    assert hasattr(TAX, "parquet_lineage_resolver")


def test_resolve_all_default_resolver_reads_from_parquet(gtdb_tab, merged_mixed, tmp_path):
    """
    End-to-end with the REAL default resolver (no injected fake): the euk user
    accession's NCBI ranks must be filled from lineage columns in the NCBI Parquet.
    """
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai_parquet(ai, [
        # the euk user genome: taxid 4932, lineage present in the Parquet
        ("GCA_900.1", "4932",
         "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
         "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"),
        # the GTDB E. coli genome (taxid 511145 from the gtdb_tab fixture)
        ("GCA_000005845.2", "511145",
         "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
         "Enterobacteriaceae;Escherichia;Escherichia coli"),
    ])
    out = TAX.resolve_all(merged_mixed, gtdb_tab, str(ai))   # default resolver
    r = out.set_index("accession")
    assert r.loc["GCA_900.1", "ncbi_genus"] == "Saccharomyces"
    assert r.loc["GCA_000005845.2", "ncbi_genus"] == "Escherichia"


# ─── dead_accession_cores (suppression screen) ─────────────────────────────

def _write_ai(path, accessions):
    _write_ai_parquet(path, [(a, "999") for a in accessions])


def test_dead_accession_cores_detects_absence(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai(ai, ["GCA_000000003.1", "GCA_000000004.1"])   # 002 absent (suppressed)
    dead = TAX.dead_accession_cores(
        ["GCA_000000002.1", "GCA_000000003.1", "GCA_000000004.1"], str(ai))
    assert dead == {"000000002"}


def test_dead_accession_cores_matches_gca_gcf_twin(tmp_path):
    # table has the GenBank (GCA) row; a GCF query must still register present
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai(ai, ["GCA_000005845.2"])
    assert TAX.dead_accession_cores(["GCF_000005845.2"], str(ai)) == set()


def test_dead_accession_cores_version_tolerant(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai(ai, ["GCA_000005845.3"])           # newer version in table
    # older version queried -- same numeric core, so still live
    assert TAX.dead_accession_cores(["GCA_000005845.1"], str(ai)) == set()


def test_dead_accession_cores_missing_file():
    assert TAX.dead_accession_cores(["GCA_000000001.1"], None) == set()
    assert TAX.dead_accession_cores(["GCA_000000001.1"], "/no/such.parquet") == set()


def test_dead_accession_cores_empty_candidates(tmp_path):
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai(ai, ["GCA_000000003.1"])
    assert TAX.dead_accession_cores([], str(ai)) == set()


def test_dead_accession_cores_spans_batches(tmp_path, monkeypatch):
    """
    A live accession must be found no matter which batch it lands in -- the
    screen ORs each batch into a running mask rather than checking one batch.
    """
    monkeypatch.setattr(TAX, "NCBI_SCREEN_BATCH_ROWS", 2)
    ai = tmp_path / "ncbi-data.parquet"
    _write_ai(ai, [f"GCA_00000000{i}.1" for i in range(1, 8)])   # 7 rows -> 4 batches
    # first, last, and a middle accession are live; 009 was never in the table
    dead = TAX.dead_accession_cores(
        ["GCA_000000001.1", "GCA_000000004.1", "GCA_000000007.1", "GCA_000000009.1"],
        str(ai))
    assert dead == {"000000009"}
