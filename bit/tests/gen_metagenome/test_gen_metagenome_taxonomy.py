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


def fake_runner(taxid_to_lineage):
    """ build a get_lineage_from_taxids-shaped runner from a taxid->lineage map. """
    def _runner(tin, tout):
        with open(tin) as fh:
            taxids = [l.strip() for l in fh if l.strip()]
        with open(tout, "w") as fh:
            fh.write("taxid\t" + "\t".join(RANKS) + "\n")
            for t in taxids:
                if t in taxid_to_lineage:
                    fh.write(t + "\t" + taxid_to_lineage[t].replace(";", "\t") + "\n")
    return _runner


def make_resolver(taxid_to_lineage):
    runner = fake_runner(taxid_to_lineage)
    return lambda taxids: TAX.resolve_lineages(taxids, _runner=runner)


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


# ─── assembly-info taxid map ───────────────────────────────────────────────

def test_assembly_info_taxid_map_reads_field_5(tmp_path):
    ai = tmp_path / "ncbi-assembly-info.tsv"
    row = ["na"] * 23
    row[0] = "GCA_900.1"
    row[5] = "4932"
    ai.write_text("\t".join(row) + "\n")
    m = TAX.assembly_info_taxid_map(["GCA_900.1"], str(ai))
    assert m["GCA_900"] == "4932"


def test_assembly_info_taxid_map_missing_file():
    assert TAX.assembly_info_taxid_map(["GCA_1.1"], None) == {}
    assert TAX.assembly_info_taxid_map(["GCA_1.1"], "/no/such/file.tsv") == {}


# ─── tiered gather_taxids ──────────────────────────────────────────────────

def test_gather_taxids_prefers_gtdb_then_assembly_info(gtdb_tab, tmp_path):
    ai = tmp_path / "ncbi-assembly-info.tsv"
    row = ["na"] * 23
    row[0] = "GCA_900.1"     # a euk, not in GTDB
    row[5] = "4932"
    ai.write_text("\t".join(row) + "\n")

    taxids = TAX.gather_taxids(
        ["GCF_000005845.2", "GCA_900.1"], gtdb_tab, str(ai),
        use_datasets_fallback=False)
    assert taxids["GCF_000005845"] == "511145"   # from GTDB table
    assert taxids["GCA_900"] == "4932"            # from assembly-info


def test_gather_taxids_no_sources_returns_empty():
    taxids = TAX.gather_taxids(["GCA_xyz.1"], None, None, use_datasets_fallback=False)
    assert taxids == {}


# ─── resolve_lineages (mocked taxonkit) ────────────────────────────────────

def test_resolve_lineages_parses_runner_output():
    runner = fake_runner({
        "511145": "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                  "Enterobacteriaceae;Escherichia;Escherichia coli"})
    lineages = TAX.resolve_lineages(["511145"], _runner=runner)
    assert lineages["511145"].split(";")[5] == "Escherichia"


def test_resolve_lineages_empty_input():
    assert TAX.resolve_lineages([]) == {}


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
    ai = tmp_path / "ncbi-assembly-info.tsv"
    row = ["na"] * 23
    row[0] = "GCA_900.1"
    row[5] = "4932"
    ai.write_text("\t".join(row) + "\n")

    resolver = make_resolver({
        "511145": "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                  "Enterobacteriaceae;Escherichia;Escherichia coli",
        "4932": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"})

    out = TAX.resolve_all(merged_mixed, gtdb_tab, str(ai),
                          use_datasets_fallback=False, _resolver=resolver)
    r = out.set_index("accession")

    # GTDB genome: gtdb already there, ncbi resolved via GTDB-table taxid
    assert r.loc["GCA_000005845.2", "gtdb_genus"] == "Escherichia"
    assert r.loc["GCA_000005845.2", "ncbi_genus"] == "Escherichia"
    # euk user genome: gtdb stays NA, ncbi resolved via assembly-info taxid
    assert r.loc["GCA_900.1", "gtdb_genus"] == "NA"
    assert r.loc["GCA_900.1", "ncbi_genus"] == "Saccharomyces"


def test_resolve_all_no_datasets_calls(gtdb_tab, merged_mixed, tmp_path, monkeypatch):
    """ in the normal path, the datasets CLI fallback is never invoked. """
    ai = tmp_path / "ncbi-assembly-info.tsv"
    row = ["na"] * 23
    row[0] = "GCA_900.1"
    row[5] = "4932"
    ai.write_text("\t".join(row) + "\n")

    called = {"datasets": False}
    def _boom(accessions):
        called["datasets"] = True
        return {}
    monkeypatch.setattr(TAX, "datasets_taxid_map", _boom)

    resolver = make_resolver({
        "511145": "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                  "Enterobacteriaceae;Escherichia;Escherichia coli",
        "4932": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"})
    TAX.resolve_all(merged_mixed, gtdb_tab, str(ai),
                    use_datasets_fallback=True, _resolver=resolver)
    assert called["datasets"] is False


# ─── present_accessions (suppression screen) ───────────────────────────────

def _write_ai(path, accessions):
    with open(path, "w") as fh:
        fh.write("#header\n")
        for a in accessions:
            fh.write(a + "\t" + "\t".join(["x"] * 22) + "\n")


def test_present_accessions_detects_absence(tmp_path):
    ai = tmp_path / "ncbi-assembly-info.tsv"
    _write_ai(ai, ["GCA_000000003.1", "GCA_000000004.1"])   # 002 absent (suppressed)
    live = TAX.present_accessions(
        ["GCA_000000002.1", "GCA_000000003.1", "GCA_000000004.1"], str(ai))
    assert live == {"GCA_000000003.1", "GCA_000000004.1"}


def test_present_accessions_matches_gca_gcf_twin(tmp_path):
    # file has the GenBank (GCA) row; a GCF query must still register present
    ai = tmp_path / "ncbi-assembly-info.tsv"
    _write_ai(ai, ["GCA_000005845.2"])
    live = TAX.present_accessions(["GCF_000005845.2"], str(ai))
    assert live == {"GCF_000005845.2"}


def test_present_accessions_version_tolerant(tmp_path):
    ai = tmp_path / "ncbi-assembly-info.tsv"
    _write_ai(ai, ["GCA_000005845.3"])           # newer version in file
    live = TAX.present_accessions(["GCA_000005845.1"], str(ai))   # older queried
    assert live == {"GCA_000005845.1"}


def test_present_accessions_missing_file():
    assert TAX.present_accessions(["GCA_000000001.1"], None) == set()
    assert TAX.present_accessions(["GCA_000000001.1"], "/no/such.tsv") == set()
