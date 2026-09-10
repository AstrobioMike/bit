from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest  # type: ignore

from bit.modules.ncbi.parse_ncbi_assembly_summary import (parse_ncbi_assembly_summary,
                                                          build_base_link,
                                                          _resolve_links,
                                                          sanitize_assembly_name)
from bit.modules.ncbi.dl_ncbi_assemblies import RunData
from bit.modules.taxonomy.tax_ranks import RANKS, accession_core
from bit.modules.taxonomy.lineage_lookup import NO_LINEAGE


# The reader now consumes the hosted NCBI Parquet, so fixtures are Parquet rather than
# TSV. Everything else -- the RunData contract, the output columns, the link logic --
# is unchanged and asserted identical to the old line-scan reader.

_PARQUET_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]


def _row(acc, assembly_name="TestAssembly_v1", taxid="12345", org="Test organism",
         infra="", version="latest", level="Chromosome", ftp_path=""):
    return {
        "assembly_accession": acc,
        "asm_name": assembly_name,
        "taxid": taxid,
        "organism_name": org,
        "infraspecific_name": infra,
        "version_status": version,
        "assembly_level": level,
        "ftp_path": ftp_path if ftp_path else "na",
    }


def _make_summary(tmp_path, rows):
    """Write a Parquet table from a list of _row() dicts."""
    cols = {c: pa.array([str(r.get(c, "")) for r in rows]) for c in _PARQUET_COLUMNS}
    path = tmp_path / "ncbi-data.parquet"
    pq.write_table(pa.table(cols), str(path))
    return path


def _make_run_data(tmp_path, wanted_accs, wanted_format=None):
    return RunData(
        wanted_accs=wanted_accs,
        num_wanted=len(wanted_accs),
        wanted_format=wanted_format,
        output_dir=str(tmp_path),
        ncbi_sub_table_path=tmp_path / "ncbi-info.tsv",
        not_found_path=tmp_path / "not-found.txt",
    )


# --- build_base_link (unchanged; the twice-broken seam) -------------------

def test_build_base_link_url_structure():
    url, dir_basename = build_base_link("GCF_000005845", "ASM584v2")
    assert url.startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/")
    assert "GCF/000/005/845/" in url
    assert url.endswith("GCF_000005845_ASM584v2/")
    assert dir_basename == "GCF_000005845_ASM584v2"


def test_build_base_link_prefix_preserved():
    gca_url, _ = build_base_link("GCA_000001405", "GRCh38")
    gcf_url, _ = build_base_link("GCF_000001405", "GRCh38")
    assert "/GCA/" in gca_url
    assert "/GCF/" in gcf_url


def test_build_base_link_path_segments():
    url, _ = build_base_link("GCA_123456789", "MyAssembly")
    assert "/GCA/123/456/789/" in url


def test_sanitize_assembly_name():
    assert sanitize_assembly_name("A B/C,D") == "A_B_C_D"
    assert sanitize_assembly_name("x[]()#") == "x"


# --- core parsing ---------------------------------------------------------

def test_parse_all_found(tmp_path):
    summary = _make_summary(tmp_path, [
        _row("GCA_000001405.29", "GRCh38p14"),
        _row("GCF_000005845.2", "ASM584v2"),
    ])
    rd = _make_run_data(tmp_path, ["GCA_000001405.29", "GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 2
    assert rd.num_not_found == 0
    assert not rd.not_found_path.exists()


def test_parse_partial_not_found(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM584v2")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2", "GCF_999999999.1"])
    parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 1
    assert rd.num_not_found == 1
    assert rd.not_found_path.read_text().strip() == "GCF_999999999.1"


def test_parse_none_found(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2")])
    rd = _make_run_data(tmp_path, ["GCF_111111111.1"])
    parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 0
    assert rd.num_not_found == 1


def test_parse_version_stripping(tmp_path):
    """A wanted acc matches regardless of the version the table happens to hold."""
    summary = _make_summary(tmp_path, [_row("GCF_000005845.7", "ASM584v2")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])       # different version
    parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 1
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert row[0] == "GCF_000005845.2"        # target_accession = what the user asked
    assert row[1] == "GCF_000005845.7"        # found_accession = what the table had


def test_over_match_prefix_is_rejected(tmp_path):
    """
    starts_with('GCF_000000001') would also hit 'GCF_0000000019'. The exact-root
    check after stripping the version must reject that.
    """
    summary = _make_summary(tmp_path, [
        _row("GCF_000000001.1", "A"),
        _row("GCF_0000000019.1", "B"),      # longer number, shares the prefix
    ])
    rd = _make_run_data(tmp_path, ["GCF_000000001.1"])
    parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 1
    accs = [l.split("\t")[1] for l in rd.ncbi_sub_table_path.read_text().splitlines()[1:]]
    assert accs == ["GCF_000000001.1"]


def test_empty_wanted_accs_writes_header_only(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2")])
    rd = _make_run_data(tmp_path, [])
    parse_ncbi_assembly_summary(summary, rd)
    assert len(rd.ncbi_sub_table_path.read_text().splitlines()) == 1     # header only
    assert rd.num_found == 0


# --- output shape ---------------------------------------------------------

def test_output_tsv_header(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    assert header == ["target_accession", "found_accession", "assembly_name", "taxid",
                      "organism_name", "infraspecific_name", "version_status",
                      "assembly_level"]


def test_output_tsv_row_values(tmp_path):
    summary = _make_summary(tmp_path, [
        _row("GCF_000005845.2", "ASM584v2", taxid="562", org="Escherichia coli",
             infra="strain=K-12", level="Complete Genome",
             ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"),
    ])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert row[2] == "ASM584v2"
    assert row[3] == "562"
    assert row[4] == "Escherichia coli"
    assert row[5] == "strain=K-12"
    assert row[7] == "Complete Genome"


def test_empty_fields_become_NA(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", infra="")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert row[5] == "NA"       # infraspecific_name was empty


# --- format-specific columns ----------------------------------------------

def test_with_format_adds_link_columns(tmp_path):
    summary = _make_summary(tmp_path, [
        _row("GCF_000005845.2", "ASM584v2",
             ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"),
    ])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    assert header[-2:] == ["target_link", "local_destination"]


def test_with_format_link_content(tmp_path):
    summary = _make_summary(tmp_path, [
        _row("GCF_000005845.2", "ASM584v2",
             ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"),
    ])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert row[-1] == f"{tmp_path}/GCF_000005845.2.fasta.gz"
    assert row[-2].endswith("GCF_000005845.2_ASM584v2_genomic.fna.gz")


# --- link resolution: ftp_path vs the build_base_link fallback ------------

# http_base_link is no longer a column, so these pin the resolver itself plus the
# target_link it feeds -- the behaviour that column used to stand in for

def test_ftp_path_used_when_present(tmp_path):
    link, _ = _resolve_links(
        "GCF_000005845.2", "ASM584v2",
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2")
    assert link.startswith("https://")            # ftp:// -> https://
    assert link.endswith("/")


def test_fallback_when_ftp_path_absent():
    """No ftp_path -> the URL is rebuilt from accession + assembly name."""
    link, _ = _resolve_links("GCF_000005845.2", "ASM584v2", "")
    assert link == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/"


def test_fallback_sanitizes_assembly_name(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM 584/v2", ftp_path="")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert "ASM_584_v2" in row[header.index("target_link")]


def test_http_base_link_column_is_gone(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    assert "http_base_link" not in header
    assert "input_accession" not in header


def test_local_destination_has_no_dot_slash_prefix(tmp_path, monkeypatch):
    """
    The default output dir is ".", which an f-string join turned into a "./" on the
    front of every destination.
    """
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM584v2")])
    monkeypatch.chdir(tmp_path)
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    rd.output_dir = "."
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    dest = row[header.index("local_destination")]
    assert dest == "GCF_000005845.2.fasta.gz"
    assert not dest.startswith("./")


def test_explicit_output_dir_still_prefixes_destination(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM584v2")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")
    rd.output_dir = "genomes"
    parse_ncbi_assembly_summary(summary, rd)
    header = rd.ncbi_sub_table_path.read_text().splitlines()[0].split("\t")
    row = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")
    assert row[header.index("local_destination")] == "genomes/GCF_000005845.2.fasta.gz"


# --- lineage columns in the info table ------------------------------------
# --add-ncbi-tax / --add-gtdb-tax are independent of each other and of --source, and
# both default off. NCBI lineage rides along on the asset scan; GTDB lineage comes
# from a map the caller builds and hangs off run_data.

_LINEAGE = ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
            "Enterobacteriaceae", "Escherichia", "Escherichia coli")


def _lineage_summary(tmp_path, with_ranks=True):
    """A one-row Parquet fixture, optionally carrying the asset's lineage columns."""
    row = _row("GCF_000005845.2", "ASM584v2",
               ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
                        "GCF_000005845.2_ASM584v2")
    cols = {c: pa.array([str(row.get(c, ""))]) for c in _PARQUET_COLUMNS}
    if with_ranks:
        for rank, value in zip(RANKS, _LINEAGE):
            cols[rank] = pa.array([value])
    path = tmp_path / "ncbi-data-lineage.parquet"
    pq.write_table(pa.table(cols), str(path))
    return path


def _lineage_run_data(tmp_path):
    return _make_run_data(tmp_path, ["GCF_000005845.2"], wanted_format="fasta")


_lineage_parse = parse_ncbi_assembly_summary


def _lineage_table(tmp_path, add_ncbi_tax=False, add_gtdb_tax=False,
                   gtdb_lineage=None, with_ranks=True):
    summary = _lineage_summary(tmp_path, with_ranks=with_ranks)
    rd = _lineage_run_data(tmp_path)
    rd.add_ncbi_tax = add_ncbi_tax
    rd.add_gtdb_tax = add_gtdb_tax
    rd.gtdb_lineage = gtdb_lineage
    _lineage_parse(summary, rd)
    lines = Path(rd.ncbi_sub_table_path).read_text().splitlines()
    header = lines[0].split("\t")
    return header, [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def test_no_lineage_columns_by_default(tmp_path):
    header, _ = _lineage_table(tmp_path)
    assert not [c for c in header if c.startswith(("ncbi_", "gtdb_"))]


def test_add_ncbi_tax_adds_prefixed_ncbi_columns(tmp_path):
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True)
    assert header[-7:] == [f"ncbi_{r}" for r in RANKS]
    assert rows[0]["ncbi_species"] == "Escherichia coli"
    assert rows[0]["ncbi_domain"] == "Bacteria"
    assert not [c for c in header if c.startswith("gtdb_")]


def test_add_gtdb_tax_adds_prefixed_gtdb_columns(tmp_path):
    mapping = {accession_core("GCF_000005845.2"): _LINEAGE}
    header, rows = _lineage_table(tmp_path, add_gtdb_tax=True, gtdb_lineage=mapping,
                                  with_ranks=False)
    assert header[-7:] == [f"gtdb_{r}" for r in RANKS]
    assert rows[0]["gtdb_species"] == "Escherichia coli"
    assert not [c for c in header if c.startswith("ncbi_")]


def test_both_taxonomies_can_be_on_at_once(tmp_path):
    """
    They're independent flags, and the prefixes are what keeps a mixed table
    unambiguous when GTDB and NCBI disagree.
    """
    gtdb = ("Bacteria", "GtdbPhylum", "GtdbClass", "GtdbOrder", "GtdbFamily",
            "GtdbGenus", "Gtdb species")
    mapping = {accession_core("GCF_000005845.2"): gtdb}
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True, add_gtdb_tax=True,
                                  gtdb_lineage=mapping)
    assert header[-14:] == ([f"ncbi_{r}" for r in RANKS] +
                            [f"gtdb_{r}" for r in RANKS])
    assert rows[0]["ncbi_phylum"] == "Pseudomonadota"
    assert rows[0]["gtdb_phylum"] == "GtdbPhylum"


def test_accession_missing_from_gtdb_gets_NA(tmp_path):
    """GTDB is bacteria/archaea only, so a miss is expected, not an error."""
    header, rows = _lineage_table(tmp_path, add_gtdb_tax=True, gtdb_lineage={},
                                  with_ranks=False)
    assert [rows[0][f"gtdb_{r}"] for r in RANKS] == list(NO_LINEAGE)


def test_lineage_columns_come_after_the_link_columns(tmp_path):
    header, _ = _lineage_table(tmp_path, add_ncbi_tax=True)
    assert header.index("local_destination") < header.index("ncbi_domain")
    assert header.index("target_link") < header.index("ncbi_domain")


def test_every_row_has_the_full_width(tmp_path):
    mapping = {accession_core("GCF_000005845.2"): _LINEAGE}
    header, rows = _lineage_table(tmp_path, add_ncbi_tax=True, add_gtdb_tax=True,
                                  gtdb_lineage=mapping)
    assert all(len(row) == len(header) for row in rows)
    assert "" not in rows[0].values()
