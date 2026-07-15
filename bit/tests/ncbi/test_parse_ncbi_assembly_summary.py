import pyarrow as pa
import pyarrow.parquet as pq
import pytest  # type: ignore

from bit.modules.ncbi.parse_ncbi_assembly_summary import (parse_ncbi_assembly_summary,
                                                          build_base_link,
                                                          sanitize_assembly_name)
from bit.modules.ncbi.dl_ncbi_assemblies import RunData


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
    assert row[0] == "GCF_000005845.2"        # input_accession = what the user asked
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
    assert header == ["input_accession", "found_accession", "assembly_name", "taxid",
                      "organism_name", "infraspecific_name", "version_status",
                      "assembly_level", "http_base_link"]


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

def test_ftp_path_used_when_present(tmp_path):
    summary = _make_summary(tmp_path, [
        _row("GCF_000005845.2", "ASM584v2",
             ftp_path="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"),
    ])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    link = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")[8]
    assert link.startswith("https://")            # ftp:// -> https://
    assert link.endswith("/")


def test_fallback_when_ftp_path_absent(tmp_path):
    """No ftp_path -> the URL is rebuilt from accession + assembly name."""
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM584v2", ftp_path="")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    link = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")[8]
    assert link == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/"


def test_fallback_sanitizes_assembly_name(tmp_path):
    summary = _make_summary(tmp_path, [_row("GCF_000005845.2", "ASM 584/v2", ftp_path="")])
    rd = _make_run_data(tmp_path, ["GCF_000005845.2"])
    parse_ncbi_assembly_summary(summary, rd)
    link = rd.ncbi_sub_table_path.read_text().splitlines()[1].split("\t")[8]
    assert "ASM_584_v2" in link
