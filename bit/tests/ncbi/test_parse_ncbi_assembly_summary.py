import pytest # type: ignore
from pathlib import Path

from bit.modules.ncbi.parse_ncbi_assembly_summary import parse_ncbi_assembly_summary, build_base_link
from bit.modules.dl_ncbi_assemblies import RunData


def _make_summary_line(acc, assembly_name="TestAssembly_v1", taxid="12345",
                       org="Test organism", infra="", version="latest",
                       level="Chromosome"):
    fields = ["NA"] * 16
    fields[0] = acc
    fields[5] = taxid
    fields[7] = org
    fields[8] = infra
    fields[10] = version
    fields[11] = level
    fields[15] = assembly_name
    return "\t".join(fields)


def _make_run_data(tmp_path, wanted_accs, wanted_format=None):
    return RunData(
        wanted_accs=wanted_accs,
        num_wanted=len(wanted_accs),
        wanted_format=wanted_format,
        output_dir=str(tmp_path),
        ncbi_sub_table_path=tmp_path / "ncbi-info.tsv",
        not_found_path=tmp_path / "not-found.txt",
    )


def test_build_base_link_url_structure():
    url = build_base_link("GCF_000005845", "ASM584v2")
    assert url.startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/")
    assert "GCF/000/005/845/" in url
    assert url.endswith("GCF_000005845_ASM584v2/")


def test_build_base_link_prefix_preserved():
    gca_url = build_base_link("GCA_000001405", "GRCh38")
    gcf_url = build_base_link("GCF_000001405", "GRCh38")
    assert "/GCA/" in gca_url
    assert "/GCF/" in gcf_url


def test_build_base_link_path_segments():
    # GCA_123456789 → p1=123, p2=456, p3=789
    url = build_base_link("GCA_123456789", "MyAssembly")
    assert "/GCA/123/456/789/" in url


def test_parse_all_found(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n" +
        _make_summary_line("GCF_000005845.2", "ASM584v2") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29", "GCF_000005845.2"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 2
    assert rd.num_not_found == 0
    assert not (tmp_path / "not-found.txt").exists()


def test_parse_partial_not_found(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29", "GCF_000005845.2"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 1
    assert rd.num_not_found == 1
    not_found_text = (tmp_path / "not-found.txt").read_text()
    assert "GCF_000005845.2" in not_found_text


def test_parse_none_found(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_999999999.1", "SomeAssembly") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 0
    assert rd.num_not_found == 1


def test_parse_version_stripping(tmp_path):
    # accession GCA_000001405.1 in wanted should match GCA_000001405.29 in file
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.1"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 1


def test_parse_output_tsv_header(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"])
    parse_ncbi_assembly_summary(summary, rd)
    header = (tmp_path / "ncbi-info.tsv").read_text().splitlines()[0].split("\t")
    assert header[0] == "input_accession"
    assert "organism_name" in header
    assert "assembly_level" in header
    assert "target_link" not in header


def test_parse_output_tsv_row_values(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14", taxid="9606",
                           org="Homo sapiens", level="Chromosome") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"])
    parse_ncbi_assembly_summary(summary, rd)
    lines = (tmp_path / "ncbi-info.tsv").read_text().splitlines()
    assert len(lines) == 2  # header + one data row
    fields = lines[1].split("\t")
    header = lines[0].split("\t")
    assert fields[header.index("input_accession")] == "GCA_000001405.29"
    assert fields[header.index("taxid")] == "9606"
    assert fields[header.index("organism_name")] == "Homo sapiens"
    assert fields[header.index("assembly_level")] == "Chromosome"


def test_parse_with_format_adds_link_columns(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    header = (tmp_path / "ncbi-info.tsv").read_text().splitlines()[0].split("\t")
    assert "target_link" in header
    assert "local_destination" in header


def test_parse_with_format_link_content(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_000001405.29", "GRCh38p14") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    lines = (tmp_path / "ncbi-info.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    target_link = row[header.index("target_link")]
    local_dest = row[header.index("local_destination")]
    assert "_genomic.fna.gz" in target_link
    assert "ftp.ncbi.nlm.nih.gov" in target_link
    assert local_dest.endswith(".fa.gz")


def test_parse_empty_summary_file(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text("")
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 0
    assert rd.num_not_found == 1
