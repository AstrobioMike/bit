import pytest # type: ignore
from pathlib import Path

from bit.modules.ncbi.parse_ncbi_assembly_summary import (parse_ncbi_assembly_summary,
                                                          build_base_link,
                                                          sanitize_assembly_name)
from bit.modules.ncbi.dl_ncbi_assemblies import RunData


def _make_summary_line(acc, assembly_name="TestAssembly_v1", taxid="12345",
                       org="Test organism", infra="", version="latest",
                       level="Chromosome", ftp_path=""):
    fields = ["NA"] * 23
    fields[0] = acc
    fields[5] = taxid
    fields[7] = org
    fields[8] = infra
    fields[10] = version
    fields[11] = level
    fields[15] = assembly_name
    fields[19] = ftp_path if ftp_path else "na"
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
    assert local_dest.endswith(".fasta.gz")


def test_parse_empty_summary_file(tmp_path):
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text("")
    rd = _make_run_data(tmp_path, ["GCA_000001405.29"])
    rd = parse_ncbi_assembly_summary(summary, rd)
    assert rd.num_found == 0
    assert rd.num_not_found == 1


def test_ftp_path_used_when_present(tmp_path):
    # the METABAT case: assembly_name has a double underscore that NCBI collapsed
    # in the real directory; ftp_path must win so the URL matches NCBI's path
    real_ftp = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/938/034/415/"
                "GCA_938034415.1_S1_Bin_METABAT_151_1")
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_938034415.1", "S1_Bin_METABAT__151_1",
                           ftp_path=real_ftp) + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_938034415.1"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    lines = (tmp_path / "ncbi-info.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    target = lines[1].split("\t")[header.index("target_link")]
    # single underscore (from ftp_path), not the double from assembly_name
    assert "METABAT_151_1_genomic.fna.gz" in target
    assert "METABAT__151" not in target


def test_fallback_sanitizes_assembly_name(tmp_path):
    # ftp_path empty -> reconstruct, and the reconstruction must sanitize the
    # double underscore the same way NCBI would
    summary = tmp_path / "assembly_summary.tsv"
    summary.write_text(
        _make_summary_line("GCA_111222333.1", "Some_MAG__bin_2", ftp_path="") + "\n"
    )
    rd = _make_run_data(tmp_path, ["GCA_111222333.1"], wanted_format="fasta")
    parse_ncbi_assembly_summary(summary, rd)
    lines = (tmp_path / "ncbi-info.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    target = lines[1].split("\t")[header.index("target_link")]
    assert "GCA/111/222/333/GCA_111222333.1_Some_MAG_bin_2/" in target
    assert "Some_MAG__bin_2" not in target


def test_sanitize_assembly_name():
    assert sanitize_assembly_name("S1_Bin_METABAT__151_1") == "S1_Bin_METABAT_151_1"
    assert sanitize_assembly_name("ASM584v2") == "ASM584v2"
    assert sanitize_assembly_name("Release 6 plus ISO1/MT") == "Release_6_plus_ISO1_MT"
    assert sanitize_assembly_name("weird#name (v2) [draft]") == "weird_name_v2_draft"


################################################################################
# header-based (slim table) parsing — equivalence with legacy positional format
################################################################################

SLIM_COLUMNS = [
    "assembly_accession", "taxid", "organism_name", "infraspecific_name",
    "version_status", "assembly_level", "asm_name", "ftp_path",
]


def _slim_line(acc, taxid="12345", org="Test organism", infra="NA",
               version="latest", level="Chromosome", asm_name="TestAssembly_v1",
               ftp_path="na"):
    """A row in the slim (8-column, header-based) layout."""
    return "\t".join([acc, taxid, org, infra, version, level, asm_name, ftp_path])


def test_parse_slim_header_format_found(tmp_path):
    summary = tmp_path / "slim.tsv"
    summary.write_text(
        "\t".join(SLIM_COLUMNS) + "\n" +
        _slim_line("GCF_000005845.2", taxid="562", asm_name="ASM584v2",
                   ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2") + "\n"
    )
    run_data = _make_run_data(tmp_path, ["GCF_000005845.2"])
    run_data = parse_ncbi_assembly_summary(summary, run_data)
    assert run_data.num_found == 1
    out = (tmp_path / "ncbi-info.tsv").read_text().splitlines()
    # header + one data row
    assert len(out) == 2
    rec = dict(zip(out[0].split("\t"), out[1].split("\t")))
    assert rec["found_accession"] == "GCF_000005845.2"
    assert rec["taxid"] == "562"
    assert rec["http_base_link"].endswith("GCF_000005845.2_ASM584v2/")


def test_parse_slim_and_legacy_give_same_result(tmp_path):
    """The header-based slim file and the legacy positional file must yield
    identical parsed output for the same underlying record."""
    ftp = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"

    # legacy positional (23-field) file, no header
    legacy = tmp_path / "legacy.tsv"
    legacy.write_text(
        _make_summary_line("GCF_000005845.2", assembly_name="ASM584v2",
                           taxid="562", org="Escherichia coli", infra="strain=K-12",
                           version="latest", level="Complete Genome",
                           ftp_path=ftp) + "\n"
    )
    rd_legacy = _make_run_data(tmp_path, ["GCF_000005845.2"])
    rd_legacy.ncbi_sub_table_path = tmp_path / "legacy-out.tsv"
    parse_ncbi_assembly_summary(legacy, rd_legacy)

    # slim header-based file with the same values
    slim = tmp_path / "slim.tsv"
    slim.write_text(
        "\t".join(SLIM_COLUMNS) + "\n" +
        _slim_line("GCF_000005845.2", taxid="562", org="Escherichia coli",
                   infra="strain=K-12", version="latest", level="Complete Genome",
                   asm_name="ASM584v2", ftp_path=ftp) + "\n"
    )
    rd_slim = _make_run_data(tmp_path, ["GCF_000005845.2"])
    rd_slim.ncbi_sub_table_path = tmp_path / "slim-out.tsv"
    parse_ncbi_assembly_summary(slim, rd_slim)

    assert (tmp_path / "legacy-out.tsv").read_text() == (tmp_path / "slim-out.tsv").read_text()


def test_parse_slim_not_found_recorded(tmp_path):
    summary = tmp_path / "slim.tsv"
    summary.write_text(
        "\t".join(SLIM_COLUMNS) + "\n" +
        _slim_line("GCF_000005845.2") + "\n"
    )
    run_data = _make_run_data(tmp_path, ["GCF_999999999.1"])
    run_data = parse_ncbi_assembly_summary(summary, run_data)
    assert run_data.num_found == 0
    assert run_data.num_not_found == 1
