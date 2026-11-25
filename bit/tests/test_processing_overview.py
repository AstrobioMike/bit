from bit.modules.processing_overview import write_processing_overview, total_reads, contigs_in


def own(tmp_path, *samples):
    """ each sample's reads mapped to its own assembly """
    return {s: (s, tmp_path / "assemblies" / f"{s}-assembly.fasta") for s in samples}


def make_sample(tmp_path, sample, pairs, mapped):
    for d in ("assemblies", "genes", "mapping", "bins", "MAGs"):
        (tmp_path / d).mkdir(exist_ok=True)
    (tmp_path / "assemblies" / f"{sample}-assembly.fasta").write_text(">c1\nACGT\n")
    (tmp_path / "genes" / f"{sample}-genes.faa").write_text(">g1\nM\n")
    (tmp_path / "mapping" / f"{sample}.bam").write_text("bam")
    (tmp_path / "mapping" / f"{sample}-mapping-info.txt").write_text(f"{pairs} reads; of these:\n  ...\n")
    (tmp_path / "mapping" / f"{sample}-mapped-reads-per-contig.tsv").write_text(
        "".join(f"{c}\t1000\t{n}\t0\n" for c, n in mapped.items()) + "*\t0\t0\t7\n")


def test_percent_of_reads_mapped_to_assembly_bins_and_MAGs(tmp_path):
    make_sample(tmp_path, "s1", 50, {"c1": 40, "c2": 30, "c3": 20})        # 100 reads (pairs x 2), 90 mapped
    (tmp_path / "bins" / "s1-bin.1.fasta").write_text(">c1 total_depth=3\nACGT\n>c2\nACGT\n")
    (tmp_path / "bins" / "s1-bin.2.fasta").write_text(">c3\nACGT\n")
    (tmp_path / "MAGs" / "s1-MAG-1.fasta").write_text(">c1 total_depth=3\nACGT\n>c2\nACGT\n")
    out = tmp_path / "overview.tsv"

    write_processing_overview(["s1"], out, own(tmp_path, "s1"), tmp_path / "mapping", True, genes_dir = tmp_path / "genes",
                              bins_dir = tmp_path / "bins", MAGs_dir = tmp_path / "MAGs")

    header, row = [l.split("\t") for l in out.read_text().splitlines()]
    values = dict(zip(header, row))
    assert values["reads"] == "100"
    assert values["percent_of_reads_mapped_to_assembly"] == "90.00"
    assert values["percent_of_reads_mapped_to_bins"] == "90.00"
    assert values["percent_of_reads_mapped_to_MAGs"] == "70.00"
    assert (values["bins_recovered"], values["MAGs_recovered"]) == ("2", "1")


def test_without_binning_the_bin_and_MAG_columns_are_left_out(tmp_path):
    make_sample(tmp_path, "s1", 10, {"c1": 5})                              # single-end: 10 reads
    out = tmp_path / "overview.tsv"

    write_processing_overview(["s1"], out, own(tmp_path, "s1"), tmp_path / "mapping", False, genes_dir = tmp_path / "genes")

    header, row = [l.split("\t") for l in out.read_text().splitlines()]
    assert header[-2:] == ["reads", "percent_of_reads_mapped_to_assembly"]
    assert row[-2:] == ["10", "50.00"]
    assert not any("MAG" in h for h in header)


def test_empty_mapping_summary_gives_NA(tmp_path):
    make_sample(tmp_path, "s1", 0, {})
    (tmp_path / "mapping" / "s1-mapping-info.txt").write_text("")
    assert total_reads(tmp_path / "mapping" / "s1-mapping-info.txt", True) == 0
    out = tmp_path / "overview.tsv"
    write_processing_overview(["s1"], out, own(tmp_path, "s1"), tmp_path / "mapping", True, genes_dir = tmp_path / "genes")
    assert out.read_text().splitlines()[1].endswith("\t0\tNA")


def test_a_co_assembly_without_gene_calls(tmp_path):
    """ magrec with one co-assembly: every sample's reads went to it, and its bins and MAGs count for each """
    for sample, mapped in (("s1", {"c1": 40, "c2": 10}), ("s2", {"c1": 5, "c2": 45})):
        make_sample(tmp_path, sample, 50, mapped)
    (tmp_path / "assemblies" / "co.fasta").write_text(">c1\nACGT\n>c2\nACGT\n")
    (tmp_path / "bins" / "co-bin.1.fasta").write_text(">c1\nACGT\n")
    (tmp_path / "MAGs" / "co-MAG-1.fasta").write_text(">c1\nACGT\n")
    out = tmp_path / "overview.tsv"

    co = {s: ("co", tmp_path / "assemblies" / "co.fasta") for s in ("s1", "s2")}
    write_processing_overview(["s1", "s2"], out, co, tmp_path / "mapping", True,
                              bins_dir = tmp_path / "bins", MAGs_dir = tmp_path / "MAGs")

    header, *rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert "gene_calls_identified" not in header
    s1, s2 = (dict(zip(header, row)) for row in rows)
    assert (s1["bins_recovered"], s2["MAGs_recovered"]) == ("1", "1")
    assert (s1["percent_of_reads_mapped_to_MAGs"], s2["percent_of_reads_mapped_to_MAGs"]) == ("40.00", "5.00")


def test_bowtie2_warnings_before_its_summary_are_skipped(tmp_path):
    info = tmp_path / "info.txt"
    info.write_text("Warning: skipping mate #1 of read 'r1/1' because length (1) <= # seed mismatches (0)\n"
                    "Warning: skipping mate #1 of read 'r1/1' because it was < 2 characters long\n"
                    "716144 reads; of these:\n  716144 (100.00%) were paired; of these:\n")
    assert total_reads(info, True) == 1432288
