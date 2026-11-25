from bit.modules.combine_coverage_tables import combine_gene_tables, combine_contig_tables, read_table

RANK_HEADER = "domain\tphylum\tclass\torder\tfamily\tgenus\tspecies"
ECOLI = "Bacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli"
ESCH = "Bacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tNA"
NONE = "\t".join(["NA"] * 7)


def write(path, lines):
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_contig_tables_summed_by_lineage(tmp_path):
    a = write(tmp_path / "a.tsv", ["contig_ID\tcoverage\t" + RANK_HEADER + "\tevidence",
                                   f"c1\t10\t{ECOLI}\tk-mers", f"c2\t5\t{ECOLI}\tk-mers",
                                   f"c3\t3\t{ESCH}\tproteins", f"c4\t2\t{NONE}\tnone"])
    b = write(tmp_path / "b.tsv", ["contig_ID\tcoverage\t" + RANK_HEADER + "\tevidence",
                                   f"c1\t4\t{ESCH}\tproteins"])
    empty = write(tmp_path / "empty.tsv", ["contig_ID\tcoverage\t" + RANK_HEADER + "\tevidence"])
    out, cpm = tmp_path / "out.tsv", tmp_path / "cpm.tsv"

    combine_contig_tables({"A": a, "B": b, "E": empty}, str(out), str(cpm))

    rows = read_table(out)
    # the sample without contigs is left out; species-level before genus-level, unclassified last
    assert list(rows[0]) == RANK_HEADER.split("\t") + ["A", "B"]
    assert [(r["species"], r["genus"], r["A"], r["B"]) for r in rows] == [
        ("Escherichia coli", "Escherichia", "15", "0"),
        ("NA", "Escherichia", "3", "4"),
        ("NA", "NA", "2", "0"),
    ]
    assert [float(r["A"]) for r in read_table(cpm)] == [750000, 150000, 100000]
    assert [float(r["B"]) for r in read_table(cpm)] == [0, 1000000, 0]


def test_gene_tables_summed_by_KO_and_lineage(tmp_path):
    header = "gene_ID\tcoverage\tKO_ID\tKO_function\t" + RANK_HEADER
    a = write(tmp_path / "a.tsv", [header, f"g1\t1\tK00001\talcohol dehydrogenase\t{ECOLI}",
                                   f"g2\t2\tK00001\talcohol dehydrogenase\t{ESCH}",
                                   f"g3\t3\tNA\tNA\t{ECOLI}"])
    b = write(tmp_path / "b.tsv", [header, f"g1\t2\tK00002\taldehyde reductase\t{NONE}",
                                   f"g2\t2\tNA\tNA\t{NONE}"])
    paths = [str(tmp_path / n) for n in ("ko.tsv", "ko-cpm.tsv", "tax.tsv", "tax-cpm.tsv")]

    combine_gene_tables({"A": a, "B": b}, *paths)

    KOs = read_table(paths[0])
    assert [(r["KO_ID"], r["KO_function"], r["A"], r["B"]) for r in KOs] == [
        ("K00001", "alcohol dehydrogenase", "3", "0"),
        ("K00002", "aldehyde reductase", "0", "2"),
        ("Not annotated", "Not annotated", "3", "2"),
    ]
    tax = read_table(paths[2])
    assert [(r["species"], r["genus"], r["A"], r["B"]) for r in tax] == [
        ("Escherichia coli", "Escherichia", "4", "0"),
        ("NA", "Escherichia", "2", "0"),
        ("NA", "NA", "0", "4"),
    ]


def test_identical_coverage_rows_are_kept_apart(tmp_path):
    """ two KOs with the same coverage in every sample stay two rows """

    header = "gene_ID\tcoverage\tKO_ID\tKO_function\t" + RANK_HEADER
    a = write(tmp_path / "a.tsv", [header, f"g1\t2\tK00001\tf1\t{ECOLI}", f"g2\t2\tK00002\tf2\t{ESCH}"])
    paths = [str(tmp_path / n) for n in ("ko.tsv", "ko-cpm.tsv", "tax.tsv", "tax-cpm.tsv")]

    combine_gene_tables({"A": a}, *paths)

    assert [r["KO_ID"] for r in read_table(paths[0])] == ["K00001", "K00002"]
    assert len(read_table(paths[2])) == 2
