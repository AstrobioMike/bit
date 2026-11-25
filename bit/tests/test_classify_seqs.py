from bit.modules.classify_seqs import (common_lineage, deepest_rank_allowed, label_gene, label_from_genes,
                                       label_from_kmers, classify)

ECOLI = ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli")
EFERG = ECOLI[:6] + ("Escherichia fergusonii",)
SALMO = ECOLI[:5] + ("Salmonella", "Salmonella enterica")
BSUBT = ("Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus", "Bacillus subtilis")


def test_common_lineage():
    assert common_lineage([ECOLI, EFERG]) == ECOLI[:6]
    assert common_lineage([ECOLI, BSUBT]) == ("Bacteria",)
    assert common_lineage([ECOLI]) == ECOLI
    assert common_lineage([]) == ()
    assert common_lineage([ECOLI, ECOLI[:3]]) == ECOLI[:3]


def test_deepest_rank_allowed():
    assert deepest_rank_allowed(99) == 7
    assert deepest_rank_allowed(95) == 7
    assert deepest_rank_allowed(94.9) == 6
    assert deepest_rank_allowed(80) == 6
    assert deepest_rank_allowed(60) == 5
    assert deepest_rank_allowed(40) == 4


def test_label_gene_uses_only_near_best_hits():
    lineage, score = label_gene([(ECOLI, 500, 99), (EFERG, 300, 90), (BSUBT, 100, 50)])
    assert lineage == ECOLI and score == 500

    lineage, _ = label_gene([(ECOLI, 500, 99), (EFERG, 460, 97)])
    assert lineage == ECOLI[:6]


def test_label_gene_identity_caps_the_rank():
    assert label_gene([(ECOLI, 500, 85)])[0] == ECOLI[:6]
    assert label_gene([(ECOLI, 500, 70)])[0] == ECOLI[:5]
    assert label_gene([(ECOLI, 500, 45)])[0] == ECOLI[:4]
    assert label_gene([]) == ((), 0.0)


def test_a_near_best_hit_to_a_root_cluster_leaves_the_gene_unlabelled():
    ROOT = ()   # a cluster whose members span domains has an empty lineage
    assert label_gene([(ECOLI, 500, 99), (ROOT, 480, 90)]) == ((), 500)
    assert label_gene([(ECOLI, 500, 99), (ROOT, 300, 60)]) == (ECOLI, 500)
    assert label_from_genes([label_gene([(ROOT, 500, 99)]), (ECOLI, 100)]) == ECOLI


def test_label_from_genes_needs_more_than_half_the_score():
    assert label_from_genes([(ECOLI, 300), (EFERG, 200), (SALMO, 100)]) == ECOLI[:6]
    assert label_from_genes([(ECOLI, 400), (SALMO, 100)]) == ECOLI
    assert label_from_genes([(ECOLI, 100), (BSUBT, 100)]) == ("Bacteria",)
    assert label_from_genes([((), 0.0)]) == ()
    assert label_from_genes([]) == ()


def test_shared_filled_ranks_are_not_shared_taxa():
    fill = ("unclassified Viruses phylum", "unclassified Viruses class", "unclassified Viruses order",
            "unclassified Viruses family", "unclassified Viruses genus")
    virus_a = ("Viruses",) + fill + ("Hubei permutotetra-like virus 5",)
    virus_b = ("Viruses",) + fill + ("Lake Sarah-associated circular molecule 2",)
    assert label_gene([(virus_a, 100, 99), (virus_b, 99, 99)]) == (("Viruses",), 100)
    assert label_gene([(virus_a, 100, 99)]) == (virus_a, 100)
    assert label_from_genes([(virus_a, 45), (virus_b, 45)]) == ("Viruses",)
    # a filled rank inside a lineage stays when a named taxon below it is shared
    euk = ("Eukaryota", "Ascomycota", "Dothideomycetes", "Pleosporales", "unclassified Pleosporales family", "Aaosphaeria")
    assert label_from_genes([(euk + ("Aaosphaeria arxii",), 45), (euk + ("Aaosphaeria sp.",), 45)]) == euk


def test_label_from_kmers():
    assert label_from_kmers([(ECOLI, 40), (SALMO, 5)], 50) == ECOLI
    # not enough of the sequence's k-mers
    assert label_from_kmers([(ECOLI, 20)], 50) == ()
    # too few shared k-mers
    assert label_from_kmers([(ECOLI, 2)], 2) == ()
    # two species equally good: no species call
    assert label_from_kmers([(ECOLI, 40), (EFERG, 39)], 50) == ()
    # several genomes of the one species
    assert label_from_kmers([(ECOLI, 40), (ECOLI, 38)], 50) == ECOLI
    assert label_from_kmers([], 50) == ()


def test_classify_prefers_kmers():
    assert classify(ECOLI, ECOLI[:6]) == (ECOLI, "k-mers")
    assert classify((), ECOLI[:6]) == (ECOLI[:6], "proteins")
    assert classify((), ()) == ((), "none")


from bit.modules.classify_seqs import read_lineages, read_protein_hits, classify_seqs, write_labels, seq_of_gene


def test_read_lineages(tmp_path):
    path = tmp_path / "lineages.tsv"
    path.write_text("17\t" + ";".join(ECOLI) + "\nGCF_1.1\t" + ";".join(BSUBT[:6]) + "\n")
    assert read_lineages(path) == {"17": ECOLI, "GCF_1.1": BSUBT[:6]}


def test_seq_of_gene_keeps_underscores_in_the_sequence_id():
    assert seq_of_gene("sample_1_contig_12_3") == "sample_1_contig_12"


def test_read_protein_hits_groups_by_gene_and_drops_unknown_proteins(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("c1_1\t17.5\t98.0\t400\nc1_1\t99.1\t60.0\t100\nc1_1\t1.3\t55.0\t90\nc1_2\t18.2\t70.5\t250\nc2_1\t99.9\t50\t80\n")
    hits = list(read_protein_hits(path, {"17": ECOLI, "18": SALMO, "1": ()}))
    assert hits == [("c1_1", [(ECOLI, 400.0, 98.0), ((), 90.0, 55.0)]), ("c1_2", [(SALMO, 250.0, 70.5)]), ("c2_1", [])]


def test_classify_seqs_and_the_table_written(tmp_path):
    protein_hits = [("c1_1", [(ECOLI, 400.0, 99.0)]), ("c1_2", [(EFERG, 100.0, 99.0)]),
                    ("c2_1", [(SALMO, 300.0, 85.0)]),
                    ("c3_1", [(BSUBT, 300.0, 99.0)])]
    kmer_matches = [("c3", "gE", 45), ("c2", "gE", 2)]
    labels = classify_seqs(["c1", "c2", "c3", "c4"], protein_hits, kmer_matches, {"c2": 50, "c3": 50}, {"gE": ECOLI})

    assert labels["c1"] == (ECOLI, "proteins")           # most of the score is on one species
    assert labels["c2"] == (SALMO[:6], "proteins")       # 85% identity: genus at most; too few k-mers shared
    assert labels["c3"] == (ECOLI, "k-mers")             # k-mers decide when they support a species
    assert labels["c4"] == ((), "none")

    out = tmp_path / "tax.tsv"
    write_labels(out, ["c1", "c2", "c3", "c4"], labels)
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0] == ["seq_id", "domain", "phylum", "class", "order", "family", "genus", "species", "evidence"]
    assert rows[2] == ["c2", *SALMO[:6], "NA", "proteins"]
    assert rows[4] == ["c4"] + ["NA"] * 7 + ["none"]


def test_label_inputs_writes_a_table_per_input(tmp_path):
    import gzip, json, zipfile
    import numpy as np # type: ignore
    from bit.modules.kmer_lookup import write_table
    from bit.modules.classify_seqs import kmer_step, label_inputs

    write_table(str(tmp_path / "kmers"), [("gE", list(range(100, 150)))])
    (tmp_path / "genome-lineages.tsv").write_text("gE\t" + ";".join(ECOLI) + "\n")
    (tmp_path / "protein-lineages.tsv").write_text("7\t" + ";".join(SALMO) + "\n")

    def sketch(path, seqs):
        records = [{"name": name + " some description", "signatures": [{"mins": mins}]} for name, mins in seqs]
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(records).encode()))

    (tmp_path / "a.fasta").write_text(">c1 len=10\nACGT\n>c2\nACGT\n")
    (tmp_path / "b.fasta").write_text(">c1\nACGT\n")
    sketch(tmp_path / "a.zip", [("c1", list(range(100, 140))), ("c2", [1, 2, 3])])
    sketch(tmp_path / "b.zip", [("c1", [5, 6])])
    (tmp_path / "a.faa").write_text(">a__c1_1\nMK\n>a__c2_1\nMK\n")
    (tmp_path / "b.faa").write_text(">b__c1_1\nMK\n")
    # the same sequence name in both inputs: kept apart by the input name on the genes
    (tmp_path / "hits.tsv").write_text("b__c1_1\t7.1\t99.0\t300\na__c2_1\t7.2\t85.0\t200\n")

    genes = {"a": tmp_path / "a.faa", "b": tmp_path / "b.faa"}
    # gene tables wanted: every gene is searched
    written = kmer_step(inputs={"a": (tmp_path / "a.fasta", tmp_path / "a.zip"), "b": (tmp_path / "b.fasta", tmp_path / "b.zip")},
                        kmer_table=str(tmp_path / "kmers"), genome_lineages_path=tmp_path / "genome-lineages.tsv",
                        labels_out=tmp_path / "kmer-labels.tsv", genes=genes, genes_out=tmp_path / "to-search.faa", all_genes=True)
    assert written == 3
    label_inputs(inputs={"a": tmp_path / "a.fasta", "b": tmp_path / "b.fasta"}, kmer_labels_path=tmp_path / "kmer-labels.tsv",
                 hits_path=tmp_path / "hits.tsv", protein_lineages_path=tmp_path / "protein-lineages.tsv",
                 out_paths={"a": tmp_path / "a-tax.tsv", "b": tmp_path / "b-tax.tsv"},
                 genes=genes, gene_out_paths={"a": tmp_path / "a-gene-tax.tsv", "b": tmp_path / "b-gene-tax.tsv"})

    genes_a = [line.split("\t") for line in (tmp_path / "a-gene-tax.tsv").read_text().splitlines()]
    assert genes_a[0] == ["gene_id", "domain", "phylum", "class", "order", "family", "genus", "species"]
    assert genes_a[1] == ["c1_1"] + ["NA"] * 7                    # called, but no protein hit
    assert genes_a[2] == ["c2_1", *SALMO[:6], "NA"]
    assert (tmp_path / "b-gene-tax.tsv").read_text().splitlines()[1].split("\t") == ["c1_1", *SALMO]

    a = [line.split("\t") for line in (tmp_path / "a-tax.tsv").read_text().splitlines()]
    b = [line.split("\t") for line in (tmp_path / "b-tax.tsv").read_text().splitlines()]
    assert a[1] == ["c1", *ECOLI, "k-mers"]
    assert a[2] == ["c2", *SALMO[:6], "NA", "proteins"]
    assert b[1] == ["c1", *SALMO, "proteins"]


def test_label_inputs_can_label_each_input_as_a_whole(tmp_path):
    """ e.g. a genome: its contigs' k-mers pooled and all its genes together, one row named by the input """
    import gzip, json, zipfile
    from bit.modules.kmer_lookup import write_table
    from bit.modules.classify_seqs import kmer_step, label_inputs

    write_table(str(tmp_path / "kmers"), [("gE", list(range(100, 150)))])
    (tmp_path / "genome-lineages.tsv").write_text("gE\t" + ";".join(ECOLI) + "\n")
    (tmp_path / "protein-lineages.tsv").write_text("7\t" + ";".join(SALMO) + "\n")

    def sketch(path, seqs):
        records = [{"name": name, "signatures": [{"mins": mins}]} for name, mins in seqs]
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(records).encode()))

    # g1: its contigs' k-mers pooled (and a k-mer in both counted once): 30 of 45 (67%) in gE
    (tmp_path / "g1.fasta").write_text(">c1\nACGT\n>c2\nACGT\n")
    sketch(tmp_path / "g1.zip", [("c1", list(range(100, 130)) + list(range(1000, 1005))), ("c2", [100] + list(range(3000, 3010)))])
    # g2: no k-mers shared; its genes (named however the gene caller names them) all count towards it
    (tmp_path / "g2.fasta").write_text(">x_1\nACGT\n>y\nACGT\n")
    sketch(tmp_path / "g2.zip", [("x_1", [1, 2]), ("y", [3])])
    (tmp_path / "hits.tsv").write_text("g2__x_1|a_b_1\t7.1\t99.0\t300\ng2__y|c_2\t7.2\t99.0\t300\n")
    (tmp_path / "g1.faa").write_text(">g1__c1_1\nMK\n>g1__c2_1\nMK\n")
    (tmp_path / "g2.faa").write_text(">g2__x_1|a_b_1\nMK\n>g2__y|c_2\nMK\n")

    # only g2 is left for genes to be called (with genes called after the k-mers)
    kmer_step(inputs={"g1": (tmp_path / "g1.fasta", tmp_path / "g1.zip"), "g2": (tmp_path / "g2.fasta", tmp_path / "g2.zip")},
              kmer_table=str(tmp_path / "kmers"), genome_lineages_path=tmp_path / "genome-lineages.tsv",
              labels_out=tmp_path / "kmer-labels.tsv", undecided_out=tmp_path / "undecided.tsv", whole_inputs=True)
    assert (tmp_path / "undecided.tsv").read_text() == f"g2\t{tmp_path / 'g2.fasta'}\n"

    # only g2's genes need searching: g1 is labelled by its k-mers
    kmer_step(inputs={"g1": (tmp_path / "g1.fasta", tmp_path / "g1.zip"), "g2": (tmp_path / "g2.fasta", tmp_path / "g2.zip")},
              kmer_table=str(tmp_path / "kmers"), genome_lineages_path=tmp_path / "genome-lineages.tsv",
              labels_out=tmp_path / "kmer-labels.tsv", genes={"g1": tmp_path / "g1.faa", "g2": tmp_path / "g2.faa"},
              genes_out=tmp_path / "to-search.faa", whole_inputs=True)
    assert [l for l in (tmp_path / "to-search.faa").read_text().splitlines() if l.startswith(">")] == [">g2__x_1|a_b_1", ">g2__y|c_2"]
    label_inputs(inputs={"g1": tmp_path / "g1.fasta", "g2": tmp_path / "g2.fasta"}, kmer_labels_path=tmp_path / "kmer-labels.tsv",
                 hits_path=tmp_path / "hits.tsv", protein_lineages_path=tmp_path / "protein-lineages.tsv",
                 out_paths={"g1": tmp_path / "g1-tax.tsv", "g2": tmp_path / "g2-tax.tsv"}, whole_inputs=True)

    g1 = [line.split("\t") for line in (tmp_path / "g1-tax.tsv").read_text().splitlines()]
    g2 = [line.split("\t") for line in (tmp_path / "g2-tax.tsv").read_text().splitlines()]
    assert g1[1:] == [["g1", *ECOLI, "k-mers"]]
    assert g2[1:] == [["g2", *SALMO, "proteins"]]


def test_only_genes_of_sequences_without_a_kmer_label_are_searched(tmp_path):
    import gzip, json, zipfile
    from bit.modules.kmer_lookup import write_table
    from bit.modules.classify_seqs import kmer_step, read_kmer_labels

    write_table(str(tmp_path / "kmers"), [("gE", list(range(100, 150)))])
    (tmp_path / "genome-lineages.tsv").write_text("gE\t" + ";".join(ECOLI) + "\n")
    records = [{"name": "c1", "signatures": [{"mins": list(range(100, 140))}]}, {"name": "c2", "signatures": [{"mins": [1, 2, 3]}]}]
    with zipfile.ZipFile(tmp_path / "a.zip", "w") as archive:
        archive.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(records).encode()))
    (tmp_path / "a.fasta").write_text(">c1\nACGT\n>c2\nACGT\n")
    (tmp_path / "a.faa").write_text(">a__c1_1\nMK\n>a__c1_2\nMKL\n>a__c2_1\nMV\n")

    written = kmer_step(inputs={"a": (tmp_path / "a.fasta", tmp_path / "a.zip")}, kmer_table=str(tmp_path / "kmers"),
                        genome_lineages_path=tmp_path / "genome-lineages.tsv", labels_out=tmp_path / "labels.tsv",
                        genes={"a": tmp_path / "a.faa"}, genes_out=tmp_path / "to-search.faa")
    assert written == 1 and (tmp_path / "to-search.faa").read_text() == ">a__c2_1\nMV\n"
    assert read_kmer_labels(tmp_path / "labels.tsv") == {("a", "c1"): ECOLI, ("a", "c2"): ()}


def test_find_fasta_inputs(tmp_path):
    import pytest # type: ignore
    from bit.modules.classify_seqs import find_fasta_inputs

    d = tmp_path / "dir"
    d.mkdir()
    for name in ("a.fasta", "b.fna", "notes.txt"):
        (d / name).write_text(">s\nACGT\n")
    single = tmp_path / "c.assembly.fa"
    single.write_text(">s\nACGT\n")

    assert find_fasta_inputs([str(single)]) == {"c.assembly": str(single)}
    assert find_fasta_inputs([str(d), str(single)]) == {"a": str(d / "a.fasta"), "b": str(d / "b.fna"), "c.assembly": str(single)}
    # a file given directly is taken whatever its extension
    odd = tmp_path / "contigs.txt"
    odd.write_text(">s\nACGT\n")
    assert find_fasta_inputs([str(odd)]) == {"contigs": str(odd)}

    empty = tmp_path / "empty"
    empty.mkdir()
    for bad, message in (([str(empty)], "No fasta files"), ([str(d), str(d / "a.fasta")], "same output table"),
                         ([str(tmp_path / "x.fasta.gz")], "gzipped"),
                         ([str(tmp_path / "a.fasta"), str(tmp_path / "a__b.fasta")], "can't be told apart")):
        with pytest.raises(ValueError, match=message):
            find_fasta_inputs(bad)
    # a name holding the separator is fine on its own
    assert find_fasta_inputs([str(tmp_path / "a__b.fasta")]) == {"a__b": str(tmp_path / "a__b.fasta")}


def test_gene_names_split_by_the_input_names_even_when_they_hold_the_separator():
    import pytest # type: ignore
    from bit.modules.classify_seqs import input_name_splitter
    split = input_name_splitter(["C__x", "plain", "s1"])
    assert split("C__x__c_C__x_1_1") == ("C__x", "c_C__x_1_1")
    assert split("plain__contig__7_2") == ("plain", "contig__7_2")
    assert split("s1__k141_5_3") == ("s1", "k141_5_3")
    with pytest.raises(ValueError, match="can't be told apart"):
        input_name_splitter(["a", "a__b"])
