"""
Taxonomic labelling rules for `bit classify-seqs` (contigs or high-quality long reads).

A lineage is a tuple of names from domain down (domain, phylum, class, order, family, genus, species),
cut off wherever the label stops; the empty tuple means no label.

Two sources of evidence per sequence:
  - DNA k-mers shared with reference genomes: used only for a species call, and only when one
    genome (or genomes of one species) holds a large share of the sequence's sampled k-mers
  - protein hits of its genes: each gene gets the lowest common ancestor of its near-best hits, capped
    by how similar the best hit is; the sequence gets the deepest taxon that holds more than half of
    the summed gene scores
"""

import os

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")

DEFAULT_TOP_PERCENT = 10
DEFAULT_IDENTITY_CUTOFFS = (95, 80, 60)   # least % identity of a gene's best hit for a species, genus, family label
DEFAULT_MIN_SHARED_KMERS = 3
DEFAULT_MIN_CONTAINMENT = 0.65
# for a whole genome (--genomes, characterize-genomes, MAGs where GTDB-Tk doesn't run): its thousands of sampled k-mers
# make a lower share a solid signal. 0.40 was the lowest that made no wrong species call for 46 benchmark genomes with
# their own species hidden (sister species reached 0.37), while calling 45 of them right with it present (the lowest
# same-species value was 0.38); 0.20 (~95% identity) mixed up close species in 11 of the 46
DEFAULT_GENOME_MIN_CONTAINMENT = 0.40
DEFAULT_KMER_WINDOW = 0.9


def common_lineage(lineages):
    """ the shared start of the given lineages """

    lineages = list(lineages)
    if not lineages:
        return ()

    shared = []
    for names in zip(*lineages):
        if len(set(names)) != 1:
            break
        shared.append(names[0])

    return tuple(shared)


UNNAMED_PREFIX = "unclassified "   # the reference fills a rank its taxonomy skips as "unclassified <taxon above> <rank>"


def without_unnamed_tail(lineage):
    """
    drops trailing filled ranks: every taxon under the same named ancestor carries the same fill, so a shared
    fill is not evidence the taxa are related below that ancestor (two viruses with only a species name share fills down to genus)
    """

    lineage = tuple(lineage)
    while lineage and lineage[-1].startswith(UNNAMED_PREFIX):
        lineage = lineage[:-1]

    return lineage


def deepest_rank_allowed(percent_identity, cutoffs=DEFAULT_IDENTITY_CUTOFFS):
    """ how many ranks (counted from domain) a protein hit of this identity may support """

    for depth, cutoff in zip((7, 6, 5), cutoffs):
        if percent_identity >= cutoff:
            return depth

    return 4


def label_gene(hits, top_percent=DEFAULT_TOP_PERCENT, cutoffs=DEFAULT_IDENTITY_CUTOFFS):
    """
    hits: (lineage, bitscore, percent identity) per protein hit of one gene; an empty lineage is a protein
    of a cluster whose members span domains (labelled root), so a near-best hit to one leaves the gene unlabelled
    returns (lineage, best bitscore); the lineage is empty when there are no hits or no shared taxon
    """

    if not hits:
        return (), 0.0

    best_score = max(hit[1] for hit in hits)
    near_best = [hit for hit in hits if hit[1] >= best_score * (1 - top_percent / 100)]
    lineage = common_lineage(hit[0] for hit in near_best)
    depth = deepest_rank_allowed(max(hit[2] for hit in near_best), cutoffs)

    return without_unnamed_tail(lineage[:depth]), best_score


def label_from_genes(gene_labels):
    """
    gene_labels: (lineage, bitscore) per gene of one sequence
    returns the deepest lineage that holds more than half of the summed scores of the labelled genes
    """

    gene_labels = [(lineage, score) for lineage, score in gene_labels if lineage]
    total = sum(score for _, score in gene_labels)
    if total <= 0:
        return ()

    weights = {}
    for lineage, score in gene_labels:
        for depth in range(1, len(lineage) + 1):
            weights[lineage[:depth]] = weights.get(lineage[:depth], 0.0) + score

    supported = [lineage for lineage, weight in weights.items() if weight > total / 2]
    if not supported:
        return ()

    return without_unnamed_tail(max(supported, key=len))


def label_from_kmers(matches, num_kmers, min_shared=DEFAULT_MIN_SHARED_KMERS,
                     min_containment=DEFAULT_MIN_CONTAINMENT, window=DEFAULT_KMER_WINDOW):
    """
    matches: (lineage, shared k-mer count) per reference genome sharing k-mers with the sequence
    num_kmers: the number of sampled k-mers in the sequence
    returns a species lineage, or the empty tuple when the k-mers don't support a species call
    """

    matches = [match for match in matches if match[0]]
    if not matches or num_kmers <= 0:
        return ()

    best = max(shared for _, shared in matches)
    if best < min_shared or best / num_kmers < min_containment:
        return ()

    lineage = common_lineage(lineage for lineage, shared in matches if shared >= best * window)

    return lineage if len(lineage) == len(RANKS) else ()


def classify(kmer_lineage, gene_lineage):
    """ returns (lineage, the evidence it came from: "k-mers", "proteins" or "none") """

    if kmer_lineage:
        return kmer_lineage, "k-mers"
    if gene_lineage:
        return gene_lineage, "proteins"

    return (), "none"


############################################################
##################### reading and writing ##################
############################################################

DIAMOND_COLUMNS = ("qseqid", "sseqid", "pident", "bitscore")   # the --outfmt 6 fields the search is run with
TABLE_HEADER = ("seq_id",) + RANKS + ("evidence",)
NO_LABEL = "NA"


def read_lineages(path):
    """
    A two-column table: reference ID (a genome of the k-mer table, or the part of a protein's name before
    its last "."), then its lineage as names from domain down joined by ";". Returns {ID: lineage tuple}.
    """

    lineages = {}
    with open(path) as table:
        for line in table:
            ref_id, _, lineage = line.rstrip("\n").partition("\t")
            lineages[ref_id] = tuple(name for name in lineage.split(";") if name)

    return lineages


def seq_of_gene(gene_id):
    """
    prodigal names the genes of a sequence <sequence ID>_<number>
    (a gene given as (input name, gene ID) gives (input name, sequence ID))
    """

    if isinstance(gene_id, tuple):
        return (gene_id[0], gene_id[1].rsplit("_", 1)[0])
    return gene_id.rsplit("_", 1)[0]


def read_protein_hits(path, protein_lineages):
    """
    Reads DIAMOND output (fields as in DIAMOND_COLUMNS, a gene's hits on consecutive lines) and yields
    (gene ID, [(lineage, bitscore, percent identity), ...]) per gene. A reference protein is named
    <its lineage's ID>.<number>; hits to proteins whose ID isn't in the lineage table are left out (root, with an
    empty lineage, is in it).
    """

    current, hits = None, []
    with open(path) as table:
        for line in table:
            gene, protein, identity, bitscore = line.rstrip("\n").split("\t")
            if gene != current:
                if current is not None:
                    yield current, hits
                current, hits = gene, []
            lineage = protein_lineages.get(protein.rsplit(".", 1)[0])
            if lineage is not None:
                hits.append((lineage, float(bitscore), float(identity)))
    if current is not None:
        yield current, hits


def kmer_labels(seq_ids, kmer_matches, kmer_counts, genome_lineages,
                min_shared=DEFAULT_MIN_SHARED_KMERS, min_containment=DEFAULT_MIN_CONTAINMENT):
    """
    kmer_matches:  (sequence ID, genome ID, shared k-mer count)
    kmer_counts:   {sequence ID: its number of sampled k-mers}
    Returns {sequence ID: species lineage from its k-mers, or () where they don't support one}.
    """

    genomes = {}
    for seq, genome, shared in kmer_matches:
        genomes.setdefault(seq, []).append((genome_lineages.get(genome, ()), shared))

    return {seq: label_from_kmers(genomes.get(seq, []), kmer_counts.get(seq, 0), min_shared, min_containment)
            for seq in seq_ids}


def combine_evidence(seq_ids, protein_hits, from_kmers, top_percent=DEFAULT_TOP_PERCENT,
                     cutoffs=DEFAULT_IDENTITY_CUTOFFS, gene_labels=None, seq_of=seq_of_gene):
    """
    from_kmers:    {sequence ID: species lineage from k-mers, or ()} (see kmer_labels)
    protein_hits:  (gene ID, hits) as from read_protein_hits (only needed for the sequences without a k-mer label)
    gene_labels:   when a dict is given, it is filled with {gene ID: lineage} for the genes with protein hits
    seq_of:        gene ID -> the ID of what it counts towards (by default its sequence, see seq_of_gene)
    Returns {sequence ID: (lineage, evidence)} for every sequence given.
    """

    genes = {}
    for gene, hits in protein_hits:
        label = label_gene(hits, top_percent, cutoffs)
        genes.setdefault(seq_of(gene), []).append(label)
        if gene_labels is not None:
            gene_labels[gene] = label[0]

    labels = {}
    for seq in seq_ids:
        kmer_label = from_kmers.get(seq, ())
        labels[seq] = classify(kmer_label, () if kmer_label else label_from_genes(genes.get(seq, [])))

    return labels


def classify_seqs(seq_ids, protein_hits, kmer_matches, kmer_counts, genome_lineages,
                  top_percent=DEFAULT_TOP_PERCENT, cutoffs=DEFAULT_IDENTITY_CUTOFFS,
                  min_shared=DEFAULT_MIN_SHARED_KMERS, min_containment=DEFAULT_MIN_CONTAINMENT,
                  gene_labels=None, seq_of=seq_of_gene):
    """ kmer_labels() then combine_evidence(), for every sequence given: {sequence ID: (lineage, evidence)} """

    from_kmers = kmer_labels(seq_ids, kmer_matches, kmer_counts, genome_lineages, min_shared, min_containment)
    return combine_evidence(seq_ids, protein_hits, from_kmers, top_percent, cutoffs, gene_labels, seq_of)


def write_labels(path, seq_ids, labels):
    """ one row per sequence: its ID, a column per rank (NA below where the label stops), and the evidence used """

    with open(path, "w") as out:
        out.write("\t".join(TABLE_HEADER) + "\n")
        for seq in seq_ids:
            lineage, evidence = labels[seq]
            names = list(lineage) + [NO_LABEL] * (len(RANKS) - len(lineage))
            out.write("\t".join([seq, *names, evidence]) + "\n")


def write_gene_labels(path, gene_ids, gene_labels):
    """ one row per gene: its ID and a column per rank (all NA for a gene without usable protein hits) """

    with open(path, "w") as out:
        out.write("\t".join(("gene_id",) + RANKS) + "\n")
        for gene in gene_ids:
            lineage = gene_labels.get(gene, ())
            out.write("\t".join([gene, *lineage, *[NO_LABEL] * (len(RANKS) - len(lineage))]) + "\n")


def read_seq_ids(fasta_path):
    """ the sequence IDs of a fasta file, in order """

    with open(fasta_path) as fasta:
        return [line[1:].split()[0] for line in fasta if line.startswith(">")]


INPUT_SEPARATOR = "__"   # genes searched together are named <input name>__<gene ID>


def input_name_splitter(names):
    """
    A function taking a searched gene's name, <input name>__<gene ID>, to (input name, gene ID), for inputs whose
    names may themselves hold the separator (e.g. a sample named C__x). Raises ValueError when two input names would
    make that ambiguous (one being another's name followed by the separator).
    """

    names = set(names)
    for name in names:
        clash = sorted(other for other in names if other != name and other.startswith(name + INPUT_SEPARATOR))
        if clash:
            raise ValueError(f"The input names '{name}' and '{clash[0]}' can't be told apart in the gene names "
                             f"(<input name>{INPUT_SEPARATOR}<gene>); please rename one of them.")

    def split(gene):
        start = 0
        while (at := gene.find(INPUT_SEPARATOR, start)) != -1:
            if gene[:at] in names:
                return gene[:at], gene[at + len(INPUT_SEPARATOR):]
            start = at + 1
        name, _, gene_id = gene.partition(INPUT_SEPARATOR)
        return name, gene_id

    return split


def kmer_step(inputs, kmer_table, genome_lineages_path, labels_out, genes=None, genes_out=None, all_genes=False,
              threads=1, min_shared=DEFAULT_MIN_SHARED_KMERS, min_containment=DEFAULT_MIN_CONTAINMENT, whole_inputs=False,
              undecided_out=None):
    """
    The first step of labelling (`bit classify-seqs`, and the taxonomy of the metagenomics workflow and
    characterize-genomes): the k-mer lookup, all the inputs in one pass over the reference table, and which genes still
    need searching against the reference proteins.
    inputs:          {input name: (fasta path, sourmash sketch zip of its sequences)}
    labels_out:      written: input name, sequence ID, and its species lineage from k-mers (";"-joined, empty for none);
                     with whole_inputs, one row per input (its name as the sequence ID)
    genes, genes_out: {input name: fasta of its genes, named <input name>__<gene ID>}, and the fasta to write of the genes
                     to search: those of the sequences (or inputs) without a k-mer label, or every gene with all_genes
                     (when each gene's own label is wanted, as for the metagenomics workflow's gene tables)
    whole_inputs:    label each input as a whole (e.g. a genome): its sequences' sampled k-mers pooled
    undecided_out:   with whole_inputs, written: the name and fasta path of each input the k-mers didn't label (for
                     calling genes on just those)
    Returns the number of genes written to genes_out.
    """

    import numpy as np
    from bit.modules.kmer_lookup import read_sketches, find_shared_kmers

    owners, hashes = [], []
    for name, (_, sketch_zip) in inputs.items():
        if os.path.getsize(sketch_zip) == 0:   # an input without sequences
            continue
        seq_names, seq_hashes = read_sketches(sketch_zip)
        if whole_inputs:
            owners.append((name, name))
            hashes.append(np.unique(np.concatenate(seq_hashes)) if seq_hashes else np.empty(0, dtype=np.uint64))
        else:
            owners += [(name, seq_name.split()[0]) for seq_name in seq_names]
            hashes += seq_hashes

    counts = {owner: len(h) for owner, h in zip(owners, hashes)}
    matches = [(owners[query], genome, shared)
               for query, genome, shared in find_shared_kmers(kmer_table, hashes, threads=threads)]
    units = [(name, name) for name in inputs] if whole_inputs else \
            [(name, seq) for name, (fasta, _) in inputs.items() for seq in read_seq_ids(fasta)]
    labels = kmer_labels(units, matches, counts, read_lineages(genome_lineages_path), min_shared, min_containment)

    with open(labels_out, "w") as out:
        for (name, seq), lineage in labels.items():
            out.write(f"{name}\t{seq}\t{';'.join(lineage)}\n")

    if undecided_out is not None:
        with open(undecided_out, "w") as out:
            for name, (fasta, _) in inputs.items():
                if not labels.get((name, name)):
                    out.write(f"{name}\t{fasta}\n")

    written = 0
    if genes_out is not None:
        with open(genes_out, "w") as out:
            for name, path in (genes or {}).items():
                keep = True
                with open(path) as fasta:
                    for line in fasta:
                        if line.startswith(">"):
                            gene_id = line[1:].split()[0][len(name) + len(INPUT_SEPARATOR):]
                            unit = (name, name) if whole_inputs else (name, seq_of_gene(gene_id))
                            keep = all_genes or not labels.get(unit)
                            written += keep
                        if keep:
                            out.write(line)
    return written


def read_kmer_labels(path):
    """ {(input name, sequence ID): lineage} as written by kmer_step """

    labels = {}
    with open(path) as table:
        for line in table:
            name, seq, lineage = line.rstrip("\n").split("\t")
            labels[(name, seq)] = tuple(n for n in lineage.split(";") if n)
    return labels


def label_inputs(inputs, kmer_labels_path, hits_path, protein_lineages_path, out_paths,
                 top_percent=DEFAULT_TOP_PERCENT, cutoffs=DEFAULT_IDENTITY_CUTOFFS,
                 genes=None, gene_out_paths=None, whole_inputs=False):
    """
    The second step (after kmer_step): each sequence's label, from its k-mers, or else from its genes' protein hits.
    inputs:          {input name: fasta path}
    kmer_labels_path: kmer_step's labels
    hits_path:       DIAMOND output for the genes that were searched (all the inputs together)
    out_paths:       {input name: table to write}
    genes, gene_out_paths: to also write a table of the genes' own labels per input, give
                     {input name: fasta of its genes as searched} and {input name: table to write}
    whole_inputs:    as given to kmer_step: one row per input, all its genes counting towards it
    """

    from_kmers = read_kmer_labels(kmer_labels_path)
    protein_lineages = read_lineages(protein_lineages_path)

    split_gene_name = input_name_splitter(inputs)

    def hits_by_seq():
        for gene, hits in read_protein_hits(hits_path, protein_lineages):
            yield split_gene_name(gene), hits

    seq_ids = {name: ([name] if whole_inputs else read_seq_ids(fasta)) for name, fasta in inputs.items()}
    all_seqs = [(name, seq) for name, ids in seq_ids.items() for seq in ids]

    gene_labels = {} if gene_out_paths else None
    labels = combine_evidence(all_seqs, hits_by_seq(), from_kmers, top_percent, cutoffs, gene_labels=gene_labels,
                              # as a whole, every gene of an input counts towards the input
                              seq_of=(lambda gene: (gene[0], gene[0])) if whole_inputs else seq_of_gene)

    for name, ids in seq_ids.items():
        write_labels(out_paths[name], ids, {seq: labels[(name, seq)] for seq in ids})

    if gene_out_paths:
        for name, genes_fasta in genes.items():
            # the genes are named <input name>__<gene ID> in the fasta that was searched
            prefix = name + INPUT_SEPARATOR
            gene_ids = [gene[len(prefix):] if gene.startswith(prefix) else gene for gene in read_seq_ids(genes_fasta)]
            write_gene_labels(gene_out_paths[name], gene_ids,
                              {gene: gene_labels[(name, gene)] for gene in gene_ids if (name, gene) in gene_labels})


FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".fsa")


def find_fasta_inputs(paths):
    """
    paths: fasta files and/or directories holding fasta files (those ending with one of FASTA_EXTENSIONS)
    Returns {input name: path}, the name being the file name without its extension (what the output
    table is named for). Raises ValueError with a message for the user when the inputs can't be used.
    """

    import os

    files = []
    for path in paths:
        if os.path.isdir(path):
            found = sorted(os.path.join(path, name) for name in os.listdir(path)
                           if name.endswith(FASTA_EXTENSIONS) and os.path.isfile(os.path.join(path, name)))
            if not found:
                raise ValueError(f"No fasta files (ending with {', '.join(FASTA_EXTENSIONS)}) were found in '{path}' :(")
            files += found
        else:
            files.append(path)

    inputs = {}
    for path in files:
        if path.endswith(".gz"):
            raise ValueError(f"The fasta files can't be gzipped, but this one seems to be: '{path}'")
        name = os.path.splitext(os.path.basename(path))[0]
        if name in inputs:
            raise ValueError(f"These two fasta files would give the same output table ('{name}-tax.tsv'): "
                             f"'{inputs[name]}' and '{path}'")
        inputs[name] = path

    input_name_splitter(inputs)   # (raises if two names can't be told apart in the searched genes' names)
    return inputs
