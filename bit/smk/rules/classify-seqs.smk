# Taxonomic classification of sequences (contigs, or high-quality long reads), shared by `bit classify-seqs`, the
# metagenomics workflow and characterize-genomes (--euk). Each sequence gets a species label from its sampled k-mers
# when a reference genome holds enough of them, and otherwise the deepest label its genes' protein hits agree on
# (bit.modules.classify_seqs).
#
# Pulled in with `include:` after rules/ref-dbs.smk and after these are defined:
#   tax_IDs                the names of the inputs (samples, fasta files, or genomes)
#   tax_seqs_path(ID)      an input's sequences (fasta)
#   tax_genes_path(ID)     its called genes (amino acids)
#   tax_tables             {ID: the table of its sequences' labels to write}
#   tax_gene_tables        {ID: the table of its genes' labels to write}, or None for no gene tables
#   tax_whole_inputs       True to label each input as a whole (e.g. a genome: all its k-mers and genes together,
#                          one row named by the input), False to label each of its sequences
#   tax_genes_after_kmers  (with tax_whole_inputs) True to call genes (prodigal) only for the inputs their k-mers
#                          didn't label, rather than taking them from tax_genes_path for all of them
#   tax_work_dir           where the intermediate files go
#   tax_intermediate(path) how intermediate files are kept (e.g. temp())
#   logs_dir, benchmarks_dir
# and config holding: threads, identity_cutoffs, min_containment, genome_min_containment, search_mem,
# classify_seqs_mem

# the reference (https://github.com/AstrobioMike/bit-tax-db, `bit data get tax-db`): proteins to search (each named
# <ID of its lineage>.<number>), the table of the reference genomes' sampled k-mers (see bit.modules.kmer_lookup), and
# the lineages of both
from bit.modules.classify_seqs import DEFAULT_TOP_PERCENT

tax_db_dir = ref_db_dir("tax-db")
tax_proteins_db = tax_db_dir + "proteins.dmnd"
tax_protein_lineages = tax_db_dir + "protein-lineages.tsv"
tax_kmer_table = tax_db_dir + "genome-kmers"
tax_genome_lineages = tax_db_dir + "genome-lineages.tsv"

# the k-mers sampled from the sequences have to be sampled the way the reference genomes' were
TAX_KMER_SIZE = 31
TAX_KMER_SCALED = 200


# genes named <input name>__<gene ID>, so all inputs' genes can be searched together
rule tax_name_genes:
    input:
        lambda wildcards: tax_genes_path(wildcards.ID)
    output:
        tax_intermediate(tax_work_dir + "{ID}-genes-for-search.faa")
    shell:
        r"""
        sed 's/^>\([^ ]*\).*/>{wildcards.ID}__\1/; s/\*$//' {input} > {output}
        """


# one search for all the inputs (loading the reference is most of the cost), of only the genes of sequences (or, as a
# whole, inputs) that their k-mers didn't already give a species; when every gene's own label is wanted (gene tables),
# of all of them. No search at all when there's nothing left to search.
rule tax_search_proteins:
    conda:
        env("classify-seqs.yaml")
    input:
        genes = tax_work_dir + "tax-genes-to-search.faa"
    output:
        tax_intermediate(tax_work_dir + "tax-protein-hits.tsv")
    params:
        db = tax_proteins_db,
        top = DEFAULT_TOP_PERCENT
    threads:
        config["threads"]
    resources:
        mem_mb = config["search_mem"],
        cpus = config["threads"]
    log:
        logs_dir + "tax-search-proteins.log"
    benchmark:
        benchmarks_dir + "tax_search_proteins-benchmarks.tsv"
    shell:
        """
        if [ -s {input.genes} ]; then
            diamond blastp --db "{params.db}" --query {input.genes} --out {output} --outfmt 6 qseqid sseqid pident bitscore \
                --top {params.top} --evalue 1e-5 --threads {threads} > {log} 2>&1
        else
            touch {output}
            printf "Every sequence got its taxonomic assignment from k-mers, so no protein search was needed.\n" > {log}
        fi
        """


# an input with no sequences gets an empty file in place of the sketch
rule tax_sketch_seqs:
    conda:
        env("classify-seqs.yaml")
    input:
        lambda wildcards: tax_seqs_path(wildcards.ID)
    output:
        tax_intermediate(tax_work_dir + "{ID}.sig.zip")
    params:
        sketch = f"k={TAX_KMER_SIZE},scaled={TAX_KMER_SCALED}"
    log:
        logs_dir + "{ID}-tax-sketch.log"
    benchmark:
        benchmarks_dir + "tax_sketch_seqs-{ID}-benchmarks.tsv"
    shell:
        """
        if [ -s "{input}" ]; then
            sourmash sketch dna -p {params.sketch} --singleton "{input}" -o {output} > {log} 2>&1
        else
            touch {output}
        fi
        """


# the k-mer lookup (one pass over the reference table for all the inputs), and which genes are left to search (or, with
# tax_genes_after_kmers, which inputs are left to call genes for)
rule tax_kmer_step:
    input:
        seqs = [tax_seqs_path(ID) for ID in tax_IDs],
        sketches = expand(tax_work_dir + "{ID}.sig.zip", ID = tax_IDs),
        genes = [] if tax_genes_after_kmers else expand(tax_work_dir + "{ID}-genes-for-search.faa", ID = tax_IDs)
    output:
        labels = tax_intermediate(tax_work_dir + "tax-kmer-labels.tsv"),
        genes = tax_intermediate(tax_work_dir + ("tax-undecided.tsv" if tax_genes_after_kmers else "tax-genes-to-search.faa"))
    threads:
        config["threads"]
    resources:
        mem_mb = config["classify_seqs_mem"],
        cpus = config["threads"]
    benchmark:
        benchmarks_dir + "tax_kmer_step-benchmarks.tsv"
    run:
        from bit.modules.classify_seqs import kmer_step
        kmer_step(
            inputs = {ID: (seqs, sketch) for ID, seqs, sketch in zip(tax_IDs, input.seqs, input.sketches)},
            kmer_table = tax_kmer_table, genome_lineages_path = tax_genome_lineages, labels_out = output.labels,
            genes = None if tax_genes_after_kmers else dict(zip(tax_IDs, input.genes)),
            genes_out = None if tax_genes_after_kmers else output.genes, all_genes = bool(tax_gene_tables),
            undecided_out = output.genes if tax_genes_after_kmers else None, threads = threads,
            min_containment = config["genome_min_containment"] if tax_whole_inputs else config["min_containment"],
            whole_inputs = tax_whole_inputs,
        )


# genes (prodigal, named <input>__<gene>) of only the inputs the k-mers didn't label, several at a time
if tax_genes_after_kmers:

    rule tax_genes_of_undecided:
        conda:
            env("classify-seqs.yaml")
        input:
            tax_work_dir + "tax-undecided.tsv"
        params:
            tmp_dir = tax_work_dir + "undecided-genes-tmp/"
        output:
            tax_intermediate(tax_work_dir + "tax-genes-to-search.faa")
        threads:
            config["threads"]
        resources:
            cpus = config["threads"]
        log:
            logs_dir + "tax-genes-of-undecided.log"
        benchmark:
            benchmarks_dir + "tax_genes_of_undecided-benchmarks.tsv"
        shell:
            r"""
            rm -rf {params.tmp_dir} && mkdir -p {params.tmp_dir}
            : > {log}
            tr '\t\n' '\000\000' < {input} | xargs -0 -r -n 2 -P {threads} sh -c                 'prodigal -p meta -q -i "$1" -a {params.tmp_dir}"$0".faa -o /dev/null 2>> {log}'
            : > {output}
            cut -f 1 {input} | while read -r name; do
                sed "s/^>\([^ ]*\).*/>${{name}}__\1/; s/\*$//" {params.tmp_dir}"$name".faa >> {output}
            done
            rm -rf {params.tmp_dir}
            """


# each sequence's label: from its k-mers, or else from its genes' protein hits
rule tax_label_seqs:
    input:
        seqs = [tax_seqs_path(ID) for ID in tax_IDs],
        genes = expand(tax_work_dir + "{ID}-genes-for-search.faa", ID = tax_IDs) if tax_gene_tables else [],
        kmer_labels = tax_work_dir + "tax-kmer-labels.tsv",
        hits = tax_work_dir + "tax-protein-hits.tsv"
    output:
        tables = [tax_tables[ID] for ID in tax_IDs],
        gene_tables = [tax_gene_tables[ID] for ID in tax_IDs] if tax_gene_tables else []
    resources:
        mem_mb = config["classify_seqs_mem"]
    benchmark:
        benchmarks_dir + "tax_label_seqs-benchmarks.tsv"
    run:
        from bit.modules.classify_seqs import label_inputs
        label_inputs(
            inputs = dict(zip(tax_IDs, input.seqs)), kmer_labels_path = input.kmer_labels,
            hits_path = input.hits, protein_lineages_path = tax_protein_lineages, out_paths = tax_tables,
            cutoffs = tuple(float(c) for c in str(config["identity_cutoffs"]).split(",")),
            genes = dict(zip(tax_IDs, input.genes)) if tax_gene_tables else None,
            gene_out_paths = tax_gene_tables,
            whole_inputs = tax_whole_inputs,
        )
