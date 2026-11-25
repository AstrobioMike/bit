# Taxonomy of a directory of genomes (or bins/MAGs) from bit-tax-db, each genome labelled as a whole, as with
# `bit classify-seqs --genomes`: from the k-mers of all its contigs together when a reference genome holds enough of them
# (--genome-min-containment), otherwise from all its genes' protein hits (bit.modules.classify_seqs). Genes are only
# called and searched for the genomes their k-mers didn't label. Used where GTDB-Tk doesn't run (macOS). The genomes are found when
# the rules run (bins and MAGs only exist once binning is done).
#
# Pulled in with `include:` after rules/ref-dbs.smk, and after the workflow has listed what it wants run:
#   tax_db_genome_runs   [dict(name, genomes_dir, extension, needs, taxonomy_tab, work_dir, threads, search_mem,
#                              classify_seqs_mem, log, benchmark), ...]
# and config holding: identity_cutoffs, genome_min_containment
# Each run writes taxonomy_tab: Assembly, then the seven ranks ("NA" below the label; just the header with no genomes).

from bit.modules.classify_seqs import DEFAULT_TOP_PERCENT

tax_db_dir = ref_db_dir("tax-db")


def genomes_in(directory, extension):
    import glob
    return {os.path.basename(path)[:-len(extension)]: path for path in sorted(glob.glob(directory + "*" + extension))}


for _run in tax_db_genome_runs:

    # each genome's k-mers, sampled the way the reference genomes' were
    rule:
        name:
            _run["name"] + "_sketches"
        conda:
            env("classify-seqs.yaml")
        input:
            _run["needs"]
        params:
            genomes_dir = _run["genomes_dir"],
            extension = _run["extension"],
            work_dir = _run["work_dir"]
        output:
            temp(touch(_run["work_dir"] + "sketches.done"))
        log:
            _run["log"]
        benchmark:
            _run["benchmark"] + "-sketches.tsv"
        shell:
            r"""
            mkdir -p {params.work_dir}
            rm -f {params.work_dir}*.sig.zip
            : > {log}
            for genome in "{params.genomes_dir}"*{params.extension}; do
                [ -e "$genome" ] || continue
                name=$(basename "$genome" {params.extension})
                sourmash sketch dna -p k=31,scaled=200 --singleton "$genome" -o {params.work_dir}"$name".sig.zip >> {log} 2>&1
            done
            """

    # the k-mer lookup, and the genomes it didn't label
    rule:
        name:
            _run["name"] + "_kmers"
        input:
            _run["work_dir"] + "sketches.done"
        params:
            genomes_dir = _run["genomes_dir"],
            extension = _run["extension"],
            work_dir = _run["work_dir"]
        output:
            labels = temp(_run["work_dir"] + "kmer-labels.tsv"),
            undecided = temp(_run["work_dir"] + "undecided.tsv")
        threads:
            _run["threads"]
        resources:
            cpus = _run["threads"],
            mem_mb = _run["classify_seqs_mem"]
        benchmark:
            _run["benchmark"] + "-kmers.tsv"
        run:
            from bit.modules.classify_seqs import kmer_step
            genomes = genomes_in(params.genomes_dir, params.extension)
            kmer_step(
                inputs = {name: (path, params.work_dir + name + ".sig.zip") for name, path in genomes.items()},
                kmer_table = tax_db_dir + "genome-kmers", genome_lineages_path = tax_db_dir + "genome-lineages.tsv",
                labels_out = output.labels, undecided_out = output.undecided, threads = threads,
                min_containment = config["genome_min_containment"], whole_inputs = True,
            )

    # genes (prodigal, named <genome>__<gene>) of only those genomes, several at a time
    rule:
        name:
            _run["name"] + "_genes"
        conda:
            env("classify-seqs.yaml")
        input:
            _run["work_dir"] + "undecided.tsv"
        params:
            tmp_dir = _run["work_dir"] + "genes-tmp/"
        output:
            temp(_run["work_dir"] + "genes-to-search.faa")
        threads:
            _run["threads"]
        resources:
            cpus = _run["threads"]
        log:
            _run["log"].replace(".log", "-genes.log")
        benchmark:
            _run["benchmark"] + "-genes.tsv"
        shell:
            r"""
            rm -rf {params.tmp_dir} && mkdir -p {params.tmp_dir}
            : > {log}
            tr '\t\n' '\000\000' < {input} | xargs -0 -r -n 2 -P {threads} sh -c 'prodigal -q -c -p meta -i "$1" -a {params.tmp_dir}"$0".faa > /dev/null 2>> {log}'
            : > {output}
            cut -f 1 {input} | while read -r name; do
                sed "s/^>\([^ ]*\).*/>${{name}}__\1/; s/\*$//" {params.tmp_dir}"$name".faa >> {output}
            done
            rm -rf {params.tmp_dir}
            """

    rule:
        name:
            _run["name"] + "_search"
        conda:
            env("classify-seqs.yaml")
        input:
            genes = _run["work_dir"] + "genes-to-search.faa"
        params:
            db = tax_db_dir + "proteins.dmnd",
            top = DEFAULT_TOP_PERCENT
        output:
            temp(_run["work_dir"] + "protein-hits.tsv")
        threads:
            _run["threads"]
        resources:
            cpus = _run["threads"],
            mem_mb = _run["search_mem"]
        log:
            _run["log"].replace(".log", "-search.log")
        benchmark:
            _run["benchmark"] + "-search.tsv"
        shell:
            r"""
            if [ -s {input.genes} ]; then
                diamond blastp --db "{params.db}" --query {input.genes} --out {output} --outfmt 6 qseqid sseqid pident bitscore \
                    --top {params.top} --evalue 1e-5 --threads {threads} > {log} 2>&1
            else
                touch {output}
                printf "Every genome got its taxonomic assignment from k-mers, so no protein search was needed.\n" > {log}
            fi
            """

    rule:
        name:
            _run["name"] + "_label"
        input:
            labels = _run["work_dir"] + "kmer-labels.tsv",
            hits = _run["work_dir"] + "protein-hits.tsv"
        params:
            genomes_dir = _run["genomes_dir"],
            extension = _run["extension"],
            work_dir = _run["work_dir"]
        output:
            _run["taxonomy_tab"]
        resources:
            mem_mb = _run["classify_seqs_mem"]
        benchmark:
            _run["benchmark"] + "-label.tsv"
        run:
            import glob
            from bit.modules.classify_seqs import label_inputs, RANKS
            genomes = genomes_in(params.genomes_dir, params.extension)
            tables = {name: params.work_dir + name + "-tax.tsv" for name in genomes}
            label_inputs(
                inputs = genomes, kmer_labels_path = input.labels, hits_path = input.hits,
                protein_lineages_path = tax_db_dir + "protein-lineages.tsv", out_paths = tables,
                cutoffs = tuple(float(c) for c in str(config["identity_cutoffs"]).split(",")),
                whole_inputs = True,
            )
            with open(output[0], "w") as out:
                out.write("\t".join(("Assembly",) + RANKS) + "\n")
                for name in genomes:
                    with open(tables[name]) as table:
                        next(table)
                        for line in table:
                            out.write("\t".join(line.rstrip("\n").split("\t")[:8]) + "\n")
            # the per-genome sketches and tables (the steps' own outputs are removed by snakemake)
            for path in glob.glob(params.work_dir + "*.sig.zip") + list(tables.values()):
                os.remove(path)
