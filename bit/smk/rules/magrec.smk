# Binning and MAG recovery, shared by the metagenomics workflow and magrec: each assembly is binned (metabat2) with the
# coverage of the samples mapped to it, the bins are estimated for quality (checkm2), those passing the cutoffs become
# MAGs, which get taxonomy (GTDB-Tk; bit-tax-db where it doesn't run), and the MAGs of all assemblies are dereplicated
# (dRep) with each sample's reads mapped to the dereplicated set for the MAG-by-sample tables.
#
# Pulled in with `include:` after rules/ref-dbs.smk and rules/mapping.smk (whose asm_of, asm_path, asm_names and reads_of
# it uses), and after these are defined: sample_ID_list, paired_end, prefix, mapping_dir, bins_dir, MAGs_dir, logs_dir,
# benchmarks_dir; and config holding: threads, min_completion, max_redundancy, derep_ani, MAG_min_read_pid,
# mapping_mem, checkm2_mem, gtdbtk_cpus, gtdbtk_scratch_dir, gtdbtk_mem, and (for bit-tax-db)
# identity_cutoffs, genome_min_containment, search_mem, classify_seqs_mem.
# Defines magrec_outputs (the tables the workflow should ask for), and binning_done.

shared_scripts_dir = os.path.join(SMK, "scripts")

binning_done = expand(mapping_dir + "{asm}-metabat-assembly-depth.tsv", asm = asm_names)
MAGs_taxonomy_tab = MAGs_dir + prefix + ("MAGs-gtdb-taxonomy.tsv" if GTDBTK_RUNS_HERE else "MAGs-taxonomy.tsv")

# the dereplicated MAGs' detection, percent of reads recruited, and mean coverage in each sample
MAG_by_sample_tables = [MAGs_dir + prefix + f"MAG-{what}-by-sample.tsv" for what in ("detection", "percent-of-reads", "coverage")]

magrec_outputs = [
    MAGs_dir + prefix + "MAGs-overview.tsv",
    bins_dir + prefix + "bins-overview.tsv",
    MAGs_dir + prefix + "dereplicated-MAGs.tsv",
    *MAG_by_sample_tables
]

for dir in (bins_dir, MAGs_dir):
    os.makedirs(dir, exist_ok = True)


rule generate_MAGs_overview_table:
    input:
        assembly_summaries = MAGs_dir + prefix + "MAG-assembly-summaries.tsv",
        checkm_results = MAGs_dir + prefix + "MAGs-checkm-out.tsv",
        taxonomy = MAGs_taxonomy_tab
    output:
        MAGs_dir + prefix + "MAGs-overview.tsv"
    benchmark:
        benchmarks_dir + "generate_MAGs_overview_table-benchmarks.tsv"
    run:
        # assembly stats, checkm2 estimates and taxonomy (GTDB-Tk's, or bit-tax-db's where GTDB-Tk doesn't run) in one
        # table, ordered by taxonomy (or the note that no MAGs were recovered)
        from bit.modules.genome_overview import write_overview
        if GTDBTK_RUNS_HERE:
            write_overview(input.assembly_summaries, output[0], checkm2_path = input.checkm_results, gtdbtk_paths = [input.taxonomy])
        else:
            write_overview(input.assembly_summaries, output[0], checkm2_path = input.checkm_results, rank_table_path = input.taxonomy)


# MAGs across samples: the MAGs of all samples are dereplicated, each sample's reads are mapped to the dereplicated set
# together (so a read goes to the MAG it matches best), and each MAG's detection, share of the reads, and coverage
# per sample are tabulated
rule dereplicate_MAGs:
    conda:
        env("drep.yaml")
    input:
        checkm_results = MAGs_dir + prefix + "MAGs-checkm-out.tsv"
    params:
        MAGs_dir = MAGs_dir,
        drep_dir = MAGs_dir + "drep-out.tmp",
        genome_info = MAGs_dir + "drep-genome-info.tmp",
        ani = config["derep_ani"] / 100
    output:
        table = MAGs_dir + prefix + "dereplicated-MAGs.tsv",
        MAGs_list = temp(mapping_dir + "dereplicated-MAGs-list.txt"),
        fasta = temp(mapping_dir + "dereplicated-MAGs.fasta")
    threads:
        config["threads"]
    resources:
        cpus = config["threads"]
    log:
        logs_dir + "dereplicate-MAGs.log"
    benchmark:
        benchmarks_dir + "dereplicate_MAGs-benchmarks.tsv"
    shell:
        r"""
        num_MAGs=$(find {params.MAGs_dir} -maxdepth 1 -name "*.fasta" | wc -l | sed 's/^ *//')

        if [ ${{num_MAGs}} -gt 1 ]; then

            # dRep gets checkm2's estimates (so it never runs CheckM), and keeps every MAG (they passed the workflow's cutoffs already)
            awk -F $'\t' ' BEGIN {{ print "genome,completeness,contamination" }} NR > 1 {{ print $1".fasta,"$2","$3 }} ' {input.checkm_results} > {params.genome_info}
            rm -rf {params.drep_dir}
            dRep dereplicate {params.drep_dir} -g {params.MAGs_dir}*.fasta --genomeInfo {params.genome_info} \
                -comp 0 -con 100 -l 1 -sa {params.ani} --S_algorithm skani -p {threads} > {log} 2>&1

            # each MAG and its representative
            printf "MAG_ID\trepresentative\n" > {output.table}
            python {shared_scripts_dir}/drep-clusters.py {params.drep_dir} >> {output.table}
            rm -rf {params.drep_dir} {params.genome_info}

        elif [ ${{num_MAGs}} -eq 1 ]; then

            MAG=$(basename {params.MAGs_dir}*.fasta .fasta)
            printf "MAG_ID\trepresentative\n${{MAG}}\t${{MAG}}\n" > {output.table}
            printf "Only one MAG was recovered, so there was nothing to dereplicate.\n" > {log}

        else

            printf "There were no MAGs recovered.\n" > {output.table}
            printf "There were no MAGs recovered, so dRep was not run.\n" > {log}

        fi

        # the representatives, for mapping to
        if [ ${{num_MAGs}} -gt 0 ]; then
            tail -n +2 {output.table} | cut -f 2 | sort -u | sed 's|^|{params.MAGs_dir}|; s|$|.fasta|' > {output.MAGs_list}
        else
            touch {output.MAGs_list}
        fi
        cat $(cat {output.MAGs_list}) /dev/null > {output.fasta}
        """


rule index_dereplicated_MAGs:
    conda:
        env("mapping.yaml")
    input:
        mapping_dir + "dereplicated-MAGs.fasta"
    params:
        index = mapping_dir + "dereplicated-MAGs-index/index"
    output:
        temp(directory(mapping_dir + "dereplicated-MAGs-index"))
    threads:
        config["threads"]
    resources:
        cpus = config["threads"],
        mem_mb = config["mapping_mem"]
    log:
        logs_dir + "index-dereplicated-MAGs.log"
    benchmark:
        benchmarks_dir + "index_dereplicated_MAGs-benchmarks.tsv"
    shell:
        r"""
        mkdir -p {output}
        if [ -s {input} ]; then
            bowtie2-build --threads {threads} {input} {params.index} > {log} 2>&1
        else
            printf "There were no MAGs recovered, so there was nothing to index.\n" > {log}
        fi
        """


rule map_to_dereplicated_MAGs:
    conda:
        env("mapping.yaml")
    input:
        index = mapping_dir + "dereplicated-MAGs-index",
        reads = lambda wildcards: reads_of(wildcards.ID)
    params:
        index = mapping_dir + "dereplicated-MAGs-index/index",
        reads_args = lambda wildcards, input: "-1 " + input.reads[0] + " -2 " + input.reads[1] if paired_end else "-U " + input.reads[0]
    output:
        bam = temp(mapping_dir + "{ID}-to-dereplicated-MAGs.bam"),
        mapping_info = temp(mapping_dir + "{ID}-to-dereplicated-MAGs-mapping-info.txt")
    threads:
        config["threads"]
    resources:
        cpus = config["threads"],
        mem_mb = config["mapping_mem"]
    log:
        logs_dir + "{ID}-mapping-to-dereplicated-MAGs.log"
    benchmark:
        benchmarks_dir + "map_to_dereplicated_MAGs-{ID}-benchmarks.tsv"
    shell:
        r"""
        if [ -n "$(ls {input.index})" ]; then
            bowtie2 --mm -q --threads {threads} -x {params.index} {params.reads_args} --no-unal 2> {output.mapping_info} | samtools sort -@ {threads} -o {output.bam} - 2> {log}
        else
            touch {output.bam} {output.mapping_info}
            printf "Mapping not performed because there were no MAGs recovered.\n" > {log}
        fi
        """


rule dereplicated_MAG_cov_stats:
    # (run with the bit running the workflow: --min-read-pid is newer than the bit in envs/bit.yaml)
    input:
        bam = mapping_dir + "{ID}-to-dereplicated-MAGs.bam",
        MAGs_list = mapping_dir + "dereplicated-MAGs-list.txt"
    params:
        out_prefix = mapping_dir + "{ID}-dereplicated-MAGs-cov-stats",
        min_read_pid = config["MAG_min_read_pid"]
    output:
        temp(mapping_dir + "{ID}-dereplicated-MAGs-cov-stats-per-ref.tsv")
    log:
        logs_dir + "{ID}-dereplicated-MAGs-cov-stats.log"
    benchmark:
        benchmarks_dir + "dereplicated_MAG_cov_stats-{ID}-benchmarks.tsv"
    shell:
        r"""
        if [ -s {input.bam} ]; then
            bit cov-stats -R {input.MAGs_list} -b {input.bam} -o {params.out_prefix} --skip-per-contig --min-read-pid {params.min_read_pid} -F > {log} 2>&1
        else
            touch {output}
            printf "There were no MAGs recovered, so cov-stats was not run.\n" > {log}
        fi
        """


rule make_MAG_by_sample_tables:
    input:
        cov_stats = expand(mapping_dir + "{ID}-dereplicated-MAGs-cov-stats-per-ref.tsv", ID = sample_ID_list),
        mapping_info = expand(mapping_dir + "{ID}-to-dereplicated-MAGs-mapping-info.txt", ID = sample_ID_list)
    output:
        MAG_by_sample_tables
    benchmark:
        benchmarks_dir + "make_MAG_by_sample_tables-benchmarks.tsv"
    run:
        from bit.modules.MAG_abundance import write_MAG_tables
        write_MAG_tables(sample_ID_list, input.cov_stats, input.mapping_info, paired_end, *output)


rule summarize_MAG_assemblies:
    conda:
        env("bit.yaml")
    input:
        trigger = MAGs_dir + prefix + "MAGs-checkm-out.tsv"
    params:
        MAGs_dir = MAGs_dir
    output:
        temp(MAGs_dir + prefix + "MAG-assembly-summaries.tsv")
    benchmark:
        benchmarks_dir + "summarize_MAG_assemblies-benchmarks.tsv"
    shell:
        r"""
        # only running if any MAGs were recovered
        if [ $(find {params.MAGs_dir} -name "*.fasta" | wc -l | sed 's/^ *//') -gt 0 ]; then

            bit summarize-assembly {params.MAGs_dir}*.fasta -o {output} -t

        else

            printf "There were no MAGs recovered.\n" > {output}

        fi
        """


# checkm2 on the bins, and GTDB-Tk on the MAGs only (the bins that pass the quality cutoffs): low-quality bins
# take GTDB-Tk the longest and what they are assigned is the least reliable (rules shared with `bit characterize-genomes`);
# where GTDB-Tk doesn't run (macOS), the MAGs' taxonomy comes from bit-tax-db (rules/genome-tax.smk)
checkm2_runs = [dict(
    name = "run_checkm2_on_bins", genomes_dir = bins_dir, extension = ".fasta", what = "bins",
    needs = binning_done,
    results_tab = bins_dir + prefix + "bins-checkm-out.tsv", out_dir = None, temp_outputs = True,
    threads = config["threads"], mem_mb = config["checkm2_mem"],
    log = logs_dir + "checkm2.log", benchmark = benchmarks_dir + "run_checkm2_on_bins-benchmarks.tsv")]

gtdbtk_runs = [dict(
    name = "gtdbtk_on_MAGs", genomes_dir = MAGs_dir, extension = ".fasta", what = "MAGs",
    needs = [MAGs_dir + prefix + "MAGs-checkm-out.tsv"],
    out_dir = MAGs_dir + "gtdbtk-out", taxonomy_tab = MAGs_dir + prefix + "MAGs-gtdb-taxonomy.tsv",
    empty_note = MAGs_dir + "No-MAGs-recovered.txt", temp_outputs = True,
    threads = config["gtdbtk_cpus"], scratch_dir = config["gtdbtk_scratch_dir"],
    mem_mb = config["gtdbtk_mem"],
    log = logs_dir + "gtdbtk-run.log", benchmark = benchmarks_dir + "run_gtdbtk_on_MAGs-benchmarks.tsv")]

include: os.path.join(SMK, "rules", "genome-quality.smk")

tax_db_genome_runs = [] if GTDBTK_RUNS_HERE else [dict(
    name = "tax_db_on_MAGs", genomes_dir = MAGs_dir, extension = ".fasta", needs = [MAGs_dir + prefix + "MAGs-checkm-out.tsv"],
    taxonomy_tab = MAGs_taxonomy_tab, work_dir = MAGs_dir + "tax-db-working-files/",
    threads = config["threads"], search_mem = config["search_mem"], classify_seqs_mem = config["classify_seqs_mem"],
    log = logs_dir + "MAGs-tax-db.log", benchmark = benchmarks_dir + "tax_db_on_MAGs-benchmarks")]

include: os.path.join(SMK, "rules", "genome-tax.smk")


rule filter_checkm_results_and_copy_MAGs:
    input:
        bins_dir + prefix + "bins-checkm-out.tsv"
    output:
        temp(MAGs_dir + prefix + "MAGs-checkm-out.tsv")
    params:
        bins_dir = bins_dir,
        MAGs_dir = MAGs_dir,
        tmp_file = MAGs_dir + "MAGs-checkm-out.tmp",
        min_est_comp = config["min_completion"],
        max_est_redund = config["max_redundancy"],
    benchmark:
        benchmarks_dir + "filtering_checkm_results_and_copying_MAGs-benchmarks.tsv"
    shell:
        r"""
        # the MAGs of an earlier run (e.g. with other cutoffs) go first, so only this run's are picked up later
        rm -f {params.MAGs_dir}*-MAG-*.fasta

        # only running if there were bins recovered
        if [ $(find {params.bins_dir} -name "*.fasta" | wc -l | sed 's/^ *//') -gt 0 ]; then

            # keeping the header line and the bins that pass the cutoffs
            awk -F $'\t' ' NR == 1 || ( $2 >= {params.min_est_comp} && $3 <= {params.max_est_redund} ) ' {input} > {params.tmp_file}

            sed 's/-bin\./-MAG-/' {params.tmp_file} > {output}

            for MAG in $(cut -f 1 {params.tmp_file} | tail -n +2)
            do
                new_ID=$(echo $MAG | sed 's/-bin\./-MAG-/')
                cp {params.bins_dir}${{MAG}}.fasta {params.MAGs_dir}${{new_ID}}.fasta
            done

            rm {params.tmp_file}

        else

            printf "There were no MAGs recovered.\n" > {output}

        fi
        """


rule generate_bins_overview_table:
    input:
        assembly_summaries = bins_dir + prefix + "bin-assembly-summaries.tsv",
        checkm_results = bins_dir + prefix + "bins-checkm-out.tsv"
    output:
        bins_dir + prefix + "bins-overview.tsv"
    benchmark:
        benchmarks_dir + "generate_bins_overview_table-benchmarks.tsv"
    run:
        # assembly stats and checkm2 estimates in one table (or the note that no bins were recovered)
        from bit.modules.genome_overview import write_overview
        write_overview(input.assembly_summaries, output[0], checkm2_path = input.checkm_results)


rule summarize_bin_assemblies:
    conda:
        env("bit.yaml")
    input:
        trigger = binning_done
    params:
        bins_dir = bins_dir
    output:
        temp(bins_dir + prefix + "bin-assembly-summaries.tsv")
    benchmark:
        benchmarks_dir + "summarize_bin_assemblies-benchmarks.tsv"
    shell:
        r"""
        # only running if any bins were recovered
        if [ $(find {params.bins_dir} -name "*.fasta" | wc -l | sed 's/^ *//') -gt 0 ]; then

            bit summarize-assembly {params.bins_dir}*.fasta -o {output} -t

        else

            printf "There were no bins recovered.\n" > {output}

        fi
        """


# each assembly is binned with the coverage of every sample mapped to it (one sample, or all of them for a co-assembly)
rule metabat_binning:
    conda:
        env("metabat.yaml")
    input:
        assembly = lambda wildcards: asm_path(wildcards.asm),
        bams = lambda wildcards: [mapping_dir + ID + ".bam" for ID in sample_ID_list if asm_of(ID) == wildcards.asm]
    params:
        bins_dir = bins_dir,
        prefix = bins_dir + "{asm}-bin"
    threads:
        config["threads"]
    resources:
        cpus = config["threads"]
    output:
        depth_file = mapping_dir + "{asm}-metabat-assembly-depth.tsv"
    log:
        logs_dir + "{asm}-bam-summarize-and-metabat.log"
    benchmark:
        benchmarks_dir + "metabat_binning-{asm}-benchmarks.tsv"
    shell:
        r"""
        # only running if the assembly produced anything
        if [ -s {input.assembly} ]; then

            # removing any bins from a previous run of this sample
            rm -f {params.prefix}.[0-9]*.fa {params.prefix}.[0-9]*.fasta

            jgi_summarize_bam_contig_depths --outputDepth {output.depth_file} --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta {input.assembly} {input.bams} > {log} 2>&1

            # only running if there are contigs with coverage information in the coverage file we just generated
            if [ $(wc -l < {output.depth_file}) -gt 1 ]; then
                metabat2 --inFile {input.assembly} --outFile {params.prefix} --abdFile {output.depth_file} -t {threads} >> {log} 2>&1
            else
                printf "\n\nThere was no coverage info generated in {output.depth_file}, so no binning with metabat was performed.\n\n" >> {log}
            fi

            # changing extensions from .fa to .fasta to match nt fasta extension elsewhere
            for bin in {params.prefix}.[0-9]*.fa
            do
                if [ -e "${{bin}}" ]; then
                    mv "${{bin}}" "${{bin%.fa}}.fasta"
                fi
            done

        else

            touch {output}
            printf "Binning not performed because the assembly didn't produce anything.\n" > {log}

        fi
        """
