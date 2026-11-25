# Mapping each sample's reads to an assembly (bowtie2), shared by the metagenomics and magrec workflows. Each assembly
# is indexed once, so with a co-assembly every sample's reads go to the one index.
#
# Pulled in with `include:` after these are defined:
#   sample_ID_list, paired_end, mapping_dir, logs_dir, benchmarks_dir
#   asm_of(ID)        the name of the assembly sample ID's reads are mapped to (its own, or the co-assembly)
#   asm_path(name)    that assembly's fasta
#   reads_of(ID)      the sample's read files ([R1, R2], or [reads] when single-end)
# and config holding: threads, mapping_mem
# Writes, per sample: <mapping_dir><ID>.bam (sorted and indexed), <ID>-mapping-info.txt (bowtie2's summary) and
# <ID>-mapped-reads-per-contig.tsv (samtools idxstats).

import re

asm_names = sorted({asm_of(ID) for ID in sample_ID_list})

wildcard_constraints:
    asm = "|".join(re.escape(name) for name in asm_names) or "[^/]+"


rule index_assembly:
    conda:
        env("mapping.yaml")
    input:
        lambda wildcards: asm_path(wildcards.asm)
    output:
        temp(directory(mapping_dir + "{asm}-index"))
    threads:
        config["threads"]
    resources:
        cpus = config["threads"],
        mem_mb = config["mapping_mem"]
    log:
        logs_dir + "{asm}-index.log"
    benchmark:
        benchmarks_dir + "index_assembly-{asm}-benchmarks.tsv"
    shell:
        r"""
        mkdir -p {output}
        # only if the assembly produced anything
        if [ -s {input} ]; then
            bowtie2-build --threads {threads} {input} {output}/index > {log} 2>&1
        else
            printf "Not indexed because the assembly didn't produce anything.\n" > {log}
        fi
        """


rule mapping:
    conda:
        env("mapping.yaml")
    input:
        index = lambda wildcards: mapping_dir + asm_of(wildcards.ID) + "-index",
        assembly = lambda wildcards: asm_path(asm_of(wildcards.ID)),
        reads = lambda wildcards: reads_of(wildcards.ID)
    params:
        reads_args = lambda wildcards, input: "-1 " + input.reads[0] + " -2 " + input.reads[1] if paired_end else "-U " + input.reads[0]
    threads:
        config["threads"]
    resources:
        cpus = config["threads"],
        mem_mb = config["mapping_mem"]
    output:
        bam = mapping_dir + "{ID}.bam",
        mapping_info = mapping_dir + "{ID}-mapping-info.txt",
        contig_reads = mapping_dir + "{ID}-mapped-reads-per-contig.tsv"
    log:
        logs_dir + "{ID}-mapping.log"
    benchmark:
        benchmarks_dir + "run_mapping-{ID}-benchmarks.tsv"
    shell:
        r"""
        # only running if the assembly produced anything
        if [ -s {input.assembly} ]; then

            bowtie2 --mm -q --threads {threads} -x {input.index}/index {params.reads_args} --no-unal 2> {output.mapping_info} | samtools sort -@ {threads} -o {output.bam} - 2> {log}

            # reads mapped to each contig (for the share of the reads the assembly, bins and MAGs account for)
            samtools index {output.bam} 2>> {log}
            samtools idxstats {output.bam} > {output.contig_reads} 2>> {log}

        else

            touch {output.bam} {output.mapping_info} {output.contig_reads}
            printf "Mapping not performed because the assembly didn't produce anything.\n" > {log}

        fi
        """
