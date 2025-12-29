import os
from pathlib import Path

reads_dict = config['reads_json']
output_dir = config['output_dir']
base_output_dir = config['base_output_dir']
log_files_dir = config['log_files_dir']
base_log_files_dir = config['base_log_files_dir']
reads_dir = config['reads_dir']
assembler = config['assembler']
threads = config['threads']
run_fastp = config['run_fastp']
run_bbnorm = config['run_bbnorm']
memory = config['memory']
if memory.is_integer():
    memory = int(memory)
min_contig_len = config['min_contig_len']
isolate = config['isolate']
conda_env = config['conda_env']

samples = list(reads_dict.keys())


rule all:
    input:
        output_dir + "/assembly-summaries.tsv"


if run_fastp:
    post_fastp_reads_dict = {
        sample: {
            'R1': f"{output_dir}/{sample}/filtered-reads/{sample}-filtered-R1.fastq.gz",
            'R2': f"{output_dir}/{sample}/filtered-reads/{sample}-filtered-R2.fastq.gz"
        } for sample in samples
    }

    rule fastp:
        conda:
            conda_env
        input:
            R1=lambda wildcards: reads_dict[wildcards.sample]['R1'],
            R2=lambda wildcards: reads_dict[wildcards.sample]['R2']
        output:
            R1=output_dir + "/{sample}/filtered-reads/{sample}-filtered-R1.fastq.gz",
            R2=output_dir + "/{sample}/filtered-reads/{sample}-filtered-R2.fastq.gz"
        params:
            output_dir=output_dir + "/{sample}/filtered-reads",
            html=output_dir + "/{sample}/filtered-reads/{sample}-fastp.html",
            json=output_dir + "/{sample}/filtered-reads/{sample}-fastp.json"
        log:
            log_files_dir + "/{sample}-fastp.log"
        shell:
            """
            mkdir -p {output_dir}

            fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
                  --html {params.html} --json {params.json} \
                  --thread 8 --detect_adapter_for_pe > {log} 2>&1
            """
else:
    post_fastp_reads_dict = reads_dict


if run_bbnorm:
    post_bbnorm_reads_dict = {
        sample: {
            'R1': f"{output_dir}/{sample}/bbnormd-reads/{sample}-bbnormd-R1.fastq.gz",
            'R2': f"{output_dir}/{sample}/bbnormd-reads/{sample}-bbnormd-R2.fastq.gz"
        } for sample in samples
    }

    rule bbnorm:
        conda:
            conda_env
        input:
            R1=lambda wildcards: post_fastp_reads_dict[wildcards.sample]['R1'],
            R2=lambda wildcards: post_fastp_reads_dict[wildcards.sample]['R2']
        output:
            R1=output_dir + "/{sample}/bbnormd-reads/{sample}-bbnormd-R1.fastq.gz",
            R2=output_dir + "/{sample}/bbnormd-reads/{sample}-bbnormd-R2.fastq.gz"
        params:
            output_dir=output_dir + "/{sample}/bbnormd-reads",
            target="100"
        log:
            log_files_dir + "/{sample}-bbnorm.log"
        shell:
            """
            mkdir -p {params.output_dir}

            bbnorm.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} \
                     threads=8 target={params.target} > {log} 2>&1
            """
else:
    post_bbnorm_reads_dict = post_fastp_reads_dict


final_reads_dict = post_bbnorm_reads_dict


if assembler == 'megahit':
    rule megahit_assembly:
        conda:
            conda_env
        input:
            R1=lambda wildcards: final_reads_dict[wildcards.sample]['R1'],
            R2=lambda wildcards: final_reads_dict[wildcards.sample]['R2']
        output:
            output_dir + "/{sample}/{sample}-assembly.fasta"
        params:
            sample=lambda wildcards: wildcards.sample,
            threads=threads,
            memory=memory,
            min_contig_len=min_contig_len,
            assembler_out_dir=output_dir + "/{sample}/megahit-out"
        log:
            log_files_dir + "/{sample}-assembly.log"
        shell:
            """
            rm -rf {params.assembler_out_dir}
            megahit -1 {input.R1} -2 {input.R2} -m {params.memory} -t {params.threads} \
                    --min-contig-len {params.min_contig_len} \
                    -o {params.assembler_out_dir} > {log} 2>&1

            cp {params.assembler_out_dir}/final.contigs.fa {output}
            """

else:
    rule spades_assembly:
        conda:
            conda_env
        input:
            R1=lambda wildcards: final_reads_dict[wildcards.sample]['R1'],
            R2=lambda wildcards: final_reads_dict[wildcards.sample]['R2']
        output:
            output_dir + "/{sample}/spades-out/contigs.fasta"
            # output_dir + "/{sample}/{sample}-assembly.fasta"
        params:
            sample=lambda wildcards: wildcards.sample,
            threads=threads,
            mode="--isolate" if isolate else "--meta",
            assembler_out_dir=output_dir + "/{sample}/spades-out"
        log:
            log_files_dir + "/{sample}-assembly.log"
        shell:
            """
            rm -rf {params.assembler_out_dir}
            spades.py -1 {input.R1} -2 {input.R2} -t {params.threads} \
                    {params.mode} -o {params.assembler_out_dir} \
                    --phred-offset 33 > {log} 2>&1

            # cp {params.assembler_out_dir}/contigs.fasta {output}
            """

    rule filter_spades_contigs:
        input:
            rules.spades_assembly.output
        output:
            output_dir + "/{sample}/{sample}-assembly.fasta"
        params:
            min_contig_len=min_contig_len
        log:
            log_files_dir + "/{sample}-bit-filter-contigs.log"
        shell:
            """
            printf "\nFiltering spades contigs for minimum length {params.min_contig_len}\n\n" >> {log}

            bit-filter-fasta-by-length -i {input} --min-length {params.min_contig_len} -o {output} >> {log}
            """


rule summarize_assemblies:
    input:
        expand(output_dir + "/{sample}/{sample}-assembly.fasta", sample=samples)
    output:
        output_dir + "/assembly-summaries.tsv"
    shell:
        """
        bit-summarize-assembly {input} -o {output}
        """
