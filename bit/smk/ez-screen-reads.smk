from pathlib import Path
from bit.modules.ez_screen import gen_reads_summary_table, combine_reads_summary_outputs

reads_dict = config['reads_json']
base_output_prefix = config['base_output_prefix']
base_output_dir = config['base_output_dir']
mapping_output_dir = config['mapping_output_dir']
log_files_dir = config['log_files_dir']
targets_path = config['targets']
targets_base_name = Path(targets_path).stem
mapping_index_prefix = f"{mapping_output_dir}/index/{targets_base_name}"
reads_dir = config['reads_dir']
min_perc_id = config['min_perc_id']
min_perc_cov = config['min_perc_cov']

samples = list(reads_dict.keys())

rule all:
    input:
        f"{base_output_prefix}-combined-summary.tsv"
        # expand(f"{mapping_output_dir}/{{sample}}/{{sample}}-ez-screen-summary.tsv", sample=samples)
        # expand(f"{mapping_output_dir}/{{sample}}/{{sample}}.mosdepth.global.dist.txt", sample=samples)
        # expand(f"{mapping_output_dir}/{{sample}}/{{sample}}.bam", sample=samples)

rule make_bwa_index:
    input:
        targets_path
    output:
        expand(
            f"{mapping_index_prefix}.{{ext}}",
            ext=["amb", "ann", "bwt", "pac", "sa"]
        )
    params:
        mapping_index_prefix = mapping_index_prefix
    log:
        f"{log_files_dir}/{targets_base_name}-make-bwa-index.log"
    shell:
        """
        bwa index -p {params.mapping_index_prefix} {targets_path} 2> {log}
        """


rule map_reads:
    input:
        expand(
            f"{mapping_index_prefix}.{{ext}}",
            ext=["amb", "ann", "bwt", "pac", "sa"]
        )
    output:
        mapping_output_dir + "/{sample}/{sample}.bam"
    params:
        R1 = lambda wildcards: reads_dict[wildcards.sample]['R1'],
        R2 = lambda wildcards: reads_dict[wildcards.sample]['R2'],
        mapping_index_prefix = mapping_index_prefix
    shell:
        """
        ### need to check if the reads_dict has full paths or not when the reads are somewhere else
        bwa mem -t 8 -a -T 0 {params.mapping_index_prefix} {reads_dir}/{params.R1} {reads_dir}/{params.R2} 2> {log_files_dir}/{wildcards.sample}-bwa-mem.log | \
            samtools view -b | \
            samtools sort -@ 8 -o {output}

        samtools index {output}
        """


rule run_mosdepth:
    input:
        mapping_output_dir + "/{sample}/{sample}.bam"
    output:
        mapping_output_dir + "/{sample}/{sample}.mosdepth.global.dist.txt"
    shell:
        """
        mosdepth --no-per-base -x {mapping_output_dir}/{wildcards.sample}/{wildcards.sample} {mapping_output_dir}/{wildcards.sample}/{wildcards.sample}.bam
        """


rule gen_reads_summary_table:
    input:
        bam = mapping_output_dir + "/{sample}/{sample}.bam",
        global_dist_tab = mapping_output_dir + "/{sample}/{sample}.mosdepth.global.dist.txt"
    output:
        summary_tab = mapping_output_dir + "/{sample}/{sample}-ez-screen-summary.tsv"
    run:
        gen_reads_summary_table(input.bam, input.global_dist_tab, output.summary_tab)


rule combine_outputs:
    input:
        expand(
            mapping_output_dir + "/{sample}/{sample}-ez-screen-summary.tsv",
            sample=samples
        )
    output:
        final_output = f"{base_output_prefix}-combined-summary.tsv"
    run:
        samples_output_summaries_dict = {
            sample: f"{mapping_output_dir}/{sample}/{sample}-ez-screen-summary.tsv"
            for sample in samples
        }
        combine_reads_summary_outputs(samples_output_summaries_dict, output.final_output)
