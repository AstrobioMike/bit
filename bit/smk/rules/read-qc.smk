# fastqc and multiqc of the raw and the filtered reads. Pulled in with `include:` after rules/reads.smk and after
# these are defined: fastqc_out_dir, raw_fastqc_dir, filtered_fastqc_dir, filtered_reads_dir, prefix, logs_dir,
# benchmarks_dir. The filtered reads are expected as filtered_reads_dir + {ID} + each of filtered_suffixes.

multiqc_config = os.path.join(SMK, "scripts", "multiqc.config")

def fastqc_names(suffix):
    """ fastqc names its outputs after the read file, minus its last two extensions (e.g. '.fastq.gz') """
    base = suffix.rsplit(".", 2)[0]
    return [base + "_fastqc.zip", base + "_fastqc.html"]


rule raw_multiqc:
    conda:
        env("qc.yaml")
    input:
        expand(raw_fastqc_dir + "{ID}{suffix}", ID = sample_ID_list, suffix = [fastqc_names(s)[0] for s in input_suffixes])
    params:
        int_out_dir = prefix + "raw_multiqc_report",
        out_filename_prefix = prefix + "raw_multiqc",
        int_out_data_dir = prefix + "raw_multiqc_data",
        int_html_file = prefix + "raw_multiqc.html",
        int_zip = prefix + "raw_multiqc_report.zip",
        config_file = multiqc_config
    output:
        final_out_zip = fastqc_out_dir + prefix + "raw_multiqc_report.zip"
    log:
        logs_dir + "raw-multiqc.log"
    benchmark:
        benchmarks_dir + "raw_multiqc-benchmarks.tsv"
    shell:
        r"""
        multiqc -q -n {params.out_filename_prefix} --force --cl-config 'max_table_rows: 99999999' --interactive --config {params.config_file} {input} > {log} 2>&1

        # making an output report directory and moving things into it
        rm -rf {params.int_out_dir} && mkdir -p {params.int_out_dir}
        mv {params.int_html_file} {params.int_out_data_dir} {params.int_out_dir}

        # zipping and removing unzipped dir
        zip -q -r {params.int_zip} {params.int_out_dir} && rm -rf {params.int_out_dir}

        # moving to final wanted location
        mv {params.int_zip} {output.final_out_zip}
        """


# (the individual fastqc outputs are intermediate files, removed by snakemake once multiqc has summarized them)
rule raw_fastqc:
    conda:
        env("qc.yaml")
    input:
        [raw_reads_dir + "{ID}" + suffix for suffix in input_suffixes]
    output:
        temp([raw_fastqc_dir + "{ID}" + name for suffix in input_suffixes for name in fastqc_names(suffix)])
    params:
        out_dir = raw_fastqc_dir
    threads:
        len(input_suffixes)
    benchmark:
        benchmarks_dir + "raw_fastqc-{ID}-benchmarks.tsv"
    shell:
        """
        fastqc {input} --outdir {params.out_dir} -t {threads} -q
        """


use rule raw_multiqc as filtered_multiqc with:
    input:
        expand(filtered_fastqc_dir + "{ID}{suffix}", ID = sample_ID_list, suffix = [fastqc_names(s)[0] for s in filtered_suffixes])
    params:
        int_out_dir = prefix + "filtered_multiqc_report",
        out_filename_prefix = prefix + "filtered_multiqc",
        int_out_data_dir = prefix + "filtered_multiqc_data",
        int_html_file = prefix + "filtered_multiqc.html",
        int_zip = prefix + "filtered_multiqc_report.zip",
        config_file = multiqc_config
    output:
        final_out_zip = fastqc_out_dir + prefix + "filtered_multiqc_report.zip"
    log:
        logs_dir + "filtered-multiqc.log"
    benchmark:
        benchmarks_dir + "filtered_multiqc-benchmarks.tsv"


use rule raw_fastqc as filtered_fastqc with:
    input:
        [filtered_reads_dir + "{ID}" + suffix for suffix in filtered_suffixes]
    output:
        temp([filtered_fastqc_dir + "{ID}" + name for suffix in filtered_suffixes for name in fastqc_names(suffix)])
    params:
        out_dir = filtered_fastqc_dir
    benchmark:
        benchmarks_dir + "filtered_fastqc-{ID}-benchmarks.tsv"
