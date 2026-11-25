# Running checkm2 and GTDB-Tk on a directory of genome fasta files (input genomes, bins, MAGs).
# Pulled in with `include:` after rules/ref-dbs.smk and after the workflow has listed what it wants run:
#
#   checkm2_runs = [dict(
#       name = "run_checkm2",              # the rule's name
#       genomes_dir = "genomes/",          # directory holding the fasta files
#       extension = ".fasta",
#       what = "genomes",                  # what they are called in messages when there are none
#       needs = [...],                     # files that must exist before it runs
#       results_tab = "checkm2-results.tsv",   # where checkm2's quality_report.tsv is copied
#       out_dir = None,                    # where checkm2's full output is kept (None: not kept)
#       temp_outputs = True,               # whether snakemake removes the outputs once they've been used
#       threads = 8, mem_mb = 50000,
#       log = "logs/checkm2.log", benchmark = "benchmarks/run_checkm2-benchmarks.tsv")]
#
#   gtdbtk_runs = [dict(
#       name = "gtdbtk_classify", genomes_dir = ..., extension = ..., what = ..., needs = [...],
#       out_dir = "gtdb-tk-output",        # GTDB-Tk's output directory
#       taxonomy_tab = "gtdb-taxonomy.tsv",    # where the genome and classification columns of its summaries are written
#       empty_note = None,                 # a file to write a note in if there are no genomes (optional)
#       temp_outputs = True, threads = 8, scratch_dir = "", mem_mb = 175000,
#       log = ..., benchmark = ...)]
#
# On macOS the GTDB-Tk rules aren't defined (see GTDBTK_RUNS_HERE in rules/ref-dbs.smk).
# If the directory holds no fasta files when a rule runs (e.g., no bins were recovered), the program isn't run
# and its outputs say so.

def _kept_or_temp(run, path):
    return temp(path) if run["temp_outputs"] else path


for _run in checkm2_runs:

    rule:
        name:
            _run["name"]
        conda:
            env("checkm2.yaml")
        input:
            needs = _run["needs"]
        params:
            genomes_dir = _run["genomes_dir"],
            extension = _run["extension"],
            what = _run["what"],
            out_dir = _run["out_dir"] or _run["results_tab"] + "-checkm2-out.tmp",
            keep_out_dir = "true" if _run["out_dir"] else "false",
            checkm2_db_file = checkm2_db_file
        threads:
            _run["threads"]
        resources:
            cpus = _run["threads"],
            mem_mb = _run["mem_mb"]
        output:
            results_tab = _kept_or_temp(_run, _run["results_tab"]),
            out_dir = [_kept_or_temp(_run, directory(_run["out_dir"]))] if _run["out_dir"] else []
        log:
            _run["log"]
        benchmark:
            _run["benchmark"]
        shell:
            r"""
            if [ $(find "{params.genomes_dir}" -maxdepth 1 -name "*{params.extension}" | wc -l | sed 's/^ *//') -gt 0 ]; then

                checkm2 predict -i "{params.genomes_dir}" -x {params.extension} -t {threads} --database_path "{params.checkm2_db_file}" --output-directory {params.out_dir} --force > {log} 2>&1

                # making copy of primary output file
                cp {params.out_dir}/quality_report.tsv {output.results_tab}

                if [ {params.keep_out_dir} != "true" ]; then
                    rm -rf {params.out_dir}
                fi

            else

                printf "There were no {params.what} recovered.\n" > {output.results_tab}
                printf "There were no {params.what} recovered, so checkm2 was not run.\n" > {log}
                if [ {params.keep_out_dir} == "true" ]; then
                    mkdir -p {params.out_dir}
                fi

            fi
            """


for _run in (gtdbtk_runs if GTDBTK_RUNS_HERE else []):

    rule:
        name:
            _run["name"]
        conda:
            env("gtdb-tk.yaml")
        input:
            needs = _run["needs"]
        params:
            genomes_dir = _run["genomes_dir"],
            extension = _run["extension"],
            what = _run["what"],
            empty_note = _run["empty_note"] or "/dev/null",
            gtdbtk_db_dir = gtdbtk_db_dir,
            scratch_arg = "--scratch_dir " + _run["scratch_dir"] if _run["scratch_dir"] else ""
        threads:
            _run["threads"]
        resources:
            cpus = _run["threads"],
            mem_mb = _run["mem_mb"]
        output:
            out_dir = _kept_or_temp(_run, directory(_run["out_dir"])),
            taxonomy_tab = _kept_or_temp(_run, _run["taxonomy_tab"])
        log:
            _run["log"]
        benchmark:
            _run["benchmark"]
        shell:
            r"""
            # gtdb-tk finds its reference data through this variable
            export GTDBTK_DATA_PATH="{params.gtdbtk_db_dir}"

            if [ $(find "{params.genomes_dir}" -maxdepth 1 -name "*{params.extension}" | wc -l | sed 's/^ *//') -gt 0 ]; then

                gtdbtk classify_wf {params.scratch_arg} --genome_dir "{params.genomes_dir}" -x {params.extension} --out_dir {output.out_dir} --cpus {threads} --pplacer_cpus 1 > {log} 2>&1

            else

                mkdir -p {output.out_dir}
                printf "There were no {params.what} recovered.\n" > {params.empty_note}
                printf "\n\nThere were no {params.what} recovered, so GTDB-tk was not run.\n\n" > {log}

            fi

            # the genome and classification columns of the bacterial and archaeal summaries, in one table
            printf "user_genome\tclassification\n" > {output.taxonomy_tab}
            for summary in {output.out_dir}/gtdbtk.*.summary.tsv; do
                if [ -f ${{summary}} ]; then
                    tail -n +2 ${{summary}} | cut -f 1,2 >> {output.taxonomy_tab}
                fi
            done
            """
