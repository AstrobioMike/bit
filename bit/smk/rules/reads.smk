# Finding the input reads and samples, for the workflows that start from reads. Pulled in with `include:`
# after rules/common.smk; uses these settings: reads_dir, single_end, sample_ids_file.
# Defines: reads_dir, paired_end, reads_dict, sample_ID_list, sample_IDs_file, raw_reads_dir, input_suffixes,
# filtered_suffixes, and the rule that links the input files into the run directory.

sample_IDs_file = "sample-IDs.txt"

reads_dir = config["reads_dir"].rstrip("/") + "/"
paired_end = not config["single_end"]

# the read files are found the same way `bit assemble` finds them: paired files by their R1/R2 designations
# (e.g., _R1_, _R1., _1.), the sample ID being what comes before it; single-end files by their extension
from bit.modules.input_parsing import get_input_reads_dict_from_dir

reads_dict = get_input_reads_dict_from_dir(reads_dir, SE = not paired_end)
if not reads_dict:
    exit_with_message(f"No gzipped fastq files were found in '{reads_dir}' :(")

# all samples found, unless a file listing the wanted ones was given
if config["sample_ids_file"]:
    sample_ID_list = read_id_list(config["sample_ids_file"])
    not_found = [ID for ID in sample_ID_list if ID not in reads_dict]
    if not_found:
        exit_with_message(f"Read files weren't found in '{reads_dir}' for these sample IDs given in '{config['sample_ids_file']}': "
                          f"{', '.join(not_found)}. The sample IDs found there are: {', '.join(sorted(reads_dict))}")
else:
    sample_ID_list = sorted(reads_dict)

# the sample IDs being processed, for the record and for the scripts that need them
with open(sample_IDs_file, "w") as f:
    f.write("".join(ID + "\n" for ID in sample_ID_list))

# the input files are linked into the run directory under these standard names (see rule link_raw_reads),
# so the rest of the workflow doesn't depend on how they were named
raw_reads_dir = "raw-reads/"

if paired_end:
    input_suffixes = ["_R1.fastq.gz", "_R2.fastq.gz"]
    filtered_suffixes = ["_R1_filtered.fastq.gz", "_R2_filtered.fastq.gz"]
else:
    input_suffixes = [".fastq.gz"]
    filtered_suffixes = ["_filtered.fastq.gz"]

localrules: link_raw_reads

rule link_raw_reads:
    input:
        lambda wildcards: [reads_dict[wildcards.ID][which] for which in (("R1", "R2") if paired_end else ("SE",))]
    output:
        [raw_reads_dir + "{ID}" + suffix for suffix in input_suffixes]
    run:
        for source, link in zip(input, output):
            os.symlink(source, link)
