# [bit](https://github.com/AstrobioMike/bit) amplicon workflow
This is a [snakemake](https://snakemake.github.io/) workflow for processing short-read amplicon/marker-gene data (16S/18S/ITS) through to ASVs, a count table, and a taxonomy table. For all workflows available with _bit_, see [here](https://github.com/AstrobioMike/bit?tab=readme-ov-file#workflows).

---

* [**Overview**](#overview)
* [**Usage**](#usage)
  * [Retrieving the workflow](#retrieving-the-workflow)
  * [Creating the input file and modifying the config.yaml](#creating-the-input-file-and-modifying-the-configyaml)
  * [Running the workflow](#running-the-workflow)
* [**Primary outputs**](#primary-outputs)
* [**Version info**](#version-info)

---

## Overview

This workflow processes short-read amplicon data (16S/18S/ITS) through to ASVs, a count table, and a taxonomy table via the following programs and reference databases:

  - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)/[multiqc](https://multiqc.info/) for read-quality assessment and summarization
  - [cutadapt](https://cutadapt.readthedocs.io/en/stable/) for primer trimming
  - [dada2](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) for ASV inference and quality filtering
  - [decipher](https://decipher.codes/) for taxonomic classification of ASVs with the following [databases](https://decipher.codes/Downloads.html)
    - SILVA SSU r128.2 for 16S
    - PR2 18S v4.13 for 18S
    - UNITE 2025 for ITS

---

## Usage
_bit_ should be installed via conda as described [here](https://github.com/AstrobioMike/bit?tab=readme-ov-file#conda-install).

### Retrieving the worklfow

```bash
bit get-workflow amplicon
```

### Creating the input file and modifying the config.yaml
Before running it, you first need to make a file holding unique portions of the filenames for all input samples, one per line in a single-column. And some variables need to be set in the config.yaml.

The primary things that need to be set in the config.yaml are designated in the first block of the config.yaml, these include things like: the unique sample ID file mentioned just above; paired-end or single-end; where the input reads are located; their expected suffix; and details about the primers.

There are many other options/settings that can be changed if wanted, all described in the config.yaml.

### Running the workflow
After variables are set in the config.yaml, here's an example of how it could be run (note that it should still be run inside the _bit_ conda environment):
 
```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 4 -p
```

- `--use-conda` – this specifies to use the conda environments included in the workflow
- `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
- `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
- `-p` – specifies to print out each command being run to the screen

See `snakemake -h` for more options and details.

---

## Primary outputs
A primary output directory is produced called "workflow-outputs", with the following sub-directories and contents:

- **`fastqc-outputs/`**
  - multiqc summaries of fastqc reports for input and filtered reads
- **`filtered-reads/`**
  - quality-filtered reads
- **`trimmed-reads/`**
  - primer-trimmed reads (if primers were trimmed by the workflow)
- **`final-outputs/`**
  - read-count-tracking.tsv
  - ASVs.fasta
  - counts.tsv
  - taxonomy.tsv
  - taxonomy-and-counts.tsv
  - taxonomy-and-counts.biom

Please feel free to post an issue or email with any confusion about any of the outputs produced!

> If you'd like a small test dataset to run, you can run `bit data get test-data amplicon` to grab 2 small amplicon samples. 

### Other outputs
A `benchmarks/` directory will hold time and resource utlization info (as described [here](https://stackoverflow.com/a/66872577) for most steps performed.

---

## Version info
Note that the workflows are version independently of the _bit_ package. When you pull one with `bit get-workflow`, the directory name will have the version, and it is also listed at the top of the Snakefile.

All versions of programs used can be found in their corresponding conda yaml file in the envs/ directory. 
