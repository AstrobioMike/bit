# [bit](https://github.com/AstrobioMike/bit) sra-download workflow
This is a [snakemake](https://snakemake.github.io/) workflow for downloading reads from [NCBI's SRA](https://www.ncbi.nlm.nih.gov/sra) in fastq format. For all workflows available with _bit_, see [here](https://github.com/AstrobioMike/bit?tab=readme-ov-file#workflows).

---

* [**Overview**](#overview)
* [**Usage**](#usage)
  * [Retrieving the workflow](#retrieving-the-workflow)
  * [Creating the input file and modifying the config.yaml](#creating-the-input-file-and-modifying-the-configyaml)
  * [Running the workflow](#running-the-workflow)
* [**Version info**](#version-info)

---

## Overview

This workflow will download reads from SRA based on input run accessions (i.e., the accessions starting with ERR..., SRR, or DRR) using prefetch and fasterq-dump.

---

## Usage
_bit_ should be installed via conda as described [here](https://github.com/AstrobioMike/bit?tab=readme-ov-file#conda-install).

### Retrieving the worklfow

```bash
bit-get-workflow sra-download
```

### Creating the input file and modifying the config.yaml
Before running it, you first need to make a file holding the target run accessions, one per line in a single-column.

The path to that file needs to be set for the "target_sra_accessions_file" variable in the config.yaml.

### Running the workflow
After the target run accessions file has been created and set in the config.yaml, here's an example of how it could be run (note that it should still be run inside the _bit_ conda environment):
 
```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 4 -p
```

- `--use-conda` – this specifies to use the conda environments included in the workflow
- `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
- `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
- `-p` – specifies to print out each command being run to the screen

See `snakemake -h` for more options and details.

---

## Version info
Note that the workflows are version independently of the _bit_ package. When you pull one with `bit-get-workflow`, the directory name will have the version, and it is also listed at the top of the Snakefile.

All versions of programs used can be found in their corresponding conda yaml file in the envs/ directory. 
