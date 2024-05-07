<a href="https://github.com/AstrobioMike/bit#conda-install"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-1,300+-blue" height="22"></a>
<br>
<a href="https://github.com/AstrobioMike/bit#citation-info"><img align="right" alt="Brief paper" src="https://img.shields.io/badge/Citation%20info-blue" height="22"></a>
<br>
<a href="https://zenodo.org/badge/latestdoi/59388885"><img align="right" src="https://zenodo.org/badge/59388885.svg" alt="DOI"></a>
<br>
<a href="https://twitter.com/AstrobioMike"><img align="right" alt="Twitter Follow" src="https://img.shields.io/twitter/follow/AstrobioMike?color=blue&style=social"></a>

# Bioinformatics Tools (bit)

* [**Overview**](#overview)
  * [Programs](#programs)
  * [Workflows](#workflows)
* [**Conda install**](#conda-install)  
* [**Citation info**](#citation-info)  
* [**Shameless plug**](#shameless-plug)  

---

## Overview 
There are of course several great and widely used packages of bioinformatics helper programs out there. Some of these include the likes of [seqkit](https://github.com/shenwei356/seqkit), [seqtk](https://github.com/lh3/seqtk), [fastX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), and [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) – all of which I use regularly and have helped me do things I was trying to get done. But there are always more tasks that crop up that may not yet have a helper program or script already written that we can find.  

[*bit*](https://doi.org/10.12688/f1000research.79530.1) is a collection of one-liners, short scripts, [programs](#programs) and [workflows](#workflows) that I have been adding to over several years. Anytime I need to write something to perform a task that has more than a one-off, ad hoc use, I consider adding it here. 

*bit* runs in a Unix-like environment and is recommended to be installed with [conda](https://conda.io/docs/) as shown [below](#conda-install).  

---

### Programs
Some of the helper programs/scripts in _bit_ include:

| Program/script | Purpose | 
| ------- | ------- |
| `bit-dl-ncbi-assemblies` | downloading NCBI assemblies in different formats by just providing accession numbers |  
| `bit-get-accessions-from-GTDB` | searching the (stellar) [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) by taxonomy and getting their NCBI accessions |  
| `bit-summarize-assembly` | quickly summarizing nucleotide assemblies |  
| `bit-summarize-column` | quickly summarizing a numeric column |  
| `bit-parse-fasta-by-headers` | splitting a fasta file based on headers |  
| `bit-rename-fasta-headers` | renaming sequences in a fasta |  
| `bit-reorder-fasta` | re-ordering a fasta file |  
| `bit-extract-seqs-by-coords` | pulling out sequences from a fasta by their coordinates |  
| `bit-genbank-to-AA-seqs`, `bit-genbank-to-fasta` | pulling amino-acid or nucleotide sequences out of a GenBank file |  
| `bit-count-bases-per-seq` | counting the number of bases per sequence in a fasta file |  
| `bit-calc-variation-in-msa` | calculating [variation](http://scikit-bio.org/docs/0.5.3/generated/skbio.alignment.TabularMSA.conservation.html) in each column of a multiple-sequence alignment |  
| `bit-filter-table` | filtering a table based on wanted IDs |  
| `bit-get-lineage-from-taxids` | getting full lineage info from a list of taxon IDs (making use of the also stellar [TaxonKit](https://bioinf.shenwei.me/taxonkit/)) |  
| `bit-filter-KOFamScan-results` | filtering [KOFamScan](https://github.com/takaram/kofam_scan) results |  
| `bit-get-go-term-info` | getting information about a specific [GO](http://geneontology.org/) term |  
| `bit-summarize-go-annotations` | summarizing GO annotations |  
| `bit-kraken2-to-taxon-summaries`, `bit-combine-kraken2-taxon-summaries` | summarizing [kraken2](https://github.com/DerrickWood/kraken2) outputs in a table with counts of full taxonomic lineages, and combining multiple samples |  
| `bit-combine-bracken-and-add-lineage` | combining [bracken](https://github.com/jenniferlu717/Bracken) outputs and adding full taxonomic lineage info |  
| `bit-gen-iToL-map`, `bit-gen-iToL-colorstrip`, `bit-gen-iToL-text-dataset`, `bit-gen-iToL-binary-dataset` | generating color/mapping/data files for use with trees being viewed on the [Interactive Tree of Life](https://itol.embl.de/) site |  
| `bit-figshare-upload` | uploading a file to figshare |  

And other just convenient things that are nice to have handy, like removing soft line wraps that some fasta files have (`bit-remove-wraps`), and printing out the column names of a TSV with numbers (`bit-colnames`) to quickly see which columns we want to provide to things like `cut` or `awk` 🙂  

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.  

---

### Workflows
The [snakemake](https://snakemake.github.io/) workflows packaged with _bit_ are retrievable with `bit-get-workflow` and currently include:

| Workflow | Purpose |  
| ------- | ------- |  
| [sra-download](workflows/sra-download-wf/README.md) | downloads sra reads via prefetch and fasterq-dump, with helper program for combining run accessions if needed (see [here](workflows/sra-download-wf/README.md) for usage details) |  
| [genome-summarize](workflows/genome-summarize-wf/README.md) | generates genome assembly stats, quality estimates, and taxonomy info (see [here](workflows/genome-summarize-wf/README.md) for usage details and overview) |
| metagenomics | processes short-read metagenomics data via assembly through to merged taxonomy and KO coverage tables, and recovers and characterizes MAGs (see here for usage details and overview)

For greater detail and usage information, see the pages linked above for each workflow.

> Note that workflows are versioned independently of the _bit_ package. When you pull one with `bit-get-workflow`, the directory name will have the version, and it is also listed at the top of the Snakefile. 

---

## Conda install

> If you are new to the wonderful world of [conda](https://conda.io/docs/) and want to learn more, one place you can start learning about it is [here](https://astrobiomike.github.io/unix/conda-intro) 🙂  

Due to increasing program restrictions as *bit* has grown, it's easiest to install it in its own environment as shown below:  

```
conda create -n bit -c astrobiomike -c conda-forge -c bioconda -c defaults bit
conda activate bit
```

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.

---

## Citation info
If you happen to find *bit* useful in your work, please be sure to cite it 🙂

> Lee M. bit: a multipurpose collection of bioinformatics tools. F1000Research 2022, 11:122. [https://doi.org/10.12688/f1000research.79530.1](https://doi.org/10.12688/f1000research.79530.1)

You can get the version you are using by running `bit-version`.  

If you are using a program in *bit* that also leverages another program, please be sure to cite them too. For instance, `bit-get-lineage-from-taxids` uses [TaxonKit](https://bioinf.shenwei.me/taxonkit/), and `bit-slim-down-go-terms` uses [goatools](https://github.com/tanghaibao/goatools). For cases where a *bit* script relies on other programs like those, it will be indicated in the help menu of the *bit* program.  

---

## Shameless plug
For [phylogenomics](https://astrobiomike.github.io/genomics/phylogenomics), checkout [GToTree](https://github.com/AstrobioMike/GToTree/wiki) 🙂  

---
