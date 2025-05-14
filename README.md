<a href="https://github.com/AstrobioMike/bit#conda-install"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-2,000+-blue" height="22"></a>
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

[*bit*](https://doi.org/10.12688/f1000research.79530.1) is a collection of one-liners, short scripts, [programs](#programs) and [workflows](#workflows) that I have been adding to over the years. Anytime I need to write something to perform a task that has more than a one-off, ad hoc use, I consider adding it here. 

*bit* runs in a Unix-like environment and is recommended to be installed with [conda](https://conda.io/docs/) as shown [below](#conda-install).  

---
**Quick start**

```
conda create -n bit -c astrobiomike -c conda-forge -c bioconda -c defaults bit
conda activate bit
```

---

### Programs
Some of the helper programs/scripts in _bit_ include:

| Program/script | Purpose | 
| ------- | ------- |
| `bit-dl-ncbi-assemblies` | download NCBI assemblies in different formats by just providing accessions |  
| `bit-get-accessions-from-GTDB` | search the (stellar) [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) by taxonomy and get their NCBI accessions |  
| `bit-summarize-assembly` | quickly summarize nucleotide assemblies |  
| `bit-ez-screen` | quickly search for nucleotide targets in nucleotide input fastas, filtered based on tunable target-coverage and percent ID thresholds, and summarized in a simple table |  
| `bit-summarize-column` | quickly summarize a numeric column |  
| `bit-mutate-seqs` | introduce point mutations (substitutions/indels) in nucleotide or amino acid fasta files |  
| `bit-count-bases-per-seq` | count the number of bases per sequence in a fasta file |  
| `bit-rename-fasta-headers` | rename sequences in a fasta |  
| `bit-parse-fasta-by-headers` | split a fasta file based on headers |  
| `bit-reorder-fasta` | re-order a fasta file |  
| `bit-extract-seqs-by-coords` | pull out sequences from a fasta by their coordinates |  
| `bit-genbank-to-cds-table` | pull out general CDS info into a tsv from a GenBank file |  
| `bit-genbank-to-AA-seqs`, `bit-genbank-to-fasta` | pull amino-acid or nucleotide sequences out of a GenBank file |  
| `bit-calc-variation-in-msa` | calculate [variation](http://scikit-bio.org/docs/0.5.3/generated/skbio.alignment.TabularMSA.conservation.html) in each column of a multiple-sequence alignment |  
| `bit-filter-table` | filter a table based on wanted IDs |  
| `bit-get-lineage-from-taxids` | get full lineage info from a list of taxon IDs (making use of the also stellar [TaxonKit](https://bioinf.shenwei.me/taxonkit/)) |  
| `bit-filter-KOFamScan-results` | filter [KOFamScan](https://github.com/takaram/kofam_scan) results |  
| `bit-get-go-term-info` | get information about a specific [GO](http://geneontology.org/) term |  
| `bit-summarize-go-annotations` | summarize GO annotations |  
| `bit-gen-kraken2-tax-plots` | generate bar plots for the most abundant taxa at each rank from a kraken2 output report file |  
| `bit-kraken2-to-taxon-summaries`, `bit-combine-kraken2-taxon-summaries` | summarize [kraken2](https://github.com/DerrickWood/kraken2) outputs in a table with counts of full taxonomic lineages, and combining multiple samples |  
| `bit-combine-bracken-and-add-lineage` | combine [bracken](https://github.com/jenniferlu717/Bracken) outputs and adding full taxonomic lineage info |  
| `bit-gen-iToL-map`, `bit-gen-iToL-colorstrip`, `bit-gen-iToL-text-dataset`, `bit-gen-iToL-binary-dataset` | generate color/mapping/data files for use with trees being viewed on the [Interactive Tree of Life](https://itol.embl.de/) site |  
| `bit-figshare-upload` | upload a file to figshare |  

And other just convenient things that are nice to have handy, like removing soft line wraps that some fasta files have (`bit-remove-wraps`), and printing out the column names of a TSV with numbers (`bit-colnames`) to quickly see which columns we want to provide to things like `cut` or `awk` 🙂  

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.  

---

### Workflows
The [snakemake](https://snakemake.github.io/) workflows packaged with _bit_ are retrievable with `bit-get-workflow` and currently include:

| Workflow | Purpose |  
| ------- | ------- |  
| [sra-download](workflows/sra-download-wf) | downloads sra reads via prefetch and fasterq-dump, with helper program for combining run accessions if needed (see [here](workflows/sra-download-wf) for usage details) |  
| [genome-summarize](workflows/genome-summarize-wf) | generates genome assembly stats, quality estimates, and taxonomy info (see [here](workflows/genome-summarize-wf) for usage details and overview) |
| [metagenomics](workflows/metagenomics-wf) | processes short-read metagenomics data via assembly through to merged taxonomy and KO coverage tables, and recovers and characterizes MAGs (see [here](workflows/metagenomics-wf) for usage details and overview)

For greater detail and usage information, see the pages linked above for each workflow.

> Note that workflows are versioned independently of the _bit_ package. When you pull one with `bit-get-workflow`, the directory name will have the version, and it is also listed at the top of the Snakefile. 

---

## Conda install

> If you are new to the wonderful world of [conda](https://conda.io/docs/) and want to learn more, one place you can start learning about it is [here](https://astrobiomike.github.io/unix/conda-intro) 🙂  

It's best to install *bit* in its own environment as shown below:  

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

If you are using a script in *bit* that also leverages another program, please be sure to cite them too. For instance, `bit-get-lineage-from-taxids` uses [TaxonKit](https://bioinf.shenwei.me/taxonkit/), and `bit-slim-down-go-terms` uses [goatools](https://github.com/tanghaibao/goatools). For cases where a *bit* script relies on other programs like those, it will be indicated in the help menu of that *bit* script.  

---

## Shameless plug
For [phylogenomics](https://astrobiomike.github.io/genomics/phylogenomics), checkout [GToTree](https://github.com/AstrobioMike/GToTree/wiki) 🙂  

---
