<a href="https://zenodo.org/badge/latestdoi/59388885"><img align="right" src="https://zenodo.org/badge/59388885.svg" alt="DOI"></a>
<br>
<a href="https://github.com/AstrobioMike/bit/edit/master/README.md#conda-install"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-750+-blue" height="23"></a>
# Bioinformatics Tools (bit)

* [**Overview**](#overview)  
* [**Conda install**](#conda-install)  
* [**Citation info**](#citation-info)  
* [**Shameless plug**](#shameless-plug)  

---

## Overview 
There are of course several great and widely used packages of bioinformatics helper programs out there. Some of these include the likes of [seqtk](https://github.com/lh3/seqtk), [fastX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), and [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) – all of which I use regularly and have helped me do the things I was trying to get done. But there are always more tasks that crop up that may not yet have a helper program or script already written that we can find.  

[*bit*](https://doi.org/10.12688/f1000research.79530.1) is a collection of one-liners, short scripts, and programs that run in a Unix-like command-line environment that I have been adding to over several years. Anytime I need to write something to perform a task that has more than a one-off, ad hoc use, I consider adding it here. This includes things like:

| Purpose | Script(s) | 
| ------- | ------- |
| quickly summarizing nucleotide assemblies | `bit-summarize-assembly` |  
| splitting a fasta file based on headers | `bit-parse-fasta-by-headers` |  
| renaming sequences in a fasta | `bit-rename-fasta-headers` |  
| re-ordering a fasta file | `bit-reorder-fasta` |  
| pulling out sequences from a fasta by their coordinates | `bit-extract-seqs-by-coords` |  
| pulling amino-acid or nucleotide sequences out of a GenBank file | `bit-genbank-to-AA-seqs`, `bit-genbank-to-fasta` |  
| counting the number of bases per sequence in a fasta file | `bit-count-bases-per-seq` |  
| calculating [variation](http://scikit-bio.org/docs/0.5.3/generated/skbio.alignment.TabularMSA.conservation.html) in each column of a multiple-sequence alignment | `bit-calc-variation-in-msa` |  
| filtering a table based on wanted IDs | `bit-filter-table` |  
| downloading NCBI assemblies in different formats by just providing accession numbers | `bit-dl-ncbi-assemblies` |  
| searching the (stellar) [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) by taxonomy and getting their NCBI accessions | `bit-get-accessions-from-GTDB` |  
| getting full lineage info from a list of taxon IDs (making use of the also stellar [TaxonKit](https://bioinf.shenwei.me/taxonkit/)) | `bit-get-lineage-from-taxids` |  
| filtering [KOFamScan](https://github.com/takaram/kofam_scan) results | `bit-filter-KOFamScan-results` |  
| getting information about a specific [GO](http://geneontology.org/) term | `bit-get-go-term-info` |  
| summarizing GO annotations | `bit-summarize-go-annotations` |  
| summarizing [kraken2](https://github.com/DerrickWood/kraken2) outputs in a table with counts of full taxonomic lineages, and combining multiple samples | `bit-kraken2-to-taxon-summaries`, `bit-combine-kraken2-taxon-summaries` |  
| combining [bracken](https://github.com/jenniferlu717/Bracken) outputs and adding full taxonomic lineage info | `bit-combine-bracken-and-add-lineage` |  
| generating color/mapping/data files for use with trees being viewed on the [Interactive Tree of Life](https://itol.embl.de/) site | `bit-gen-iToL-map`, `bit-gen-iToL-colorstrip`, `bit-gen-iToL-text-dataset`, `bit-gen-iToL-binary-dataset` |  
| upload a file to figshare | `bit-figshare-upload` |  


And other just convenient things that are nice to have handy, like removing soft line wraps that some fasta files have (`bit-remove-wraps`), and printing out the column names of a TSV with numbers (`bit-colnames`) to quickly see which columns we want to provide to things like `cut` or `awk` 🙂  

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.  

*bit* runs in a Unix-like environment and is recommended to be installed with [conda](https://conda.io/docs/) as shown below.  

---

## Conda install

> If you are new to the wonderful world of [conda](https://conda.io/docs/) and want to learn more, one place you can start learning about it is [here](https://astrobiomike.github.io/unix/conda-intro) 🙂  

Due to increasing program restrictions as *bit* has grown, it's easiest to install it in its own environment as shown below:  

```
conda create -n bit -c conda-forge -c bioconda -c defaults -c astrobiomike bit
conda activate bit
```

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.

---

## Citation info
If you happen to find *bit* useful in your work, please be sure to cite it 🙂

> Lee M. bit: a multipurpose collection of bioinformatics tools. F1000Research 2022, 11:122. [https://doi.org/10.12688/f1000research.79530.1](https://doi.org/10.12688/f1000research.79530.1)

You can get the version you are using by running `bit-version`.  

If you are using a program in *bit* that also leverages another program, please be sure to cite them too. For instance, `bit-get-lineage-from-taxids` uses [TaxonKit](https://bioinf.shenwei.me/taxonkit/), and `bit-slim-down-go-terms` used [goatools](https://github.com/tanghaibao/goatools). For cases where a *bit* script relies on other programs like those, it will be indicated in the help menu of the *bit* program.  

---

## Shameless plug
For [phylogenomics](https://astrobiomike.github.io/genomics/phylogenomics), checkout [GToTree](https://github.com/AstrobioMike/GToTree/wiki) 🙂  

---
