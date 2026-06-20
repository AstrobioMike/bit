<a href="https://github.com/AstrobioMike/bit#conda-install"><img align="right" alt="Conda installs" src="https://img.shields.io/badge/Conda%20installs-2,500+-blue" height="22"></a>
<br>
<a href="https://github.com/AstrobioMike/bit#citation-info"><img align="right" alt="Brief paper" src="https://img.shields.io/badge/Citation%20info-blue" height="22"></a>
<br>
<a href="https://doi.org/10.5281/zenodo.3383647"><img align="right" alt="DOI" src="https://img.shields.io/badge/DOI-blue" height="22"></a>
<br>
<a href="https://twitter.com/AstrobioMike"><img align="right" alt="Twitter Follow" src="https://img.shields.io/twitter/follow/AstrobioMike?color=blue&style=social"></a>

# bit: a multipurpose collection of bioinformatics tools

* [**Overview**](#overview)
  * [**Programs**](#programs)
    * [NCBI/GTDB-related](#ncbigtdb-related)
    * [Coverage/mapping-related](#coveragemapping-related)
    * [Sequence manipulation / read generation](#sequence-manipulation--read-generation)
    * [Sequence searching/comparing](#sequence-searchingcomparing)
    * [Fasta utilities](#fasta-utilities)
    * [Assembly-related](#assembly-related)
    * [GenBank-format utilities](#genbank-format-utilities)
    * [Taxonomy and lineage helpers](#taxonomy-and-lineage-helpers)
    * [Table utilities](#table-utilities)
    * [Functional-annotation helpers](#functional-annotation-helpers)
    * [iToL helpers](#itol-helpers)
    * [bit-data management](#bit-data-management)
  * [**Workflows**](#workflows)
    * [SRA-download](workflows/sra-download-wf)
    * [Genome-summarize](workflows/genome-summarize-wf)
    * [Metagenomics](workflows/metagenomics-wf)
    * [Amplicon](workflows/amplicon-wf)
* [**Conda install**](#conda-install)
* [**Citation info**](#citation-info)
* [**Shameless plug**](#shameless-plug)

---

## Overview 
There are of course several great and widely used packages of bioinformatics helper programs out there. Some of these include the likes of [seqkit](https://github.com/shenwei356/seqkit), [seqtk](https://github.com/lh3/seqtk), [fastX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), and [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) – all of which I use regularly and have helped me do things I was trying to get done. But there are always more tasks that crop up that may not yet have a helper program or script already written that we can find.  

[*bit*](https://doi.org/10.12688/f1000research.79530.1) is a collection of [programs](#programs) and [workflows](#workflows) that I have been adding to over the years. Anytime I need to write something to perform a task that has more than a one-off, ad hoc use, I consider adding it here. 

*bit* runs in a Unix-like environment and is recommended to be installed with [conda](https://conda.io/docs/) as shown [below](#conda-install).  

---
**Quick start**

```
conda create -n bit -c astrobiomike -c conda-forge -c bioconda -c defaults bit
conda activate bit
```

---

### Programs

You can see an overview of available programs like below by running `bit` by itself at the command line once it's installed. 

#### NCBI/GTDB-related

| Program | Purpose |
| ------- | ------- |
| `bit get-accs-from-gtdb` | search the [GTDB](https://gtdb.ecogenomic.org/) by taxonomy and retrieve NCBI accessions |
| `bit get-accs-from-ncbi` | search the NCBI by taxonomy or taxid and retrieve NCBI accessions |
| `bit dl-ncbi-assemblies` | download NCBI assemblies in different formats given input accessions |

---

#### Simulation / sequence manipulation

| Program | Purpose |
| ------- | ------- |
| `bit gen-metagenome` | generate reads from fasta files |
| `bit gen-reads` | generate reads from fasta files |
| `bit mutate-seqs` | introduce point mutations (substitutions/indels) into nucleotide or amino-acid fasta files |
| `bit add-insertion` | add insertions into nucleotide or amino-acid fasta sequences |

---

#### Coverage/mapping-related 

| Program | Purpose |
| ------- | ------- |
| `bit cov-analyzer` | analyze coverage patterns from a bam + reference fasta to identify regions of relatively higher or lower coverage |
| `bit cov-stats` | get detection, coverage, and mean percent ID for single or multiple references given fasta(s) and a bam |
| `bit mapped-read-stats` | get percent ID and other information for mapped reads in a bam |

---

#### Sequence searching/comparing

| Program | Purpose |
| ------- | ------- |
| `bit aa-diff` | compare a query sequence to an amino-acid reference and report differences |

##### Program: `bit ez-screen`
| Subcommand | Purpose |
| ------- | ------- |
| `assembly` | runs blast-based screening of targets in assemblies |
| `reads` | runs mapping-based screening of reads against targets |

---

#### Fasta utilities

##### Program: `bit fasta`
| Subcommand | Purpose |
| ---------- | ------- |
| `calc-gc` | calculate GC content per sequence or for the full file |
| `calc-var-in-msa` | calculate [variation](https://scikit.bio/docs/dev/generated/skbio.alignment.TabularMSA.conservation.html) in each column of a multiple-sequence alignment |
| `count` | count and summarize bases or sequences |
| `extract-by-coords` | extract sequences by genomic coordinates |
| `extract-by-headers` | extract sequences by header names |
| `extract-by-primers` | extract sequences based on primer sequences |
| `filter-by-length` | filter sequences by minimum/maximum length |
| `modify-headers` | rename or reformat sequence headers |
| `remove-wraps` | remove soft line wraps |
| `to-bed` | convert fasta to BED format |
| `to-genbank` | convert fasta to GenBank format |

---

#### Assembly-related

| Program | Purpose |
| ------- | ------- |
| `bit assemble` | simple wrapper for assembly with optional quality trimming and normalization |
| `bit summarize-assembly` | quickly summarize nucleotide assemblies |

---


#### GenBank-format utilities

##### Program: `bit genbank`
| Subcommand | Purpose |
| ---------- | ------- |
| `to-AA-seqs` | extract amino acid sequences |
| `to-cds-tsv` | extract CDS info to a TSV |
| `to-cds-seqs` | extract CDS nucleotide sequences |
| `to-fasta` | extract nucleotide sequences |

---

#### Taxonomy and lineage helpers

##### Program: `bit kraken2`
| Subcommand | Purpose |
| ---------- | ------- |
| `tax-plots` | generate standard taxonomy barplots from kraken2/bracken outputs |
| `tax-summary` | generate summary tables from kraken2/bracken outputs |

##### Program: `bit lineage`
| Subcommand | Purpose |
| ---------- | ------- |
| `from-taxids` | get full lineage info from a list of NCBI taxon IDs |
| `to-tsv` | reformat lineage info to a TSV |


---

#### Table utilities

##### Program: `bit table`
| Subcommand | Purpose |
| ---------- | ------- |
| `colnames` | print column names with numbers (handy for `cut`/`awk`) |
| `filter` | filter a table based on wanted strings |
| `normalize` | normalize to CPM or with the DESeq2 median-ratio method |
| `summarize-column` | summarize a numeric column |

---

#### Functional-annotation helpers

| Program | Purpose |
| ------- | ------- |
| `bit filter-ko-results` | filter [KOFamScan](https://github.com/takaram/kofam_scan) results |


##### Program: `bit go`
Subcommand | Purpose |
| ---------- | ------- |
| `get-term-info` | look up GO term info |
| `summarize-annotations` | summarize GO annotations |
| `combine-summaries` | combine GO summary outputs |
| `slim-terms` | slim GO terms to a specified ontology |

---

#### iTOL-helpers

##### Program: `bit itol`
| Subcommand | Purpose |
| ---------- | ------- |
| `binary-dataset` | generate a binary dataset annotation file |
| `colorstrip` | generate a color strip annotation file |
| `map` | generate a mapping/connection file |
| `text-dataset` | generate a text label dataset file |

---

#### bit-data management

##### Program: `bit data`
| Subcommand | Purpose |
| ---------- | ------- |
| `get` | download or update bit-utilized reference databases, or grab test data |
| `locations` | check or set data-location environment variables |


---

Each subcommand has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. You can see an overview of all available subcommands by running `bit` by itself. 

---

### Workflows
The [snakemake](https://snakemake.github.io/) workflows packaged with _bit_ are retrievable with `bit get-workflow` and currently include:

| Workflow | Purpose |  
| ------- | ------- |  
| [sra-download](workflows/sra-download-wf) | downloads sra reads via prefetch and fasterq-dump, with helper program for combining run accessions if needed (see [here](workflows/sra-download-wf) for usage details) |  
| [genome-summarize](workflows/genome-summarize-wf) | generates genome assembly stats, quality estimates, and taxonomy info (see [here](workflows/genome-summarize-wf) for usage details and overview) |
| [metagenomics](workflows/metagenomics-wf) | processes short-read metagenomics data via assembly through to merged taxonomy and KO coverage tables, and recovers and characterizes MAGs (see [here](workflows/metagenomics-wf) for usage details and overview)
| [amplicon](workflows/amplicon-wf) | processes short-read amplicon data through to ASVs, a count table, and a taxonomy table largely with dada2 (see [here](workflows/amplicon-wf) for usage details and overview)

For greater detail and usage information, see the pages linked above for each workflow.

> Note that workflows are versioned independently of the _bit_ package. When you pull one with `bit get-workflow`, the directory name will have the version, and it is also listed at the top of the Snakefile. 

---

## Conda install

> If you are new to the wonderful world of [conda](https://conda.io/docs/) and want to learn more, one place you can start learning about it is [here](https://astrobiomike.github.io/unix/conda-intro) 🙂  

It's best to install *bit* in its own environment as shown below:  

```
conda create -n bit -c astrobiomike -c conda-forge -c bioconda -c defaults bit
conda activate bit
```

Once installed, each subcommand has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. You can see an overview of all available subcommands by running `bit` by itself.

---

## Citation info
If you happen to find *bit* useful in your work, please be sure to cite it 🙂

> Lee M. bit: a multipurpose collection of bioinformatics tools. F1000Research 2022, 11:122. [https://doi.org/10.12688/f1000research.79530.1](https://doi.org/10.12688/f1000research.79530.1)

You can get the version you are using by adding the `-v` or `--version` flag to any of the programs.

If you are using a script in *bit* that also leverages another program, please be sure to cite them too. For instance, `bit-get-lineage-from-taxids` uses [TaxonKit](https://bioinf.shenwei.me/taxonkit/), and the `bit-go` subcommands use [goatools](https://github.com/tanghaibao/goatools). For cases where *bit* relies on other programs like those, it will be indicated in the help menu of that *bit* program.  

---

## Shameless plug
For [phylogenomics](https://astrobiomike.github.io/genomics/phylogenomics), checkout [GToTree](https://github.com/AstrobioMike/GToTree/wiki) 🙂  

---
