<a href="https://zenodo.org/badge/latestdoi/59388885"><img align="right" src="https://zenodo.org/badge/59388885.svg" alt="DOI"></a>
## Bioinformatics Tools (bit)
These are a collection of one-liners and short scripts I use frequently enough that it's been worth it for me to have them instantly available anywhere. This includes things like: 
1. quickly summarizing nucleotide assemblies (`bit-summarize-assembly`)
2. downloading NCBI assemblies in different formats by just providing accession numbers (`bit-dl-ncbi-assemblies`)
3. calculating variation in each column of a multiple sequence alignment (`bit-calc-variation-in-msa`)
4. summarizing GO annotations (`bit-summarize-go-annotations`)
5. pulling out sequences by their coordinates (`bit-extract-seqs-by-coords`)
6. splitting a fasta file based on headers (`bit-parse-fasta-by-headers`)
7. renaming sequences in a fasta (`bit-rename-fasta-headers`)
8. pulling amino acid or nucleotide sequences out of a GenBank file (`bit-genbank-to-AA-seqs` / `bit-genbank-to-fasta`)


And other just convenient things that are nice to have handy like removing those annoying soft line wraps that some fasta files have (`bit-remove-wraps`) and printing out the column names of a TSV with numbers (`bit-colnames`) to quickly see which columns need to be provided to things like `cut` or `awk` 🙂

## Conda install

In current conda environment:

```
conda install -c conda-forge -c bioconda -c defaults -c astrobiomike bit
```

Or if wanted or needed in a separate environment, that can be done with:

```
conda create -n bit -c conda-forge -c bioconda -c defaults -c astrobiomike bit
conda activate bit
```

Each command has a help menu accessible by either entering the command alone or by providing `-h` as the only argument. Once installed, you can see all available commands by entering `bit-` and pressing tab twice.

## Citation info
If you happen to find this useful in your work and are inclined to cite this toolset, thank you 🙂, and please do so with the following: 

>Lee, MD. Bioinformatics Tools (bit) 2018. doi.org/10.5281/zenodo.4287770

You can get the version you are using with `bit-version`.
