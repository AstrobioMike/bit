#!/usr/bin/env bash

### Helper script for `bit-cov-summarize-go-annots-with-domains`; for version info, run or see `bit-version` ###
### generates a table grouping all taxids from Euks, Bacteria, Archaea, and viruses

# getting domain info for all taxids in ncbi files and storing in taxonkit data dir
cut -f 1 ${TAXONKIT_DB}/nodes.dmp | taxonkit lineage | taxonkit reformat -r NA | cut -f 1,3 | tr ";" "\t" | cut -f 1,2 | grep -v "NA" > ${TAXONKIT_DB}/taxids-and-domains.tsv
