#!/usr/bin/env python

import pandas as pd
import re
import argparse
import sys

## contact: Michael D. Lee (Mike.Lee@nasa.gov)

parser = argparse.ArgumentParser(description='This script combines the outputs in our GeneLab genome standard processing.')

required = parser.add_argument_group('required arguments')

required.add_argument("-s", "--input-summary-tsv", help="Input assembly summary stats file", action="store", required=True)
required.add_argument("-c", "--input-eukcc-tsv", help="Input eukcc summaries", action="store", required=True)
required.add_argument("-t", "--input-tax-tsv", help="Input CAT taxonomies", action="store", required=True)

parser.add_argument("-o", "--output-tsv", help='Output table filename (default: "Genomes-summaries.tsv")', action="store", default="Genome-summaries.tsv")

args = parser.parse_args()

if len(sys.argv)==1:
  parser.print_help(sys.stderr)
  sys.exit(0)

# reading in summary stats
stats_df = pd.read_csv(args.input_summary_tsv, sep="\t", index_col=0)

# slimming down to those we want
wanted_summary_stats = ["Total contigs", "Total length", "Ambiguous characters", "GC content", "Maximum contig length", "Minimum contig length", "N50", "L50"]
wanted_stats_df = stats_df.loc[wanted_summary_stats, ]

# transposing
trans_df = wanted_stats_df.T

## for the life of me i can't figure out how to do this the easy way, but formatting the numbers
trans_df["Total contigs"] = trans_df["Total contigs"].map('{:,.0f}'.format)
trans_df["Total length"] = trans_df["Total length"].map('{:,.0f}'.format)
trans_df["Ambiguous characters"] = trans_df["Ambiguous characters"].map('{:,.0f}'.format)
trans_df["Maximum contig length"] = trans_df["Maximum contig length"].map('{:,.0f}'.format)
trans_df["Minimum contig length"] = trans_df["Minimum contig length"].map('{:,.0f}'.format)
trans_df["N50"] = trans_df["N50"].map('{:,.0f}'.format)
trans_df["L50"] = trans_df["L50"].map('{:,.0f}'.format)

# reading in checkm results
eukcc_df = pd.read_csv(args.input_eukcc_tsv, sep="\t", index_col=0)

# slimming down to those we want
wanted_eukcc_cols = ["Est. Comp.", "Est. Redund."]
eukcc_df = eukcc_df.loc[:, wanted_eukcc_cols]

# renaming columns
eukcc_df.columns = ["Est. Completeness (%)", "Est. Redundancy (%)"]
eukcc_df.index.names = ["Assembly"]

# merging those two
combined_df = trans_df.merge(eukcc_df, left_index=True, right_index=True)

# creating a dictionary to hold lineage info from CAT
tax_dict = {}

ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# iterating through that input file
with open(args.input_tax_tsv, "r") as tax:
    for line in tax:
        if line.strip().startswith("Assembly"):
            continue

        line = line.strip().split("\t")
        ID = line[0]
        tax_list = line[1:8]

        curr_dict = dict(zip(iter(ranks), iter(tax_list)))

        tax_dict[ID] = curr_dict

# creating dataframe from our tax dictionary
tax_df = pd.DataFrame.from_dict(tax_dict, orient="index")
tax_df.index.names = ["Assembly"]

# merging with summary stats table
final_df = combined_df.merge(tax_df, left_index=True, right_index=True)
final_df.index.names = ["Assembly"]

# changing empties to "Not Assigned"
final_df.replace({"": "Not Assigned"}, inplace=True)

# writing out
with open(args.output_tsv, "w") as out:
    out.write(final_df.to_csv(index=True, sep="\t"))
