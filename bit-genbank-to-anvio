#!/usr/bin/env python

from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(description="This script takes a genbank file and converts it into the format required for importing external gene calls and functions into anvi'o.")

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input_gb", help='input Genbank file (e.g. "*.gbk", "*.gb", "*.gbff")', action="store", dest="input_gb", required=True)
parser.add_argument("-s", "--annotation_source", help='Annotation source (default: "NCBI_PGAP")', action="store", dest="source", default="NCBI_PGAP")
parser.add_argument("-v", "--annotation_version", help='Annotation source version (default: "v4.6")', action="store", dest="version", default="v4.6")
parser.add_argument("-o", "--output_gene_calls_tsv", help='Output tsv file (default: "external_gene_calls.tsv")', action="store", dest="output_gene_calls_tsv", default="external_gene_calls.tsv")
parser.add_argument("-a", "--output_functions_tsv", help='Output functions file (default: "functions.tsv")', action="store", dest="output_functions_tsv", default="functions.tsv")
parser.add_argument("-f", "--output_fasta", help='Output fasta file with matching, simplified headers to be ready for `anvi-gen-contigs-db` (default: "clean.fa")', action="store", dest="output_fasta", default="clean.fa")

args = parser.parse_args()

input_gb = open(args.input_gb, "r")

output_fasta = open(args.output_fasta, "w")

source = args.source
version = args.version

output_gene_calls = open(args.output_gene_calls_tsv, "w")
output_gene_calls.write("gene_callers_id" + "\t" + "contig" + "\t" + "start" + "\t" + "stop" + "\t" + "direction" + "\t" + "partial" + "\t" + "source" + "\t" + "version" + "\n")

output_functions = open(args.output_functions_tsv, "w")
output_functions.write("gene_callers_id" + "\t" + "source" + "\t" + "accession" + "\t" + "function" + "\t" + "e_value" + "\n")

recs = [rec for rec in SeqIO.parse(input_gb, "genbank")]

num = 0 # iterator for assigning "gene_callers_id"

note_terms_to_exclude = ["frameshifted", "internal stop", "incomplete"] # dumping gene if noted as these in the "note" section of the call
location_terms_to_exclude = ["join", "<", ">"] # dumping gene if "location" section contains any of these: "join" means the gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig

for rec in recs:

  output_fasta.write(">" + rec.name  + "\n" + str(rec.seq) + "\n") # writing out new fasta with clean headers ready for anvi'o

  genes = [gene for gene in rec.features if gene.type =="CDS"] # focusing on features annotated as "CDS" by NCBI's PGAP

  for gene in genes:

    location = str(gene.location)

      # dumping gene if "location" section contains any of these terms set above: "join" means the gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig
    if any(exclusion_term in location for exclusion_term in location_terms_to_exclude):
      continue

    if "note" in gene.qualifiers:
      note = str(gene.qualifiers["note"][0])

        # dumping gene if noted as any of these in the "note" section set above
      if any(exclusion_term in note for exclusion_term in note_terms_to_exclude):
        continue

        # dumping if overlapping translation frame
    if "transl_except" in gene.qualifiers:
      continue


    num += 1 # iteration for unique "gene_callers_id" for each

      # cleaning up gene coordinates to more easily parse:
    location = location.replace("[", "")
    location = re.sub('](.*)', '', location)
    location = location.split(":")

    start = location[0] # start coordinate
    end = location[1] # end coordinate


      # setting strand to "f" or "r":
    if gene.strand == 1:
      strand="f"
    else:
      strand="r"

      # for accession, storing protein id if it has one, else the the locus tag, else "None"
    if "protein_id" in gene.qualifiers:
      acc = gene.qualifiers["protein_id"][0]
    elif "locus_tag" in gene.qualifiers:
      acc = gene.qualifiers["locus_tag"][0]
    else:
      acc = "None"

      # storing gene product annotation
    function = gene.qualifiers["product"][0]

      # if present, adding gene name to product annotation:
    if "gene" in gene.qualifiers:
      gene_name=str(gene.qualifiers["gene"][0])
      function = function + " (" + gene_name + ")"

    output_gene_calls.write(str(num) + "\t" + rec.name + "\t" + str(start) + "\t" + str(end) + "\t" + str(strand) + "\t" + "0" + "\t" + str(source) +"\t" + str(version) + "\n")
    output_functions.write(str(num) + "\t" + str(source) + "\t" + acc + "\t" + function + "\t" + "0" + "\n")

input_gb.close()
output_fasta.close()
output_gene_calls.close()
output_functions.close()
