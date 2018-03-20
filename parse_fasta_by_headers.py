#!/usr/bin/env python

### for parsing a fasta file by pulling out seqs with desired headers
########## example usage: python parse_fasta_by_header.py -i Sequences.fa -w Wanted_sequence_headers.txt
########## generates output file: Wanted_Sequences.fa


from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="Original fasta file", action="store", dest="input_fasta", default=True)
parser.add_argument("-w", "--wanted_sequence_headers", help="Single-column file with desired sequence headers", action="store", dest="wanted_headers", default=True)

args = parser.parse_args()

fasta_in = open(args.input_fasta, "r") # starting fasta
wanted_seqs = open(args.wanted_headers, "r") # single column text file with headers of wanted seqs

wanted_set = set(line.strip() for line in wanted_seqs)

fasta_out = open("Wanted_Sequences.fa", "w")

for seq_record in SeqIO.parse(fasta_in, "fasta"):
  if seq_record.id in wanted_set:
    fasta_out.write(">" + seq_record.id + "\n")
    fasta_out.write(str(seq_record.seq) + "\n")
    
fasta_in.close()
wanted_seqs.close()
fasta_out.close()
