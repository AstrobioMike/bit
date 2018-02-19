#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="Starting fasta file", action="store", dest="input_fasta", default=True)
parser.add_argument("-w", "--desired_name", help="Name to give seqs (will append a number to each)", action="store", dest="wanted_name", default=True)

args = parser.parse_args()


in_fasta = open(args.input_fasta, "r")
new_header = args.wanted_name
out_fasta = open("Renamed_" + args.input_fasta, "w")

n = 0

for seq_record in SeqIO.parse(in_fasta, "fasta"):
	n = n + 1
	out_fasta.write(">" + new_header  + "_" + str(n) + "\n")
	out_fasta.write(str(seq_record.seq) + "\n")

in_fasta.close()
out_fasta.close()

