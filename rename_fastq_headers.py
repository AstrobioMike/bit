#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fastq", help="Starting fastq file", action="store", dest="input_fastq", default=True)
parser.add_argument("-w", "--desired_name", help="Name to give seqs (will append a number to each)", action="store", dest="wanted_name", default=True)

args = parser.parse_args()


in_fastq = open(args.input_fastq, "r")
new_header = args.wanted_name
out_fastq = open("Renamed_" + args.input_fastq, "w")

n = 0

for seq_record in SeqIO.parse(in_fastq, "fastq"):
	n = n + 1
	seq_record.id = new_header + "_" + str(n)
	seq_record.description = ""
	
	out_fastq.write(seq_record.format("fastq"))


in_fastq.close()
out_fastq.close()

