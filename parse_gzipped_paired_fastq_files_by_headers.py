#!/usr/bin/env python

from Bio import SeqIO, bgzf
import sys
import argparse
from gzip import open as gzopen

parser = argparse.ArgumentParser(description="A script for pulling out reads from paired, gzipped, fastq files based on headers. \
                                Will create gzipped output files with same names but 'wanted_' prepended to them. \
                                Example usage: parse_gzipped_paired_fastq_files_by_headers.py -f forward_reads.fq.gz \
                                -r reverse_reads.fq.gz -w wanted_seqs.txt")

parser.add_argument("-f", "--in_forward_fastq", help="Original forward read fastq file", action="store", dest="f_input_fastq", default=True, required=True)
parser.add_argument("-r", "--in_reverse_fastq", help="Original reverse read fastq file", action="store", dest="r_input_fastq", default=True, required=True)
parser.add_argument("-w", "--wanted_sequence_headers", help="Single-column file with desired sequence headers", action="store", dest="wanted_headers", default=True, required=True)

args = parser.parse_args()

wanted_seqs = open(args.wanted_headers, "r")
wanted_set = set(line.strip() for line in wanted_seqs)

wanted_seqs.close()

f_fastq_in = gzopen(args.f_input_fastq, "r")

f_fastq_out = gzopen("wanted_" + args.f_input_fastq, "w")

for seq_record in SeqIO.parse(f_fastq_in, "fastq"):
    if seq_record.id in wanted_set:
        seq_record.description = ""
        f_fastq_out.write(seq_record.format("fastq"))
        
f_fastq_in.close()
f_fastq_out.close()


r_fastq_in = gzopen(args.r_input_fastq, "r")
r_fastq_out = gzopen("wanted_" + args.r_input_fastq, "w")


for seq_record in SeqIO.parse(r_fastq_in, "fastq"):
    if seq_record.id in wanted_set:
        seq_record.description = ""            
        r_fastq_out.write(seq_record.format("fastq"))

r_fastq_in.close()
r_fastq_out.close()
