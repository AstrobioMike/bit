#!/usr/bin/env python3

from pybedtools import BedTool
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="Starting fasta file", action="store", dest="input_fasta", default=True)
parser.add_argument("-b", "--bed_file", help="Bed file of desired contigs and coordinates (3 columns – contig, start, end – no header, 0-based counting", dest="bed_file", default=True)
parser.add_argument("-o", "--output_fasta", help="Name of output fasta file", action="store", dest="output_fasta", default=True)

args = parser.parse_args()

coordinates_file = BedTool(args.bed_file)
fasta = BedTool(args.input_fasta)
seq = coordinates_file.sequence(fi=fasta)

out_fasta = open(args.output_fasta, "w")
out_fasta.write(open(seq.seqfn).read())

out_fasta.close()
