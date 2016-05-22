from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="Original fasta", action="store", dest="input_fasta", default=True)
parser.add_argument("-s", "--wanted_seqs", help="Text file where each line is header of desired seqs", action="store", dest="wanted_seqs", default=True)
parser.add_argument("-o", "--output_fasta", help="Name of fasta of selected seqs.", action="store", dest="output_fasta", default=True)

args = parser.parse_args()


in_fasta = open(args.input_fasta, "r") # original fasta
wanted_seqs = []  # creating empty list
wanted_seqs_file = open(args.wanted_seqs, 'r') # wanted seqs flat file
  
for line in wanted_seqs_file:
  wanted_seqs.append(line.strip())  # adding wanted seqs to list
  
out_fasta = open(args.output_fasta, "w")

for seq_record in SeqIO.parse(in_fasta, "fasta"):     # iterating through original fasta
  if any (seq in seq_record.id for seq in wanted_seqs):   # pulling out seqs with headers that match list
    out_fasta.write(">" + seq_record.id + "\n")
    out_fasta.write(str(seq_record.seq) + "\n")

wanted_seqs_file.close()
in_fasta.close()
out_fasta.close()