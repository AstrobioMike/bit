### Script that spits out a txt file of one column of the number of bases in each sequence of a fasta file ###
### Useful for instance if calculating RPKM ###

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="Original fasta file", action="store", dest="input_fasta", default=True)
parser.add_argument("-o", "--output_txt_file", help="Name of output txt file", action="store", dest="output_file", default=True)

args = parser.parse_args()

in_fasta = open(args.input_fasta, "r")
out_file = open(args.output_file, "w")

for seq_record in SeqIO.parse(in_fasta, "fasta"):
  out_file.write(str(len(seq_record.seq)) + "\n")
  
in_fasta.close()
out_file.close()


