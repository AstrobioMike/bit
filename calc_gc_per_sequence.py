### this is for nucleotide multifastas and will return a tab-delimited file with 
### 3 columns: header, sequence length, and gc 

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", help="fasta file", action="store", dest="input_fasta", default=True)
parser.add_argument("-o", "--output_txt_file", help="Name of output txt file", action="store", dest="output_file", default=True)

args = parser.parse_args()

in_fasta = open(args.input_fasta, "r")
out_file = open(args.output_file, "w")

out_file.write("header" + "\t" + "length" + "\t" + "gc" + "\n")

for cur_record in SeqIO.parse(in_fasta, "fasta"):
  gene_name = cur_record.name
  A_count = cur_record.seq.count('A')
  C_count = cur_record.seq.count('C')
  G_count = cur_record.seq.count('G')
  T_count = cur_record.seq.count('T')
  length = len(cur_record.seq)
  gc_percentage = float(G_count + C_count) / length
  gc_percentage = round(gc_percentage,2)
  out_file.write(str(gene_name)+"\t"+str(length)+"\t"+str(gc_percentage)+"\n")

in_fasta.close()
out_file.close()
