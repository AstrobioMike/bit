########## for parsing a fasta file by pulling out seqs with desired headers #########
########## example: python parse_fasta_by_header.py Sequences.fasta Wanted_sequence_headers.txt
########## generates output file: Wanted_Sequences.fasta


from Bio import SeqIO
import sys

fasta_in = open(sys.argv[1]) # starting fasta
wanted_seqs = open(sys.argv[2]) # single column text file with headers of wanted seqs

wanted_set = set(line.strip() for line in wanted_seqs)

fasta_out = open("Wanted_" + sys.argv[1], "w")

for seq_record in SeqIO.parse(fasta_in, "fasta"):
  if seq_record.id in wanted_set:
    fasta_out.write(">" + seq_record.id + "\n")
    fasta_out.write(str(seq_record.seq) + "\n")
    
fasta_in.close()
wanted_seqs.close()
fasta_out.close()