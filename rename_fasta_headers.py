from Bio import SeqIO
import sys

if sys.argv[1] == '-h':
  print("   First positional argument should be the original fasta. Second positional argument should be desired name of seqs.")
  print("   Example: python Rename_fasta_headers.py Original.fasta Seq_name")

else:

  in_fasta = open(sys.argv[1]) # original fasta

  new_header = sys.argv[2] # name you'd like for the seqs

  out_fasta = open(new_header + "_" + sys.argv[1], "w")

  n = 0
  for seq_record in SeqIO.parse(in_fasta, "fasta"):
    n = n + 1
    out_fasta.write(">" + new_header  + "_" + str(n) + "\n")
    out_fasta.write(str(seq_record.seq) + "\n")

  out_fasta.close()