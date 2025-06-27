from Bio import SeqIO

def calc_gc_per_seq(input_fasta):
    """
    Parses a nucleotide multifasta and returns a list of dicts with each holding:
      - header: the record name
      - length: sequence length
      - gc: GC fraction (rounded to two decimals)
    """
    stats = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq
        length = len(seq)
        gc = round((seq.count("G") + seq.count("C")) / length, 2) if length else 0.0
        stats.append({
            "header": record.name,
            "length": length,
            "gc": gc
        })
    return stats
