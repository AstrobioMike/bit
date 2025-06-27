from Bio import SeqIO

def calc_gc_per_seq(input_fasta):
    """
    Parses a nucleotide multifasta and returns a list of dicts with each holding:
      - header: the record name
      - length: sequence length
      - gc: GC fraction (rounded to two decimals)
    """
    gc_stats = []
    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq
            length = len(seq)
            gc = round((seq.count("G") + seq.count("C")) / length, 2) if length else 0.0
            gc_stats.append({
                "header": record.name,
                "length": length,
                "gc": gc
            })
    return gc_stats


def calc_gc_sliding_window(input_fasta, window, step):
    """
    Parses a nucleotide multifasta and returns a list of dicts with each holding:
      - header: the record name
      - length: sequence length
      - gc: GC fraction (rounded to two decimals)
      - gc_windows: list of GC fractions for each sliding window
    """

    window_gc_stats = []

    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq
            length = len(seq)
            gc = round((seq.count("G") + seq.count("C")) / length, 2) if length else 0.0

            gc_of_windows = []
            for i in range(0, length - window + 1, step):
                window_seq = seq[i:i + window]
                window_gc = round((window_seq.count("G") + window_seq.count("C")) / window, 2) if window else 0.0
                gc_of_windows.append(window_gc)

            window_gc_stats.append({
                "header": record.name,
                "length": length,
                "gc": gc,
                "gc_of_windows": gc_of_windows
            })

    return window_gc_stats
