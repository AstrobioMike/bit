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


def filter_fasta_by_length(in_fasta, out_fasta, min_length, max_length):

    num_initial_seqs = 0
    num_seqs_retained = 0
    num_initial_bases = 0
    num_bases_retained = 0

    with open(in_fasta, "r") as in_file, open(out_fasta, "w") as out_file:

        for seq_record in SeqIO.parse(in_file, "fasta"):

            num_initial_seqs += 1
            num_initial_bases += len(seq_record.seq)

            if (min_length is None or len(seq_record.seq) >= min_length) and len(seq_record.seq) <= max_length:

                num_seqs_retained += 1
                num_bases_retained += len(seq_record.seq)
                out_file.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")

    return (num_initial_seqs, num_seqs_retained, num_initial_bases, num_bases_retained)
