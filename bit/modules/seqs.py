from Bio import SeqIO
from skbio import TabularMSA, DNA, Protein
import pandas as pd
import gzip
from bit.modules.general import is_gzipped
from statistics import mean, median

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


def calc_variation_in_msa(args):

    msa = TabularMSA.read(args.input_alignment_fasta, constructor=eval(args.type), lowercase=True)

    list_of_cleaned_seqs = []

    # converting degenerate bases to gaps
    for seq in msa:

        seq = seq.replace(seq.degenerates(), "-")
        list_of_cleaned_seqs.append(seq)

    clean_msa = TabularMSA(list_of_cleaned_seqs)

    conserved = clean_msa.conservation(gap_mode=args.gap_treatment)
    indexes = list(range(1,clean_msa.shape[1] + 1))

    df = pd.DataFrame({"position": indexes, "variation":1 - conserved, "conservation": conserved})

    return df


def check_for_fastq_dup_headers(input_fastq):

    headers_dict = {}
    seq_count = 0

    if is_gzipped(input_fastq):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(input_fastq, "rt") as fastq_in:

        for seq_record in SeqIO.parse(fastq_in, "fastq"):
            seq_count += 1
            if seq_record.id in headers_dict:
                headers_dict[seq_record.id] += 1
            else:
                headers_dict[seq_record.id] = 1

    dup_keys = [k for k,v in headers_dict.items() if v > 1]

    return dup_keys, seq_count


def parse_fasta_lengths(input_fasta):
    lengths_list = []
    seq_lengths = {}

    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        length = len(seq_record.seq)
        seq_lengths[seq_record.id] = length
        lengths_list.append(length)

    stats = {
        "n_seqs": len(lengths_list),
        "min": min(lengths_list) if lengths_list else 0,
        "max": max(lengths_list) if lengths_list else 0,
        "mean": round(mean(lengths_list), 2) if lengths_list else 0,
        "median": round(median(lengths_list)) if lengths_list else 0,
        "total_length": sum(lengths_list) if lengths_list else 0
    }

    return {
        "lengths": seq_lengths,
        "stats": stats
    }
