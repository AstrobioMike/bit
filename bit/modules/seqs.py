from Bio import SeqIO # type: ignore
import pandas as pd # type: ignore
import gzip
from bit.modules.general import is_gzipped
from statistics import mean, median
import random

def calc_gc_per_seq(input_fasta):

    gc_stats = []
    fasta_handle = input_fasta if hasattr(input_fasta, "read") else open(input_fasta, "r")
    with fasta_handle as fasta_file:
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

    window_gc_stats = []

    fasta_handle = input_fasta if hasattr(input_fasta, "read") else open(input_fasta, "r")
    with fasta_handle as fasta_file:
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


def modify_fasta_headers(args):

    with open(args.input_fasta, "r") as in_file, open(args.output_fasta, "w") as out_file:

        if args.wanted_text:
            counter = 0
            for seq_record in SeqIO.parse(in_file, "fasta"):
                counter += 1
                out_file.write(f">{args.wanted_text}_{counter}\n")
                out_file.write(str(seq_record.seq) + "\n")

        else:

            prefix = args.prefix or ""
            suffix = args.suffix or ""

            for seq_record in SeqIO.parse(in_file, "fasta"):

                # using description so it preserves the full header
                original_header = seq_record.description

                new_header = f"{prefix}{original_header}{suffix}"
                out_file.write(f">{new_header}\n{seq_record.seq}\n")



def calc_variation_in_msa(input_alignment, type, gap_treatment="include"):

    from skbio.sequence import GrammaredSequence # type: ignore
    from skbio.util import classproperty # type: ignore
    from skbio import TabularMSA, DNA, Protein # type: ignore

    THREEDI_CHARS = set('ACDEFGHIKLMNPQRSTVWY')

    class ThreeDi(GrammaredSequence):

        @classproperty
        def definite_chars(cls):
            return THREEDI_CHARS

        @classproperty
        def degenerate_map(cls):
            return {'X': THREEDI_CHARS}

        @classproperty
        def default_gap_char(cls):
            return '-'

        @classproperty
        def gap_chars(cls):
            return set('-.')

    SEQUENCE_TYPES = {"DNA": DNA, "Protein": Protein, "3Di": ThreeDi}

    msa = TabularMSA.read(input_alignment, constructor=SEQUENCE_TYPES[type], lowercase=True)

    list_of_cleaned_seqs = []

    # converting degenerate bases to gaps
    for seq in msa:

        seq = seq.replace(seq.degenerates(), "-")
        list_of_cleaned_seqs.append(seq)

    clean_msa = TabularMSA(list_of_cleaned_seqs)

    conserved = clean_msa.conservation(gap_mode=gap_treatment)
    indexes = list(range(1,clean_msa.shape[1] + 1))

    df = pd.DataFrame({"position": indexes, "variation":1 - conserved, "conservation": conserved})

    return df


def get_gap_fracs_per_col(input_alignment):
    import numpy as np # type: ignore

    seqs = [str(r.seq).upper() for r in SeqIO.parse(input_alignment, "fasta")]
    if not seqs:
        return np.array([])

    num_seqs = len(seqs)
    num_cols = len(seqs[0])
    gap_fracs = np.array(
        [sum(1 for seq in seqs if seq[col] in ("-", ".")) / num_seqs for col in range(num_cols)]
    )
    return gap_fracs


def check_fastq_for_dup_headers(input_fastq):

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


def mutate_seq(seq, molecule_type, available_substitutions, mutation_rate, ti_tv_ratio, indel_rate):

    # converting sequence to a list for mutability
    seq_list = list(seq)
    seq_length = len(seq_list)

    # calculating mutation parameters
    total_num_mutations = round(seq_length * mutation_rate)
    num_indels = round(total_num_mutations * indel_rate)
    num_substitutions = total_num_mutations - num_indels

    # tracking counts of insertions and deletions (if any)
    num_insertions = 0
    num_deletions = 0
    num_transitions = 0
    num_transversions = 0

    substitution_indices = random.sample(range(seq_length), k=num_substitutions)

    if molecule_type == "NT":

        transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
        transversions = {
            'A': ['C', 'T'],
            'G': ['C', 'T'],
            'C': ['A', 'G'],
            'T': ['A', 'G']
        }

    for index_to_change in substitution_indices:

        original_char = seq_list[index_to_change]

        if molecule_type == "NT" and original_char in transitions:

            prob_ts = ti_tv_ratio / (ti_tv_ratio + 1)

            if random.random() < prob_ts:
                # transition
                new_char = transitions[original_char]
                num_transitions += 1
            else:
                # transversion
                new_char = random.choice([b for b in transversions[original_char] if b != original_char])
                num_transversions += 1
        else:
            # any base except original
            new_char = random.choice([char for char in available_substitutions if char != original_char])

        seq_list[index_to_change] = new_char

    # incorporating indels next, if any
    # (note this will not keep track if the same spot happens to be added and then removed as currently written)

    # limiting attempts in case of extreme edge cases
    indel_attempts = 0
    max_attempts = 1000

    while (num_insertions + num_deletions) < num_indels and indel_attempts < max_attempts:
        indel_attempts += 1

        # safeguard in the (very) off chance that all characters were deleted already
        if not seq_list:
            break

        process = random.choice(["insertion", "deletion"])
        index_to_change = random.randint(0, len(seq_list) - 1)

        if process == "insertion":
            insertion_char = random.choice(available_substitutions)
            seq_list.insert(index_to_change, insertion_char)
            num_insertions += 1
        elif process == "deletion" and len(seq_list) > 1:
            del seq_list[index_to_change]
            num_deletions += 1

    # converting list back into string
    mutated_seq = ''.join(seq_list)

    return (mutated_seq, total_num_mutations, num_substitutions,
            num_transitions, num_transversions, num_indels, num_insertions,
            num_deletions)


def dedupe_fasta_headers(input_fasta, output_fasta):

    ids = {}

    with open(input_fasta, "r") as in_fasta, open(output_fasta, "w") as out_fasta:

        for seq_record in SeqIO.parse(in_fasta, "fasta"):

            if seq_record.id not in ids:
                ids[seq_record.id] = 1
                out_fasta.write(">" + seq_record.id + "\n")
                out_fasta.write(str(seq_record.seq) + "\n")

            else:
                count = ids[seq_record.id] + 1
                ids[seq_record.id] = count
                out_fasta.write(">" + seq_record.id + "_" + str(count - 1) + "\n")
                out_fasta.write(str(seq_record.seq) + "\n")


def fasta_to_bed(input_fasta):

    bed_records = []

    fasta_handle = input_fasta if hasattr(input_fasta, "read") else open(input_fasta, "r")
    with fasta_handle as in_fasta:

        for record in SeqIO.parse(in_fasta, "fasta"):
            name = record.name
            end = len(record.seq) - 1
            bed_records.append((name, 0, end))

    return bed_records


def fasta_to_genbank(input_fasta):

    with open(input_fasta, "r") as in_fasta:

        sequences = list(SeqIO.parse(in_fasta, "fasta"))

        for seq in sequences:
            seq.annotations["molecule_type"] = "DNA"

    return sequences


def remove_wraps(input_fasta, output):

    fasta_handle = input_fasta if hasattr(input_fasta, "read") else open(input_fasta, "r")

    with fasta_handle as f:
        n = ""
        for line in f:
            if line.startswith(">"):
                output.write(n + line)
                n = ""
            else:
                output.write(line.rstrip("\n"))
                n = "\n"
        output.write(n)


def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def read_fasta(path):

    with open(path) as f:

        header = None
        seq_parts = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

        if header:
            yield header, "".join(seq_parts)
