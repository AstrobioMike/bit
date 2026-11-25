from Bio import SeqIO # type: ignore
import pandas as pd # type: ignore
import gzip
import os
import sys
from bit.modules.general import is_gzipped
from statistics import mean, median
import random


def _fasta_path(input_fasta):
    """
    A path the compiled core can read for `input_fasta`: the path itself, or for an open handle (e.g. stdin), a
    temporary copy of what it holds (removed at exit). (The core reads plain or gzipped text, any line endings.)
    """
    if isinstance(input_fasta, (str, os.PathLike)):
        return str(input_fasta)
    import atexit, tempfile
    source = input_fasta.buffer if hasattr(input_fasta, "buffer") else input_fasta
    fd, path = tempfile.mkstemp(suffix=".fasta")
    with os.fdopen(fd, "wb") as out:
        while block := source.read(1 << 20):
            out.write(block.encode() if isinstance(block, str) else block)
    atexit.register(lambda: os.path.exists(path) and os.remove(path))
    return path


def _fasta_pass(fn_name, *args):
    """ a FASTA pass of the compiled core (bit._native; see rust/src/fasta.rs); a problem with the input stops the run """
    from bit import _native
    try:
        return getattr(_native, fn_name)(*args)
    except ValueError as e:
        print(f"\n    {e}\n")
        sys.exit(1)


def calc_gc_per_seq(input_fasta):

    return [{"header": name, "length": length, "gc": round(gc / length, 2) if length else 0.0}
            for name, length, gc in _fasta_pass("fasta_records_summary", _fasta_path(input_fasta))]


def calc_gc_sliding_window(input_fasta, window, step):

    # the same round() of each window's count, done once per distinct count
    rounded = {}
    def window_gc(c):
        v = rounded.get(c)
        if v is None:
            v = rounded[c] = round(c / window, 2)
        return v

    return [{"header": name,
             "length": length,
             "gc": round(gc / length, 2) if length else 0.0,
             "gc_of_windows": [window_gc(c) for c in counts]}
            for name, length, gc, counts in _fasta_pass("fasta_gc_windows", _fasta_path(input_fasta), window, step)]


def filter_fasta_by_length(in_fasta, out_fasta, min_length, max_length):

    in_path = _fasta_path(in_fasta)
    if os.path.exists(out_fasta) and os.path.samefile(in_path, out_fasta):
        print(f"\n    The output file can't be the input file ('{out_fasta}').\n")
        sys.exit(1)
    return _fasta_pass("fasta_filter_by_length", in_path, out_fasta, min_length, max_length)


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



def _msa_seq_types():
    """The skbio sequence types calc-var-in-msa reads alignments as (skbio imported on use)."""

    from skbio.sequence import GrammaredSequence # type: ignore
    from skbio.util import classproperty # type: ignore
    from skbio import DNA, Protein # type: ignore

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

    return {"DNA": DNA, "Protein": Protein, "3Di": ThreeDi}


def _column_msa_class():
    """
    A TabularMSA whose iter_positions(ignore_metadata=True) builds each column's Sequence
    straight from the rows' strings, rather than (as skbio does) from a one-character
    Sequence object made for every cell. Each column comes out as exactly the Sequence skbio
    makes (Sequence of the column's characters, top to bottom), so skbio's conservation(),
    which walks the columns this way, computes exactly what it otherwise would, just without
    most of the per-cell overhead.
    """
    from skbio import TabularMSA, Sequence # type: ignore

    class ColumnMSA(TabularMSA):

        def iter_positions(self, reverse=False, ignore_metadata=False):
            if not ignore_metadata:
                return super().iter_positions(reverse=reverse, ignore_metadata=ignore_metadata)
            rows = [str(s) for s in self._seqs]
            columns = ("".join(chars) for chars in zip(*rows))
            if reverse:
                columns = reversed(list(columns))
            return (Sequence(c) for c in columns)

    return ColumnMSA


def _conservation_of_rows(type, rows, gap_treatment):
    """skbio's TabularMSA.conservation() of an alignment given as row strings (also run in
    worker processes, hence plain inputs)."""
    constructor = _msa_seq_types()[type]
    ColumnMSA = _column_msa_class()
    return ColumnMSA([constructor(r) for r in rows]).conservation(gap_mode=gap_treatment)


def calc_variation_in_msa(input_alignment, type, gap_treatment="include", jobs=1):

    from skbio import TabularMSA # type: ignore

    SEQUENCE_TYPES = _msa_seq_types()

    msa = TabularMSA.read(input_alignment, constructor=SEQUENCE_TYPES[type], lowercase=True)

    list_of_cleaned_seqs = []

    # converting degenerate bases to gaps
    for seq in msa:

        seq = seq.replace(seq.degenerates(), "-")
        list_of_cleaned_seqs.append(seq)

    clean_msa = TabularMSA(list_of_cleaned_seqs)

    conserved = _conservation_by_distinct_columns(clean_msa, type, gap_treatment, jobs)
    indexes = list(range(1,clean_msa.shape[1] + 1))

    df = pd.DataFrame({"position": indexes, "variation":1 - conserved, "conservation": conserved})

    return df


def _conservation_by_distinct_columns(clean_msa, type, gap_treatment, jobs=1):
    """
    skbio's clean_msa.conservation(gap_mode=gap_treatment), computed once per distinct column.

    skbio scores each column on its own, from just that column's characters (and the sequence
    type), so identical columns get identical values; this has skbio score an alignment of
    only the distinct columns (as they first appear), optionally split across `jobs` worker
    processes, and spreads the values back out to every column. Most of skbio's time goes
    to building objects for each cell as it walks the columns, so fewer columns is faster.
    """
    import numpy as np # type: ignore

    n_seqs, n_cols = clean_msa.shape
    if n_seqs == 0 or n_cols == 0:
        return clean_msa.conservation(gap_mode=gap_treatment)

    grid = np.vstack([np.frombuffer(str(s).encode("ascii"), dtype=np.uint8) for s in clean_msa])
    _, first_idx, inverse = np.unique(grid.T, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    # distinct columns in the order they first appear
    order = np.argsort(first_idx, kind="stable")
    distinct = grid[:, first_idx[order]]
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))

    n_distinct = distinct.shape[1]
    jobs = max(1, min(jobs, n_distinct))
    bounds = [round(i * n_distinct / jobs) for i in range(jobs + 1)]
    blocks = [[row[a:b].tobytes().decode("ascii") for row in distinct]
              for a, b in zip(bounds, bounds[1:])]

    if jobs == 1:
        parts = [_conservation_of_rows(type, blocks[0], gap_treatment)]
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(jobs, mp_context=mp.get_context("spawn")) as ex:
            parts = list(ex.map(_conservation_of_rows, [type] * jobs, blocks,
                                [gap_treatment] * jobs))

    distinct_values = np.concatenate(parts)
    return distinct_values[rank[inverse]]


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

    for seq_id, length, _ in _fasta_pass("fasta_records_summary", _fasta_path(input_fasta)):
        seq_lengths[seq_id] = length
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

    # the compiled core writes to the output's file descriptor (stdout's included, flushed around it)
    output.flush()
    _fasta_pass("fasta_remove_wraps", _fasta_path(input_fasta), output.fileno())


# includes IUPAC ambiguity codes, each mapping to the code for its complementary base set
# (e.g. R = A/G -> Y = T/C); S, W, and N are their own complements
_REVCOMP_TABLE = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                               "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq):
    return seq.translate(_REVCOMP_TABLE)[::-1]


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


def identify_seq_type(fasta_path):
    """
    this tries to detect if a fasta contains amino-acid or nucleotide sequences

    it reads the first sequence in the file and applies two checks:
      1. if any character is exclusive to amino acids (E, F, I, J, L, O, P, Q, Z),
         it returns "aa".
      2. elif ≥ 90 % of characters are valid IUPAC nucleotide codes,
         returns "nt"
      3. else returns "aa"
    """
    from bit.modules.general import report_message, notify_premature_exit

    AA_EXCLUSIVE = frozenset("EFIJLOPQZ")
    NT_CHARS = frozenset("ACGTUNRYSWKMBDHV")

    seq = ""
    with open(fasta_path, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            seq = str(record.seq)
            break

    clean_seq = seq.upper().replace("-", "").replace(".", "").replace("*", "")

    if not clean_seq:
        report_message(f"No sequence content found in '{fasta_path}' after removing '-', '.', and '*'.")
        notify_premature_exit()

    if any(c in AA_EXCLUSIVE for c in clean_seq):
        return "aa"

    nt_frac = sum(1 for c in clean_seq if c in NT_CHARS) / len(clean_seq)

    return "nt" if nt_frac >= 0.9 else "aa"
