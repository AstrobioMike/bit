from bit.modules.general import get_package_path
import bit.modules.seqs as seqs
from Bio import SeqIO
import gzip
from types import SimpleNamespace
from io import StringIO
import random

test_targets_fasta = get_package_path("tests/data/ez-screen-targets.fasta")

def test_calc_gc_per_seq():

    expected_gc_stats = [
        {"header": "yopE", "length": 660, "gc": 0.51},
        {"header": "yopK", "length": 549, "gc": 0.34},
    ]

    gc_stats = seqs.calc_gc_per_seq(input_fasta=test_targets_fasta)
    assert gc_stats == expected_gc_stats, "GC stats do not match expected values"


def test_calc_gc_sliding_window():

    expected_window_gc_stats = [
        {
            "header": "yopE",
            "length": 660,
            "gc": 0.51,
            "gc_of_windows": [0.47, 0.55, 0.57]
        },
        {
            "header": "yopK",
            "length": 549,
            "gc": 0.34,
            "gc_of_windows": [0.31, 0.36, 0.25]
        }
    ]

    window_gc_stats = seqs.calc_gc_sliding_window(input_fasta=test_targets_fasta,
                                                  window=100,
                                                  step=200)
    assert window_gc_stats == expected_window_gc_stats, "Sliding window GC stats do not match expected values"


def test_filter_fasta_by_length(tmp_path):

    in_fasta = tmp_path / "input.fasta"
    in_fasta.write_text(""">seq1
ATGC
>seq2
ATGCGT
>seq3
ATGCGTAA
""")

    out_fasta = tmp_path / "filtered.fasta"

    result = seqs.filter_fasta_by_length(
        in_fasta = str(in_fasta),
        out_fasta = str(out_fasta),
        min_length = 5,
        max_length = 8,
    )

    assert result == (3, 2, 18, 14)

    records = list(SeqIO.parse(out_fasta, "fasta"))
    assert len(records) == 2
    assert records[0].id == "seq2"
    assert str(records[0].seq) == "ATGCGT"
    assert records[1].id == "seq3"
    assert str(records[1].seq) == "ATGCGTAA"


def test_calc_variation_in_msa(tmp_path):

    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(""">seq1
ATGCATGC
>seq2
ATGCATGA
""")

    output_file = tmp_path / "variation.tsv"

    # mocking args
    args = SimpleNamespace(
        input_alignment_fasta=str(fasta_file),
        output_tsv=str(output_file),
        type="DNA",
        gap_treatment="ignore"
    )

    df = seqs.calc_variation_in_msa(args)

    assert set(df.columns) == {"position", "variation", "conservation"}
    assert len(df) == 8

    for i, row in df.iterrows():
        assert abs(row["variation"] + row["conservation"] - 1) < 1e-6

    row_8 = df[df["position"] == 8].iloc[0]
    assert abs(row_8["variation"] - 0.5) < 1e-6


def test_check_for_fastq_dup_headers(tmp_path):

    fastq_file = tmp_path / "test.fastq"
    fastq_gz_file = tmp_path / "test.fastq.gz"
    fastq_text = """@seq1
ATGC
+
IIII
@seq1
ATGC
+
IIII
"""

    fastq_file.write_text(fastq_text)

    dup_keys, seq_count = seqs.check_for_fastq_dup_headers(str(fastq_file))

    assert seq_count == 2
    assert dup_keys == ["seq1"]
    assert len(dup_keys) == 1

    # gzipped version
    with gzip.open(fastq_gz_file, "wt") as f:
        f.write(fastq_text)

    dup_keys_gz, seq_count_gz = seqs.check_for_fastq_dup_headers(str(fastq_gz_file))
    assert seq_count_gz == 2
    assert dup_keys_gz == ["seq1"]
    assert len(dup_keys_gz) == 1


def test_parse_fasta_lengths():

    fasta_content = """>seq1
ATGCATGCATGC
>seq2
ATGCATGC
>seq3
ATGC
"""

    fasta_io = StringIO(fasta_content)

    result = seqs.parse_fasta_lengths(fasta_io)

    expected_lengths = {
        "seq1": 12,
        "seq2": 8,
        "seq3": 4
    }

    assert result["lengths"] == expected_lengths
    assert result["stats"]["n_seqs"] == 3
    assert result["stats"]["min"] == 4
    assert result["stats"]["max"] == 12
    assert result["stats"]["mean"] == round((12 + 8 + 4) / 3, 2)
    assert result["stats"]["median"] == 8
    assert result["stats"]["total_length"] == 24


def test_mutate_seq():

    random.seed(9)

    seq = "AAAAAAAAAA"
    available_substitutions = ["A", "T", "C", "G"]
    mutation_rate = 0.8
    indel_rate = 0.5

    expected_mutated_seq = "ACATAGTA"
    expected_total_changes = 8
    expected_subtitutions = 4
    expected_indels = 4
    expected_insertions = 1
    expected_deletions = 3

    (mutated_seq, total_num_changes, num_subs, num_indels,
     num_insertions, num_deletions) = seqs.mutate_seq(seq, available_substitutions, mutation_rate, indel_rate)

    assert total_num_changes == expected_total_changes
    assert num_subs == expected_subtitutions
    assert num_indels == expected_indels
    assert num_insertions == expected_insertions
    assert num_deletions == expected_deletions
    assert mutated_seq == expected_mutated_seq


def test_dedupe_fasta_headers(tmp_path):
    input_fasta_content = """>seq
ACTG
>another_seq
GGTT
>seq
ACTG
"""

    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(input_fasta_content)

    output_fasta = tmp_path / "output.fasta"

    seqs.dedupe_fasta_headers(input_fasta, output_fasta)

    records = list(SeqIO.parse(output_fasta, "fasta"))

    assert len(records) == 3

    expected_ids = ["seq", "another_seq", "seq_1"]
    actual_ids = [r.id for r in records]
    assert actual_ids == expected_ids

    expected_seqs = ["ACTG", "GGTT", "ACTG"]
    actual_seqs = [str(r.seq) for r in records]
    assert actual_seqs == expected_seqs


def test_fasta_to_bed():

    bed_records = seqs.fasta_to_bed(test_targets_fasta)

    expected_bed_records = [
        ("yopE", 0, 659),
        ("yopK", 0, 548)
    ]

    assert bed_records == expected_bed_records


def test_fasta_to_genbank():

    sequences = seqs.fasta_to_genbank(test_targets_fasta)

    buffer = StringIO()
    SeqIO.write(sequences, buffer, "genbank")
    buffer.seek(0)

    records = list(SeqIO.parse(buffer, "genbank"))

    assert len(records) == 2
    assert records[0].id == "yopE"
    assert str(records[0].seq).startswith("ATGAAAATATCA")
    assert records[1].id == "yopK"
    assert str(records[1].seq).startswith("ATGTTTATTAAA")
