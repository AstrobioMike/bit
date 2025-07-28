from bit.modules.general import get_package_path
import bit.modules.seqs as seqs
from Bio import SeqIO


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
