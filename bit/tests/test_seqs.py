from bit.modules.general import get_package_path
import bit.modules.seqs as seqs

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
