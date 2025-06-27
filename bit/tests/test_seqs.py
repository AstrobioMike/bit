import os
import pytest
from bit.modules.general import get_package_path
import bit.modules.seqs as seqs

test_fasta = get_package_path("tests/data/ez-screen-targets.fasta")

expected_gc_stats = [
    {"header": "yopE", "length": 660, "gc": 0.51},
    {"header": "yopK", "length": 549, "gc": 0.34},
]

def test_calc_gc_per_seq():

    gc_stats = seqs.calc_gc_per_seq(input_fasta=test_fasta)
    assert gc_stats == expected_gc_stats, "GC stats do not match expected values"
