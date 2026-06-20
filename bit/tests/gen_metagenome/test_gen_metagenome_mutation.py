import os
import numpy as np
import pytest
from Bio import SeqIO  # type: ignore

from bit.modules.gen_metagenome.mutation import (
    assign_rates,
    mutate_genome,
    run_mutation,
)


ACCS = ["GCF_A", "GCF_B", "GCF_C"]


@pytest.fixture
def genomes(tmp_path):
    """ three small 2-contig genomes with known random sequence. """
    rng = np.random.default_rng(0)
    paths = {}
    for acc in ACCS:
        p = tmp_path / f"{acc}.fasta"
        with open(p, "w") as f:
            for c in range(2):
                seq = "".join(rng.choice(list("ACGT"), size=10_000))
                f.write(f">{acc}_contig{c}\n{seq}\n")
        paths[acc] = str(p)
    return paths


# ─── rate assignment ───────────────────────────────────────────────────────

def test_assign_rates_off():
    assert all(v == 0.0 for v in assign_rates(ACCS, mode="off").values())


def test_assign_rates_uniform():
    rates = assign_rates(ACCS, mode="uniform", mutation_rate=0.02)
    assert all(v == 0.02 for v in rates.values())


def test_assign_rates_distributed_varies_in_range():
    rates = assign_rates(ACCS, mode="distributed", rate_min=0.01, rate_max=0.05, seed=1)
    assert all(0.01 <= v <= 0.05 for v in rates.values())
    assert len(set(rates.values())) == len(ACCS)   # distinct per genome


def test_assign_rates_invalid_mode_raises():
    with pytest.raises(ValueError):
        assign_rates(ACCS, mode="bogus")


# ─── mutation application ──────────────────────────────────────────────────

def test_uniform_observed_rate_matches_requested(genomes, tmp_path):
    rates = assign_rates(ACCS, mode="uniform", mutation_rate=0.03)
    results = run_mutation(genomes, rates, str(tmp_path / "out"), seed=42)
    for acc in ACCS:
        r = results[acc]
        observed = r["num_total_changes"] / r["genome_size"]
        assert abs(observed - 0.03) < 0.005


def test_distributed_observed_tracks_assigned(genomes, tmp_path):
    rates = assign_rates(ACCS, mode="distributed", rate_min=0.01, rate_max=0.08, seed=7)
    results = run_mutation(genomes, rates, str(tmp_path / "out"), seed=42)
    for acc in ACCS:
        r = results[acc]
        observed = r["num_total_changes"] / r["genome_size"]
        assert abs(observed - rates[acc]) < 0.006


def test_indels_split_out(genomes, tmp_path):
    results = run_mutation({"GCF_A": genomes["GCF_A"]}, {"GCF_A": 0.05},
                           str(tmp_path / "out"), indel_rate=0.2, seed=42)
    r = results["GCF_A"]
    frac_indel = r["num_indels"] / (r["num_substitutions"] + r["num_indels"])
    assert 0.12 < frac_indel < 0.28
    assert r["num_insertions"] + r["num_deletions"] == r["num_indels"]


def test_off_mode_passes_through_input_paths(genomes, tmp_path):
    results = run_mutation(genomes, assign_rates(ACCS, mode="off"),
                           str(tmp_path / "out_off"))
    for acc in ACCS:
        assert results[acc]["mutated_fasta"] == genomes[acc]
        assert results[acc]["rate"] == 0.0


def test_mutate_genome_writes_all_contigs(genomes, tmp_path):
    out = tmp_path / "mutated.fasta"
    agg = mutate_genome(genomes["GCF_A"], str(out), rate=0.02, seed=1)
    recs = list(SeqIO.parse(str(out), "fasta"))
    assert len(recs) == 2                      # both contigs preserved
    # substitutions-only: length unchanged, agg size == sum of written records
    assert agg["genome_size"] == sum(len(r.seq) for r in recs)
    assert agg["num_total_changes"] > 0


def test_mutate_genome_size_reflects_indels(genomes, tmp_path):
    # with indels, the mutated genome size must differ from the original and
    # must equal the actual written sequence length (regression for used_genome_size)
    out = tmp_path / "indel.fasta"
    agg = mutate_genome(genomes["GCF_A"], str(out), rate=0.05, indel_rate=0.5, seed=3)
    written = sum(len(r.seq) for r in SeqIO.parse(str(out), "fasta"))
    assert agg["genome_size"] == written            # agg tracks mutated length
    assert agg["genome_size"] != 20_000             # indels changed the length
    assert agg["genome_size"] == 20_000 + agg["num_insertions"] - agg["num_deletions"]


def test_rate_zero_copies_sequence_unchanged(genomes, tmp_path):
    out = tmp_path / "copy.fasta"
    agg = mutate_genome(genomes["GCF_A"], str(out), rate=0.0)
    orig = [str(r.seq) for r in SeqIO.parse(genomes["GCF_A"], "fasta")]
    copied = [str(r.seq) for r in SeqIO.parse(str(out), "fasta")]
    assert orig == copied
    assert agg["num_total_changes"] == 0
