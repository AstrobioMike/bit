"""
gen-metagenome mutation layer.

Assigns a per-genome mutation rate and applies it, producing mutated genome fastas
and a per-genome mutation summary for the truth table. Drives bit's own
`mutate_seq` per record so behavior matches `bit mutate-seqs` exactly.

modes:
  off          - no mutation; genomes used as-is (returns rate 0.0, no new files)
  uniform      - every genome mutated at the same --mutation-rate
  distributed  - each genome drawn from [min,max] (uniform draw) so genomes sit at
                 varied ANI from their reference

Returns per-genome records: accession, mutation_rate, num_substitutions,
num_indels, ... aggregated across that genome's contigs.
"""
import os
import numpy as np
from Bio import SeqIO  # type: ignore
from bit.modules.seqs import mutate_seq


NT_SUBS = ['A', 'T', 'C', 'G']


def assign_rates(accessions, mode="off", mutation_rate=0.01,
                 rate_min=0.001, rate_max=0.05, seed=None):
    """ return dict accession -> rate, per the chosen mode. """
    rng = np.random.default_rng(seed)
    if mode == "off":
        return {a: 0.0 for a in accessions}
    if mode == "uniform":
        return {a: float(mutation_rate) for a in accessions}
    if mode == "distributed":
        return {a: float(rng.uniform(rate_min, rate_max)) for a in accessions}
    raise ValueError(f"mutation mode must be off/uniform/distributed, got '{mode}'")


def mutate_genome(in_fasta, out_fasta, rate, ti_tv_ratio=1.0, indel_rate=0.0,
                  seed=None):
    """
    Mutate every contig in in_fasta at `rate`, writing out_fasta. Returns a dict
    of aggregated counts across contigs. rate==0 copies sequences unchanged.
    Seeding: mutate_seq uses the global `random`; caller seeds once upstream for
    reproducibility across the run.
    """
    agg = dict(genome_size=0, num_substitutions=0, num_transitions=0,
               num_transversions=0, num_indels=0, num_insertions=0,
               num_deletions=0, num_total_changes=0)

    with open(in_fasta) as fin, open(out_fasta, "w") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            if rate > 0:
                (seq, tot, subs, ti, tv, indels, ins, dels) = mutate_seq(
                    rec.seq, "NT", NT_SUBS, rate, ti_tv_ratio, indel_rate)
            else:
                seq = str(rec.seq)
                tot = subs = ti = tv = indels = ins = dels = 0
            fout.write(f">{rec.id}\n{seq}\n")
            agg["genome_size"] += len(seq)
            agg["num_substitutions"] += subs
            agg["num_transitions"] += ti
            agg["num_transversions"] += tv
            agg["num_indels"] += indels
            agg["num_insertions"] += ins
            agg["num_deletions"] += dels
            agg["num_total_changes"] += tot
    return agg


def run_mutation(genome_paths, rates, out_dir, ti_tv_ratio=1.0, indel_rate=0.0,
                 seed=None, progress=None):
    """
    genome_paths: dict accession -> input fasta path
    rates:        dict accession -> rate
    Returns dict accession -> {mutated_fasta, rate, **agg_counts}.
    When all rates are 0 (mode off), input paths are returned unchanged.
    """
    import random
    if seed is not None:
        random.seed(seed)

    os.makedirs(out_dir, exist_ok=True)
    results = {}
    for acc, in_fasta in genome_paths.items():
        rate = rates.get(acc, 0.0)
        if rate == 0:
            results[acc] = dict(mutated_fasta=in_fasta, rate=0.0)
            if progress:
                progress.update(1)
            continue
        out_fasta = os.path.join(out_dir, f"{acc}.mutated.fasta")
        agg = mutate_genome(in_fasta, out_fasta, rate, ti_tv_ratio, indel_rate)
        results[acc] = dict(mutated_fasta=out_fasta, rate=rate, **agg)
        if progress:
            progress.update(1)
    return results
