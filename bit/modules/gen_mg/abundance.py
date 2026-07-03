"""
gen-mg abundance/coverage layer

Given a normalized genome table (needs 'accession', 'used_genome_size', and the pinned
override columns), assign each genome a relative abundance, a target coverage, and
a read count, sticking to:
  - mode: 'relative' (distribution describes rel. abundance; needs total_reads) or
          'coverage' (distribution describes per-genome fold-coverage directly)
  - dist: 'lognormal' (sigma) or 'even'
  - per-genome pins from the user TSV (pinned_rel_abundance / pinned_coverage)

Pure logic: returns a new dataframe with added columns
  assigned_rel_abundance, assigned_coverage, assigned_reads
and a warnings list
"""
import numpy as np # type: ignore
import pandas as pd # type: ignore

PAIRED_TYPES = {"paired-end"}


def _draw_weights(n, dist, sigma=1.0, rng=None):
    """ draw n positive weights from the chosen distribution (unnormalized). """
    rng = rng or np.random.default_rng()
    if dist == "even":
        return np.ones(n)
    if dist == "lognormal":
        return rng.lognormal(mean=0.0, sigma=sigma, size=n)
    raise ValueError(f"unknown distribution '{dist}'")


def _reads_for_coverage(coverage, genome_size, read_length, read_type, fragment_size):
    """ reads needed to reach `coverage` over genome_size, matching gen-reads' own math. """
    if read_type in PAIRED_TYPES:
        bases_per_fragment = min(2 * read_length, fragment_size)
        fragments = round(coverage * genome_size / bases_per_fragment)
        return 2 * fragments
    return round(coverage * genome_size / read_length)


def _apportion_exact(rel, total_reads, paired=False):
    """
    Distribute total_reads across genomes in proportion to `rel` so the per-genome
    counts sum EXACTLY to total_reads (largest-remainder / Hamilton method), rather
    than rounding each independently (which drifts off the target).

    For paired-end, counts must be even (reads come in pairs), so apportionment is
    done in fragment units (total_reads // 2) and then doubled — the per-genome
    counts stay even and still sum to total_reads (assuming total_reads is even,
    which gen_reads guarantees for paired mode).
    """
    rel = np.asarray(rel, dtype=float)
    unit = 2 if paired else 1
    target = int(total_reads) // unit              # apportion in fragments (paired) or reads

    exact = rel * target
    floor = np.floor(exact).astype(int)
    leftover = target - int(floor.sum())           # how many units the floors left unallocated

    if leftover > 0:
        # hand the leftover units to the genomes with the largest fractional parts
        frac = exact - floor
        order = np.argsort(-frac)                  # largest remainder first
        for i in order[:leftover]:
            floor[i] += 1

    return floor * unit


def _coverage_from_reads(reads, genome_size, read_length):
    """ realized fold-coverage from a read count (bases sequenced / genome size). """
    if pd.isna(genome_size) or genome_size in (0, None):
        return float("nan")
    return reads * read_length / genome_size


def assign_abundance(df, mode="relative", dist="lognormal", sigma=1.0,
                     total_reads=1_000_000, read_length=150, read_type="paired-end",
                     fragment_size=500, median_coverage=30.0, seed=None):
    """
    Returns (df_with_assignments, warnings).
    """
    warnings = []
    rng = np.random.default_rng(seed)
    out = df.copy().reset_index(drop=True)
    n = len(out)
    if n == 0:
        return out, ["No genomes to assign abundance to."]

    sizes = pd.to_numeric(out["used_genome_size"], errors="coerce")
    if sizes.isna().any():
        warnings.append(f"{int(sizes.isna().sum())} genome(s) have no usable genome size; "
                        "their coverage/reads can't be computed and will be NA.")

    if mode == "relative":
        # start from drawn weights, then apply pins
        weights = _draw_weights(n, dist, sigma=sigma, rng=rng)
        pinned = pd.to_numeric(out["pinned_rel_abundance"], errors="coerce")

        # genomes with a pinned rel-abundance take that fraction off the top;
        # remaining genomes split the leftover proportionally to their drawn weights
        pin_mask = pinned.notna()
        pin_total = float(pinned[pin_mask].sum()) if pin_mask.any() else 0.0
        free_mask = ~pin_mask.values
        n_free = int(free_mask.sum())
        if pin_total > 1.0:
            # over-sum: pins are renormalized down. (With --num-genomes this is
            # caught earlier in preflight; this branch remains for accessions-only
            # runs, where all or some genomes are pinned and none are drawn from
            # --num-genomes.) Word it for whether un-pinned genomes actually exist.
            if n_free > 0:
                warnings.append(f"Pinned relative abundances sum to {pin_total:.3f} "
                                f"(>1); they will be renormalized and the {n_free} "
                                "un-pinned genome(s) will get ~0.")
            else:
                warnings.append(f"Pinned relative abundances sum to {pin_total:.3f} "
                                "(>1); they will be renormalized to sum to 1.")
        elif pin_mask.any() and n_free == 0 and abs(pin_total - 1.0) > 1e-6:
            # all genomes pinned but the pins don't sum to 1: renormalize (treat
            # the pins as relative weights) and tell the user it was rescaled.
            warnings.append(f"All genomes have pinned relative abundances summing to "
                            f"{pin_total:.3f} (not 1); renormalizing them to sum to 1.")
        elif n_free > 0 and abs(pin_total - 1.0) <= 1e-6:
            # pins already fill the whole budget, so unpinned genomes get ~0
            warnings.append(f"Pinned relative abundances sum to ~1; the {n_free} "
                            "un-pinned genome(s) will get ~0 relative abundance.")
        rel = np.zeros(n)
        rel[pin_mask.values] = pinned[pin_mask].values
        leftover = max(0.0, 1.0 - pin_total)
        if free_mask.any() and leftover > 0:
            w = weights[free_mask]
            rel[free_mask] = leftover * (w / w.sum())
        # normalize defensively
        rel = rel / rel.sum()
        out["assigned_rel_abundance"] = rel

        reads = _apportion_exact(rel, total_reads, paired=(read_type == "paired-end"))
        out["assigned_reads"] = reads
        out["assigned_coverage"] = [
            _coverage_from_reads(r, s, read_length) for r, s in zip(reads, sizes)
        ]

    elif mode == "coverage":
        # distribution sets per-genome fold coverage, scaled so the median sits at
        # median_coverage; pins override. drawn weights center on 1.0 (even -> all
        # ones; lognormal -> median e^0 = 1), so multiplying by median_coverage
        # gives 'even' exactly median_coverage per genome, and 'lognormal' a spread
        # whose median is median_coverage (mean runs higher due to the right tail)
        weights = _draw_weights(n, dist, sigma=sigma, rng=rng)
        cov = weights * median_coverage
        pinned_cov = pd.to_numeric(out["pinned_coverage"], errors="coerce")
        cov = np.where(pinned_cov.notna().values, pinned_cov.values, cov)
        out["assigned_coverage"] = cov

        reads = np.array([
            _reads_for_coverage(c, s, read_length, read_type, fragment_size)
            if not pd.isna(s) else 0
            for c, s in zip(cov, sizes)
        ], dtype=int)
        out["assigned_reads"] = reads
        total = reads.sum()
        out["assigned_rel_abundance"] = reads / total if total else np.nan

    else:
        raise ValueError(f"mode must be 'relative' or 'coverage', got '{mode}'")

    return out, warnings


def write_coverage_tsv(df, fasta_path_col, out_path):
    """
    Write the per-fasta coverage TSV that `gen-reads -c <file>` consumes:
    column 1 = path to fasta, column 2 = target coverage. Skips NA-coverage rows.
    """
    with open(out_path, "w") as fh:
        for _, r in df.iterrows():
            cov = r["assigned_coverage"]
            if pd.isna(cov):
                continue
            fh.write(f"{r[fasta_path_col]}\t{cov:.6f}\n")
