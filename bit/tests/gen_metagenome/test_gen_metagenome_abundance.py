import numpy as np
import pandas as pd
import pytest

from bit.modules.gen_mg.abundance import (
    assign_abundance,
    write_coverage_tsv,
    _reads_for_coverage,
    _coverage_from_reads,
    _apportion_exact,
)


def base(sizes, pins_ra=None, pins_cov=None):
    n = len(sizes)
    return pd.DataFrame({
        "accession": [f"g{i}" for i in range(n)],
        "used_genome_size": sizes,
        "pinned_rel_abundance": pins_ra if pins_ra else [pd.NA] * n,
        "pinned_coverage": pins_cov if pins_cov else [pd.NA] * n,
    })


# ─── relative mode ─────────────────────────────────────────────────────────

def test_relative_even_equal_sizes_equal_everything():
    out, _ = assign_abundance(base([4_000_000] * 4), mode="relative", dist="even",
                              total_reads=4_000_000, seed=1)
    assert np.allclose(out["assigned_rel_abundance"], 0.25)
    assert out["assigned_reads"].sum() == 4_000_000
    assert np.allclose(out["assigned_coverage"], out["assigned_coverage"].iloc[0])


def test_relative_even_unequal_sizes_inverse_coverage():
    out, _ = assign_abundance(base([6_000_000, 4_000_000, 2_000_000, 1_000_000]),
                              mode="relative", dist="even", total_reads=4_000_000,
                              read_length=150, seed=1)
    assert np.allclose(out["assigned_rel_abundance"], 0.25)
    covs = out["assigned_coverage"].tolist()
    # tiny genome (6x smaller) gets ~6x the coverage of the big one
    assert covs[3] > covs[0]
    assert abs(covs[3] / covs[0] - 6.0) < 0.05


def test_relative_pinned_abundance_respected():
    out, _ = assign_abundance(base([4_000_000] * 4, pins_ra=[0.5, pd.NA, pd.NA, pd.NA]),
                              mode="relative", dist="even", total_reads=1_000_000, seed=1)
    assert abs(out["assigned_rel_abundance"].iloc[0] - 0.5) < 1e-9
    for i in (1, 2, 3):
        assert abs(out["assigned_rel_abundance"].iloc[i] - 0.5 / 3) < 1e-9


def test_lognormal_sigma_controls_spread():
    df = base([4_000_000] * 200)
    o1, _ = assign_abundance(df, mode="relative", dist="lognormal", sigma=0.5,
                             total_reads=10 ** 7, seed=2)
    o2, _ = assign_abundance(df, mode="relative", dist="lognormal", sigma=2.0,
                             total_reads=10 ** 7, seed=2)
    cv1 = o1["assigned_rel_abundance"].std() / o1["assigned_rel_abundance"].mean()
    cv2 = o2["assigned_rel_abundance"].std() / o2["assigned_rel_abundance"].mean()
    assert cv2 > cv1


def test_relative_ignores_median_coverage():
    # median_coverage only applies in coverage mode; relative output must be identical
    o1, _ = assign_abundance(base([4_000_000] * 5), mode="relative", dist="even",
                             total_reads=10_000, median_coverage=1.0, seed=3)
    o2, _ = assign_abundance(base([4_000_000] * 5), mode="relative", dist="even",
                             total_reads=10_000, median_coverage=999.0, seed=3)
    assert (o1["assigned_reads"].values == o2["assigned_reads"].values).all()


# ─── exact total reads (relative mode) ─────────────────────────────────────
# independent per-genome rounding used to drift off the requested total (e.g.
# 10,000,000 -> 9,999,998); largest-remainder apportionment must hit it exactly.

def test_relative_exact_total_uneven_sizes_paired():
    # uneven genome sizes => uneven per-genome reads => the case that used to drift
    out, _ = assign_abundance(base([6_000_000, 4_100_000, 2_300_000, 999_000, 5_555_555]),
                              mode="relative", dist="even", total_reads=10_000_000,
                              read_type="paired-end", seed=1)
    assert out["assigned_reads"].sum() == 10_000_000
    # paired => every per-genome count must be even (reads come in pairs)
    assert (out["assigned_reads"].values % 2 == 0).all()


def test_relative_exact_total_lognormal_many_genomes_paired():
    out, _ = assign_abundance(base([4_000_000] * 137), mode="relative",
                              dist="lognormal", sigma=1.5, total_reads=10_000_000,
                              read_type="paired-end", seed=4)
    assert out["assigned_reads"].sum() == 10_000_000
    assert (out["assigned_reads"].values % 2 == 0).all()


def test_relative_exact_total_single_end_odd_total():
    # single-end: no even constraint, exact odd totals must be preserved
    out, _ = assign_abundance(base([6_000_000, 4_100_000, 2_300_000]),
                              mode="relative", dist="even", total_reads=9_999_999,
                              read_type="single-end", seed=1)
    assert out["assigned_reads"].sum() == 9_999_999


def test_apportion_exact_helper():
    rng = np.random.default_rng(0)
    rel = rng.dirichlet(np.ones(50))
    paired = _apportion_exact(rel, 10_000_000, paired=True)
    assert paired.sum() == 10_000_000
    assert (paired % 2 == 0).all()
    single = _apportion_exact(rel, 9_999_999, paired=False)
    assert single.sum() == 9_999_999


# ─── coverage mode ─────────────────────────────────────────────────────────

def test_coverage_even_gives_median_coverage_to_every_genome():
    out, _ = assign_abundance(base([6_000_000, 2_000_000]), mode="coverage",
                              dist="even", median_coverage=30, read_length=150,
                              read_type="paired-end", fragment_size=500, seed=1)
    # even -> every genome gets exactly median_coverage
    assert np.allclose(out["assigned_coverage"], 30.0)
    # reads scale with genome size (3x bigger genome -> 3x reads at equal coverage)
    ratio = out["assigned_reads"].iloc[0] / out["assigned_reads"].iloc[1]
    assert abs(ratio - 3.0) < 0.01


def test_coverage_lognormal_median_tracks_median_coverage():
    # with lognormal, the MEDIAN of assigned coverage should sit near median_coverage
    # (the mean runs higher because of the right tail — that's expected)
    out, _ = assign_abundance(base([4_000_000] * 400), mode="coverage",
                              dist="lognormal", sigma=1.0, median_coverage=30, seed=5)
    med = np.median(out["assigned_coverage"])
    assert abs(med - 30) < 4                      # median near target
    assert out["assigned_coverage"].mean() > med  # mean pulled up by the tail


def test_coverage_mode_read_count_matches_gen_reads_formula():
    size, cov, rl, frag = 5_000_000, 10, 150, 500
    expected = 2 * round(cov * size / min(2 * rl, frag))
    df = base([size], pins_cov=[cov])
    out, _ = assign_abundance(df, mode="coverage", dist="even", read_length=rl,
                              read_type="paired-end", fragment_size=frag, seed=1)
    assert out["assigned_reads"].iloc[0] == expected


def test_coverage_pinned_overrides_distribution():
    # a pinned coverage wins over the median_coverage scaling; unpinned gets median_coverage
    out, _ = assign_abundance(base([4_000_000, 4_000_000], pins_cov=[100, pd.NA]),
                              mode="coverage", dist="even", median_coverage=30, seed=1)
    assert out["assigned_coverage"].iloc[0] == 100
    assert out["assigned_coverage"].iloc[1] == 30


# ─── edge cases ────────────────────────────────────────────────────────────

def test_na_genome_size_handled():
    out, warnings = assign_abundance(base([4_000_000, pd.NA]), mode="relative",
                                     dist="even", total_reads=1_000_000, seed=1)
    assert pd.isna(out["assigned_coverage"].iloc[1])
    assert any("genome size" in w for w in warnings)


def test_empty_input():
    out, warnings = assign_abundance(base([]), mode="relative", seed=1)
    assert len(out) == 0
    assert warnings


def test_invalid_mode_raises():
    with pytest.raises(ValueError):
        assign_abundance(base([4_000_000]), mode="bogus")


def test_invalid_dist_raises():
    with pytest.raises(ValueError):
        assign_abundance(base([4_000_000]), mode="relative", dist="bogus")
    # exponential was trimmed; lognormal + even are the supported shapes
    with pytest.raises(ValueError):
        assign_abundance(base([4_000_000]), mode="relative", dist="exponential")


# ─── helpers ───────────────────────────────────────────────────────────────

def test_reads_for_coverage_paired_and_single():
    # paired: 2 * round(cov*size / min(2*rl, frag))
    assert _reads_for_coverage(10, 1_000_000, 150, "paired-end", 500) == \
        2 * round(10 * 1_000_000 / min(300, 500))
    # single: round(cov*size / rl)
    assert _reads_for_coverage(10, 1_000_000, 150, "single-end", 500) == \
        round(10 * 1_000_000 / 150)


def test_coverage_from_reads_and_na():
    assert _coverage_from_reads(1000, 150_000, 150) == 1.0
    assert np.isnan(_coverage_from_reads(1000, pd.NA, 150))


def test_write_coverage_tsv_skips_na(tmp_path):
    df = pd.DataFrame({
        "fasta_path": ["a.fasta", "b.fasta"],
        "assigned_coverage": [12.5, float("nan")],
    })
    out = tmp_path / "cov.tsv"
    write_coverage_tsv(df, "fasta_path", str(out))
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("a.fasta\t12.5")
