//! Fold-back reads and read pairs from whole-genome amplification (PTA/MDA) for
//! `bit filter-foldbacks`.
//!
//! An amplified molecule that switched to the opposite strand nearby reads as a stretch and then
//! its own reverse complement (a hairpin). Both arms are real sequence; the join between them is
//! not, and whatever sits between the paired arms (often some untemplated bases) can't be trusted.
//!
//! Within one read, a fold shows as k-mers whose reverse complement is also in the read, all
//! agreeing on where the read turns (positions i and j with i + j the same). A pair whose fold is
//! in the unsequenced middle of the fragment shows as read 2 holding stretches of read 1 in the
//! same orientation: normal mates are on opposite strands, so all they can share (when the
//! fragment is shorter than the two reads) is read 1 against read 2's reverse complement.

use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// (position, forward code, reverse-complement code) of every k-mer of only A, C, G, T (k <= 31)
fn kmers(seq: &[u8], k: usize) -> Vec<(usize, u64, u64)> {
    let mut out = Vec::with_capacity(seq.len());
    if k == 0 || seq.len() < k {
        return out;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k - 1) as u64;
    let (mut fw, mut rc, mut valid) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        match code(b) {
            Some(c) => {
                fw = ((fw << 2) | c) & mask;
                rc = (rc >> 2) | ((3 - c) << shift);
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= k {
            out.push((i + 1 - k, fw, rc));
        }
    }
    out
}

/// a fold within one read: (end of the left arm's paired stretch, start of the right arm's paired
/// stretch, k-mer pairs agreeing on it), or None
fn fold(seq: &[u8], k: usize, min_shared: usize, tolerance: i64) -> Option<(usize, usize, usize)> {
    let km = kmers(seq, k);
    if km.len() < 2 * min_shared {
        return None;
    }
    let mut at: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for &(p, f, _) in &km {
        at.entry(f).or_default().push(p);
    }
    // (i + j, i, j) for each k-mer at i whose reverse complement is at j > i
    let mut pairs: Vec<(i64, usize, usize)> = Vec::new();
    for &(j, _, r) in &km {
        if let Some(ps) = at.get(&r) {
            for &i in ps {
                if i < j {
                    pairs.push(((i + j) as i64, i, j));
                }
            }
        }
    }
    if pairs.len() < min_shared {
        return None;
    }
    pairs.sort_unstable();
    let (mut best, mut best_at, mut s) = (0usize, 0usize, 0usize);
    for e in 0..pairs.len() {
        while pairs[e].0 - pairs[s].0 > 2 * tolerance {
            s += 1;
        }
        if e + 1 - s > best {
            best = e + 1 - s;
            best_at = s;
        }
    }
    if best < min_shared {
        return None;
    }
    let chosen = &pairs[best_at..best_at + best];
    let left_end = chosen.iter().map(|x| x.1).max().unwrap() + k;
    let right_start = chosen.iter().map(|x| x.2).min().unwrap();
    Some((left_end, right_start, best))
}

/// (start, end) of the read to keep: the longer arm of each fold, cut back to the end of its
/// paired stretch (so the join and anything unpaired beside it go), again on what's kept; and
/// how many folds were cut
fn trim(seq: &[u8], k: usize, min_shared: usize, tolerance: i64) -> (usize, usize, u8) {
    let (mut start, mut end, mut folds) = (0usize, seq.len(), 0u8);
    while folds < 4 {
        let Some((left_end, right_start, _)) = fold(&seq[start..end], k, min_shared, tolerance)
        else {
            break;
        };
        let left = left_end.min(end - start);
        let right = (end - start).saturating_sub(right_start);
        if left >= right {
            end = start + left;
        } else {
            start += right_start;
        }
        folds += 1;
    }
    (start, end, folds)
}

/// distinct k-mers of read 2 found in read 1 in the same orientation, and reverse complemented
fn shared(r1: &[u8], r2: &[u8], k: usize) -> (usize, usize) {
    let one: FxHashSet<u64> = kmers(r1, k).into_iter().map(|x| x.1).collect();
    let (mut same, mut rc) = (FxHashSet::default(), FxHashSet::default());
    for (_, f, r) in kmers(r2, k) {
        if one.contains(&f) {
            same.insert(f);
        }
        if one.contains(&r) {
            rc.insert(r);
        }
    }
    (same.len(), rc.len())
}

/// For each read (single-end) or pair: (read 1 start, end, folds cut; read 2 start, end, folds
/// cut; distinct k-mers read 2 shares with read 1 in the same orientation; and reverse
/// complemented). Read 2's fields are 0 for single-end input.
#[pyfunction]
#[pyo3(signature = (r1, r2=None, k=17, min_shared=20, tolerance=20))]
#[allow(clippy::type_complexity)]
pub fn scan_foldbacks(
    py: Python<'_>,
    r1: Vec<PyBackedStr>,
    r2: Option<Vec<PyBackedStr>>,
    k: usize,
    min_shared: usize,
    tolerance: i64,
) -> PyResult<Vec<(u32, u32, u8, u32, u32, u8, u32, u32)>> {
    if !(5..=31).contains(&k) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be between 5 and 31",
        ));
    }
    if let Some(r2) = &r2
        && r2.len() != r1.len()
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "read 1 and read 2 batches differ in length",
        ));
    }
    Ok(py.detach(|| {
        (0..r1.len())
            .into_par_iter()
            .map(|i| {
                let a = r1[i].as_bytes();
                let (s1, e1, f1) = trim(a, k, min_shared, tolerance);
                match &r2 {
                    None => (s1 as u32, e1 as u32, f1, 0, 0, 0, 0, 0),
                    Some(r2) => {
                        let b = r2[i].as_bytes();
                        let (s2, e2, f2) = trim(b, k, min_shared, tolerance);
                        let (same, rc) = shared(a, b, k);
                        (
                            s1 as u32,
                            e1 as u32,
                            f1,
                            s2 as u32,
                            e2 as u32,
                            f2,
                            same as u32,
                            rc as u32,
                        )
                    }
                }
            })
            .collect()
    }))
}

// Long reads (`--hifi`): the same fold finding and cutting as `shifia --pta`
// (github.com/AstrobioMike/shifia, rust/src/pta.rs), on the same sampled anchors (see
// anchors.rs), so the two cut reads at the same places. A long read can hold whole arms and the
// unpaired middle between them (PTA's is ~50-400 bp), and both arms of a fold share many anchors.
// Unlike short reads, a long read is mostly NOT fold, so a fold also has to pair a fair share of
// the anchors on the read's shorter side of it (min_fraction): a genomic inverted repeat inside
// long flanks pairs only a few.

/// one fold in seq: (end of the left arm's paired stretch, start of the right arm's paired
/// stretch, unpaired middle), from the pairs of anchors seen both ways that agree on the fold
fn long_fold(
    seq: &[u8],
    k: usize,
    threshold: u64,
    min_pairs: usize,
    min_fraction: f64,
    tolerance: i64,
) -> Option<(usize, usize, i64)> {
    let mut by_hash: FxHashMap<u64, Vec<(i64, bool)>> = FxHashMap::default();
    let mut positions: Vec<i64> = Vec::new();
    crate::anchors::for_each_anchor(seq, k, threshold, |p, h, f| {
        by_hash.entry(h).or_default().push((p as i64, f));
        positions.push(p as i64);
    });
    if positions.len() < 2 * min_pairs {
        return None;
    }
    // pairs (lo, hi) of one k-mer seen both ways; their fold is (lo + hi + k) / 2, kept doubled
    let mut pairs: Vec<(i64, i64, i64)> = Vec::new();
    for occ in by_hash.values() {
        if occ.len() < 2 {
            continue;
        }
        for &(a, fa) in occ {
            if !fa {
                continue;
            }
            for &(b, fb) in occ {
                if fb {
                    continue;
                }
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                pairs.push((lo + hi + k as i64, lo, hi));
            }
        }
    }
    if pairs.len() < min_pairs {
        return None;
    }
    pairs.sort_unstable();
    // the most agreeing window of folds (2 * tolerance wide, in doubled units)
    let w = 2 * tolerance;
    let (mut best, mut best_at, mut j) = (0usize, 0usize, 0usize);
    for i in 0..pairs.len() {
        while pairs[i].0 - pairs[j].0 > w {
            j += 1;
        }
        if i + 1 - j > best {
            best = i + 1 - j;
            best_at = j;
        }
    }
    if best < min_pairs {
        return None;
    }
    let chosen = &pairs[best_at..best_at + best];
    let fold2 = chosen[best / 2].0;
    let left_end = chosen.iter().map(|x| x.1).max().unwrap() + k as i64;
    let right_start = chosen.iter().map(|x| x.2).min().unwrap();
    // enough of the shorter side's anchors paired (a stray inverted repeat pairs a few)
    let shorter = if fold2 <= seq.len() as i64 {
        positions
            .iter()
            .filter(|&&p| 2 * (p + k as i64) <= fold2)
            .count()
    } else {
        positions.iter().filter(|&&p| 2 * p >= fold2).count()
    };
    if (best as f64) < min_fraction * shorter.max(1) as f64 {
        return None;
    }
    let middle = (right_start - left_end).max(0);
    Some((
        left_end.max(0) as usize,
        right_start.max(0) as usize,
        middle,
    ))
}

/// (start, end) of the long read to keep, folds cut, a sharp fold left whole (0/1, only with
/// min_middle > 0), and the first fold's unpaired middle (-1 for no fold found)
#[allow(clippy::too_many_arguments)]
fn long_trim(
    b: &[u8],
    k: usize,
    threshold: u64,
    min_middle: i64,
    min_pairs: usize,
    min_fraction: f64,
    tolerance: i64,
    max_folds: usize,
) -> (u32, u32, u8, u8, i32) {
    let (mut start, mut end) = (0usize, b.len());
    let (mut cuts, mut sharp, mut first) = (0u8, 0u8, -1i32);
    for _ in 0..max_folds {
        let Some((left_end, right_start, middle)) = long_fold(
            &b[start..end],
            k,
            threshold,
            min_pairs,
            min_fraction,
            tolerance,
        ) else {
            break;
        };
        if first < 0 {
            first = middle as i32;
        }
        if middle < min_middle {
            sharp = 1;
            break;
        }
        let left = left_end.min(end - start);
        let right = (end - start).saturating_sub(right_start);
        if left >= right {
            end = start + left;
        } else {
            start += right_start;
        }
        cuts += 1;
    }
    (start as u32, end as u32, cuts, sharp, first)
}

/// For each long read: (start, end) of what to keep, folds cut, a sharp fold left whole (0/1),
/// and the first fold's unpaired middle (-1 for no fold found). Anchors are k-mers sampled at
/// about 1 in `density`; defaults are shifia's.
#[pyfunction]
#[pyo3(signature = (seqs, k=21, density=4, min_middle=0, min_pairs=5, min_fraction=0.25,
                    tolerance=50, max_folds=4))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn trim_long_folds(
    py: Python<'_>,
    seqs: Vec<PyBackedStr>,
    k: usize,
    density: u64,
    min_middle: i64,
    min_pairs: usize,
    min_fraction: f64,
    tolerance: i64,
    max_folds: usize,
) -> PyResult<Vec<(u32, u32, u8, u8, i32)>> {
    if !(5..=31).contains(&k) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be between 5 and 31",
        ));
    }
    let threshold = crate::anchors::threshold(density);
    Ok(py.detach(|| {
        seqs.par_iter()
            .map(|s| {
                long_trim(
                    s.as_bytes(),
                    k,
                    threshold,
                    min_middle,
                    min_pairs,
                    min_fraction,
                    tolerance,
                    max_folds,
                )
            })
            .collect()
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 32) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn kmer_codes_match_revcomp() {
        let s = random(60, 1);
        let r = rc(&s);
        let a = kmers(&s, 11);
        let b = kmers(&r, 11);
        // the reverse complement of the first k-mer is the last k-mer of the reverse complement
        assert_eq!(a[0].2, b[b.len() - 1].1);
    }

    #[test]
    fn fold_in_read_cut_to_longer_arm() {
        let t = random(400, 2);
        // 100 bp first arm to 100, back from 60 (40 bp unpaired), down to 10
        let mut read = t[..100].to_vec();
        read.extend(rc(&t[10..60]));
        let (s, e, f) = trim(&read, 17, 20, 20);
        assert_eq!(f, 1);
        assert_eq!(s, 0);
        assert!((55..=65).contains(&e), "{e}");
    }

    #[test]
    fn plain_read_untouched() {
        let r = random(150, 3);
        assert_eq!(trim(&r, 17, 20, 20), (0, 150, 0));
    }

    #[test]
    fn pair_from_fold_shares_same_orientation() {
        let t = random(800, 4);
        // fragment: t[100..400] then back over t[150..350] reverse complemented
        let mut frag = t[100..400].to_vec();
        frag.extend(rc(&t[150..350]));
        let r1 = frag[..150].to_vec();
        let r2 = rc(&frag[frag.len() - 150..]);
        let (same, rcs) = shared(&r1, &r2, 17);
        assert!(same >= 20 && rcs == 0, "{same} {rcs}");
        // a normal pair from the same template shares nothing
        let r2n = rc(&t[300..450]);
        assert_eq!(shared(&t[100..250], &r2n, 17), (0, 0));
    }

    fn long_cut(read: &[u8]) -> (u32, u32, u8, u8, i32) {
        long_trim(read, 21, crate::anchors::threshold(4), 0, 5, 0.25, 50, 4)
    }

    #[test]
    fn long_fold_cut_to_longer_arm() {
        let t = random(6000, 5);
        // 2,000 bp to 2,000, back from 1,800 (200 bp unpaired middle) down to 1,000
        let mut read = t[..2000].to_vec();
        read.extend(rc(&t[1000..1800]));
        let (s, e, f, _, middle) = long_cut(&read);
        assert_eq!((s, f), (0, 1));
        assert!((1750..=1800).contains(&e), "{e}");
        assert!((150..=250).contains(&middle), "{middle}");
    }

    #[test]
    fn long_fold_keeps_the_longer_right_arm() {
        let t = random(6000, 6);
        // 600 bp out to 1,600, back from 1,500 down to -1,000 (past where the read began)
        let mut read = t[1000..1600].to_vec();
        read.extend(rc(&t[0..1500]));
        let (s, e, f, _, _) = long_cut(&read);
        assert_eq!(f, 1);
        assert_eq!(e as usize, read.len());
        // the kept piece starts where the right arm's paired stretch starts: after the 600 bp
        // left arm, whose last 100 bp (1,500-1,600) are the unpaired middle
        assert!((595..=650).contains(&s), "{s}");
    }

    #[test]
    fn long_read_without_fold_untouched() {
        let r = random(12000, 7);
        assert_eq!(long_cut(&r), (0, 12000, 0, 0, -1));
    }

    #[test]
    fn genomic_inverted_repeat_in_long_flanks_untouched() {
        let t = random(9000, 8);
        // 3 kb, a 300 bp arm, 150 bp spacer, the arm reverse complemented, 3 kb
        let mut read = t[..3450].to_vec();
        read.extend(rc(&t[3000..3300]));
        read.extend(&t[5000..8000]);
        assert_eq!(long_cut(&read).2, 0);
    }
}
