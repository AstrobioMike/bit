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
}
