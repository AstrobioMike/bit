//! Sampled canonical k-mers ("anchors") of a read, for the long-read (`--hifi`) mode of
//! `bit filter-foldbacks` and for `bit normalize-reads`.
//!
//! The same hashing as shifia's (github.com/AstrobioMike/shifia, rust/src/hash.rs), so the
//! long-read mode cuts reads exactly where `shifia --pta` does: a polynomial rolling hash
//! mod 2^64 over A=1, C=2, G=3, T=4 for both the k-mer and its reverse complement, each passed
//! through a splitmix64 finalizer; the canonical hash is the smaller of the two, and a k-mer is
//! an anchor when that is below a threshold (about 1 in `density` k-mers, see threshold()).
//! K-mers holding anything but A, C, G or T are skipped.

const BASE: u64 = 0x100000001B3;

const fn inverse(b: u64) -> u64 {
    // Newton iteration for the inverse of an odd number mod 2^64
    let mut x = b;
    let mut i = 0;
    while i < 6 {
        x = x.wrapping_mul(2u64.wrapping_sub(b.wrapping_mul(x)));
        i += 1;
    }
    x
}

const BASE_INV: u64 = inverse(BASE);

#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xBF58476D1CE4E5B9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94D049BB133111EB);
    x ^= x >> 31;
    x
}

#[inline(always)]
fn is_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[inline(always)]
fn code(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 3,
        b'T' | b't' => 4,
        _ => 1,
    }
}

/// hashes below this are anchors, about 1 in `density` k-mers (the canonical hash is the
/// smaller of two roughly uniform hashes, so P(canonical < t) ~ 2t)
pub fn threshold(density: u64) -> u64 {
    u64::MAX / (2 * density.max(1))
}

/// calls `f(position, canonical hash, is_forward)` for every anchor of `seq`, in order
pub fn for_each_anchor<F: FnMut(usize, u64, bool)>(seq: &[u8], k: usize, threshold: u64, mut f: F) {
    let n = seq.len();
    if k == 0 || n < k {
        return;
    }
    let mut bk1: u64 = 1;
    for _ in 0..k - 1 {
        bk1 = bk1.wrapping_mul(BASE);
    }
    let (mut fwd, mut rev, mut p) = (0u64, 0u64, 1u64);
    let mut last_bad: isize = -1; // the last non-ACGT base seen
    for (i, &b) in seq[..k].iter().enumerate() {
        if !is_base(b) {
            last_bad = i as isize;
        }
        let c = code(b);
        fwd = fwd.wrapping_mul(BASE).wrapping_add(c);
        rev = rev.wrapping_add((5 - c).wrapping_mul(p));
        p = p.wrapping_mul(BASE);
    }
    let mut emit = |s: usize, fwd: u64, rev: u64| {
        let (mf, mr) = (mix(fwd), mix(rev));
        let h = mf.min(mr);
        if h < threshold {
            f(s, h, mf <= mr);
        }
    };
    if last_bad < 0 {
        emit(0, fwd, rev);
    }
    for (s, (&out, &inc)) in seq[..n - k].iter().zip(&seq[k..]).enumerate() {
        fwd = fwd
            .wrapping_sub(code(out).wrapping_mul(bk1))
            .wrapping_mul(BASE)
            .wrapping_add(code(inc));
        rev = rev
            .wrapping_sub(5 - code(out))
            .wrapping_mul(BASE_INV)
            .wrapping_add((5 - code(inc)).wrapping_mul(bk1));
        if !is_base(inc) {
            last_bad = (s + k) as isize;
        }
        // the k-mer starting at s + 1 holds a non-ACGT base
        if last_bad > s as isize {
            continue;
        }
        emit(s + 1, fwd, rev);
    }
}

/// the anchors of `seq` as (position, hash), in order, into `out` (cleared first)
pub fn anchors_into(seq: &[u8], k: usize, threshold: u64, out: &mut Vec<(u32, u64)>) {
    out.clear();
    for_each_anchor(seq, k, threshold, |p, h, _| out.push((p as u32, h)));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    }

    #[test]
    fn canonical_on_both_strands() {
        let s = b"ACGTTGCAAGGCTTACCGATAGGCTAACGTTAGCATGCAAGTCCGATTACGATCGGATCAAT";
        let rc = revcomp(s);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for_each_anchor(s, 21, u64::MAX, |_, h, _| a.push(h));
        for_each_anchor(&rc, 21, u64::MAX, |_, h, _| b.push(h));
        b.reverse();
        assert_eq!(a, b);
        assert_eq!(a.len(), s.len() - 20);
    }

    #[test]
    fn skips_kmers_with_n() {
        let s = b"ACGTTGCAAGGCTTACCGATAGNCTAACGTTAGCATGCAAGTCCGATTACG";
        let mut at = Vec::new();
        for_each_anchor(s, 21, u64::MAX, |p, _, _| at.push(p));
        assert!(at.iter().all(|&p| p + 21 <= 22 || p > 22));
    }
}
