//! The read-matching core of `bit match-reads`: a port of ExactMatcher, EditMatcher and
//! process_batch() in bit/modules/match_reads.py, following them step for step so results
//! are identical (the Python classes stay as the reference, and tests compare the two).
//! ReadMatcher.assign_files (lib.rs) also reads the input (readio.rs), keeps the counts and
//! writes the hits table (SortedRows, as SortedRowWriter in match_reads.py) and the reads.

use crate::edlib;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

/// a k-mer as a lookup key: packed 2 bits per base when it's all ACGT and fits in a u64
/// (the usual case, and much smaller to hold), otherwise its bytes
enum Key<'a> {
    Packed(u64),
    Bytes(&'a [u8]),
    Owned(Vec<u8>),
}

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn key(kmer: &[u8]) -> Key<'_> {
    if kmer.len() <= 32 {
        let mut v = 0u64;
        let mut ok = true;
        for &b in kmer {
            match code(b) {
                Some(c) => v = (v << 2) | c,
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Key::Packed(v);
        }
    }
    Key::Bytes(kmer)
}

/// as Python's revcomp(): ACGTacgt complemented, anything else left as is, then reversed
fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        o => o,
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// as Python's canon(): the lesser of the k-mer and its reverse complement
fn canon(kmer: &[u8]) -> Vec<u8> {
    let rc = revcomp(kmer);
    if kmer <= rc.as_slice() { kmer.to_vec() } else { rc }
}

/// the canonical k-mer's key without building any strings: for an all-ACGT k-mer of up to
/// 32 bases, the lesser of its packed form and its reverse complement's (A<C<G<T, so this is
/// the same choice as comparing the strings); otherwise canon()'s bytes
fn canon_key(kmer: &[u8]) -> Key<'static> {
    if kmer.len() <= 32 {
        let (mut fwd, mut rc) = (0u64, 0u64);
        let mut ok = true;
        for (i, &b) in kmer.iter().enumerate() {
            match code(b) {
                Some(c) => {
                    fwd = (fwd << 2) | c;
                    rc |= (3 - c) << (2 * i);
                }
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Key::Packed(fwd.min(rc));
        }
    }
    Key::Owned(canon(kmer))
}

// packed k-mers are held in this many hash maps, chosen by key, so the maps can be filled on
// several threads at once
const SHARDS: usize = 64;

fn shard(key: u64) -> usize {
    (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 58) as usize
}

/// k-mer -> V, split by key kind (see Key)
struct KmerMap<V> {
    packed: Vec<FxHashMap<u64, V>>,
    bytes: FxHashMap<Box<[u8]>, V>,
}

impl<V: Send> KmerMap<V> {
    /// from packed keys (each once) with their values, and the other keys' map
    fn build(packed: Vec<(u64, V)>, bytes: FxHashMap<Box<[u8]>, V>) -> Self {
        let mut buckets: Vec<Vec<(u64, V)>> = (0..SHARDS).map(|_| Vec::new()).collect();
        for (k, v) in packed {
            buckets[shard(k)].push((k, v));
        }
        let packed = buckets
            .into_par_iter()
            .map(|b| {
                let mut m = FxHashMap::with_capacity_and_hasher(b.len(), Default::default());
                m.extend(b);
                m
            })
            .collect();
        KmerMap { packed, bytes }
    }

    fn get(&self, kmer: &[u8]) -> Option<&V> {
        self.get_key(key(kmer))
    }

    fn get_key(&self, k: Key) -> Option<&V> {
        match k {
            Key::Packed(v) => self.packed[shard(v)].get(&v),
            Key::Bytes(b) => self.bytes.get(b),
            Key::Owned(b) => self.bytes.get(b.as_slice()),
        }
    }
}

/// k-mer -> its (sequence, position) entries, all held in one list (a list per k-mer took a
/// third of the time to build an index, and as long again to free)
struct PosIndex {
    map: KmerMap<(u64, u32)>, // (start in `entries`, how many)
    entries: Vec<(u32, u32)>,
}

impl PosIndex {
    /// from gather_kmers()'s output
    fn build(packed: &[(u64, u32, u32)], other: Vec<(Box<[u8]>, u32, u32)>) -> Self {
        let mut entries = Vec::with_capacity(packed.len() + other.len());
        let mut keys = Vec::new();
        for (km, run) in key_runs(packed) {
            keys.push((km, (entries.len() as u64, run.len() as u32)));
            entries.extend(run.iter().map(|&(_, i, p)| (i, p)));
        }
        let mut grouped: FxHashMap<Box<[u8]>, Vec<(u32, u32)>> = FxHashMap::default();
        for (km, i, p) in other {
            grouped.entry(km).or_default().push((i, p));
        }
        let mut bytes = FxHashMap::default();
        for (km, v) in grouped {
            bytes.insert(km, (entries.len() as u64, v.len() as u32));
            entries.extend(v);
        }
        PosIndex { map: KmerMap::build(keys, bytes), entries }
    }

    fn get(&self, kmer: &[u8]) -> Option<&[(u32, u32)]> {
        self.map
            .get(kmer)
            .map(|&(start, n)| &self.entries[start as usize..start as usize + n as usize])
    }
}

/// Python's `a % b` for b > 0 (never negative)
fn pymod(a: i64, b: i64) -> i64 {
    a.rem_euclid(b)
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || memchr::memmem::find(hay, needle).is_some()
}

/// sets of ref indices, interned so each distinct set is held once (as Python shares
/// identical frozensets)
struct SetPool {
    sets: Vec<Vec<u32>>,
    ids: FxHashMap<Vec<u32>, u32>,
    // the ids of one-sequence sets, by the sequence (u32::MAX: none yet), the usual case
    single: Vec<u32>,
}

impl SetPool {
    fn new() -> Self {
        SetPool {
            sets: Vec::new(),
            ids: FxHashMap::default(),
            single: Vec::new(),
        }
    }

    fn intern(&mut self, mut s: Vec<u32>) -> u32 {
        s.sort_unstable();
        s.dedup();
        self.intern_sorted(&s)
    }

    /// as intern(), for a set already sorted and without repeats (nothing is allocated if
    /// the set is already held)
    fn intern_sorted(&mut self, s: &[u32]) -> u32 {
        if let [i] = s {
            let i = *i as usize;
            if i < self.single.len() && self.single[i] != u32::MAX {
                return self.single[i];
            }
        }
        if let Some(&id) = self.ids.get(s) {
            return id;
        }
        let id = self.sets.len() as u32;
        self.sets.push(s.to_vec());
        self.ids.insert(s.to_vec(), id);
        if let [i] = s {
            let i = *i as usize;
            if self.single.len() <= i {
                self.single.resize(i + 1, u32::MAX);
            }
            self.single[i] = id;
        }
        id
    }
}

/// The k-mers of every sequence, gathered on the threads: `at(i)` gives sequence i's text and
/// the positions of its k-mers in it. K-mers that pack (see Key) come back as (key, sequence,
/// position) sorted by key, then sequence, then position; those that don't, as (key bytes,
/// sequence, position) in sequence then position order. The key is the canonical k-mer's with
/// `canonical`, else the k-mer's own. Either way, each k-mer's (sequence, position) pairs come
/// in the order a pass over the sequences in turn would meet them, so maps built from these
/// hold the same as ones filled k-mer by k-mer.
#[allow(clippy::type_complexity)]
fn gather_kmers(
    n_seqs: usize,
    k: usize,
    canonical: bool,
    at: impl Fn(usize) -> Option<(Vec<u8>, Vec<usize>)> + Sync,
) -> (Vec<(u64, u32, u32)>, Vec<(Box<[u8]>, u32, u32)>) {
    let per_seq: Vec<(Vec<(u64, u32, u32)>, Vec<(Box<[u8]>, u32, u32)>)> = (0..n_seqs)
        .into_par_iter()
        .map(|i| {
            let (mut packed, mut other) = (Vec::new(), Vec::new());
            if let Some((text, positions)) = at(i) {
                packed.reserve(positions.len());
                for p in positions {
                    let kmer = &text[p..p + k];
                    let key = if canonical { canon_key(kmer) } else { key(kmer) };
                    match key {
                        Key::Packed(v) => packed.push((v, i as u32, p as u32)),
                        Key::Bytes(b) => other.push((b.into(), i as u32, p as u32)),
                        Key::Owned(b) => other.push((b.into_boxed_slice(), i as u32, p as u32)),
                    }
                }
            }
            (packed, other)
        })
        .collect();
    let total = per_seq.iter().map(|(p, _)| p.len()).sum();
    let mut packed = Vec::with_capacity(total);
    let mut other = Vec::new();
    for (p, o) in per_seq {
        packed.extend(p);
        other.extend(o);
    }
    packed.par_sort_unstable();
    (packed, other)
}

/// runs of equal keys in gather_kmers()'s sorted packed k-mers: (key, its entries)
fn key_runs(sorted: &[(u64, u32, u32)]) -> impl Iterator<Item = (u64, &[(u64, u32, u32)])> {
    sorted.chunk_by(|a, b| a.0 == b.0).map(|run| (run[0].0, run))
}

/// the sequences (sorted, each once) of a run of gather_kmers() entries, into `buf`
fn run_seqs<'a>(run: &[(u64, u32, u32)], buf: &'a mut Vec<u32>) -> &'a [u32] {
    buf.clear();
    for &(_, i, _) in run {
        if buf.last() != Some(&i) {
            buf.push(i);
        }
    }
    buf
}

pub struct Refs {
    seqs: Vec<Vec<u8>>,
    circular: bool,
}

impl Refs {
    /// _RefMatcher._matches_at(): `read` exactly matches ref i beginning at `start`
    fn matches_at(&self, i: usize, read: &[u8], start: i64) -> bool {
        let r = &self.seqs[i];
        let starts_with = |hay: &[u8], needle: &[u8], at: usize| -> bool {
            at <= hay.len() && hay[at..].starts_with(needle)
        };
        if start < 0 {
            // (only reached for linear refs; Python's startswith() with a negative start
            // counts from the end, but callers never pass one)
            return false;
        }
        let start = start as usize;
        if !self.circular {
            return starts_with(r, read, start);
        }
        let (l, n) = (r.len(), read.len());
        let to_end = l - start;
        if n <= to_end {
            return starts_with(r, read, start);
        }
        if !starts_with(r, &read[..to_end], start) {
            return false;
        }
        let mut pos = to_end;
        while n - pos >= l {
            if !read[pos..].starts_with(r) {
                return false;
            }
            pos += l;
        }
        pos == n || r.starts_with(&read[pos..])
    }

    /// _RefMatcher._contains(): `read` occurs in ref i (wrapping around if circular)
    fn contains(&self, i: usize, read: &[u8]) -> bool {
        let r = &self.seqs[i];
        if contains(r, read) {
            return true;
        }
        if !self.circular {
            return false;
        }
        let (l, n) = (r.len(), read.len());
        if n <= l {
            let mut junction = r[l - (n - 1)..].to_vec();
            junction.extend_from_slice(&r[..n - 1]);
            return contains(&junction, read);
        }
        contains(&r.repeat(n / l + 2), read)
    }

    /// EditMatcher._window(): `length` bases of ref i beginning at `start` (wrapping if circular)
    fn window(&self, i: usize, start: i64, length: usize) -> Vec<u8> {
        let r = &self.seqs[i];
        let l = r.len();
        if !self.circular {
            let s = (start as usize).min(l);
            let e = (s + length).min(l);
            return r[s..e].to_vec();
        }
        let s = pymod(start, l as i64) as usize;
        let e = s + length;
        if e <= l {
            return r[s..e].to_vec();
        }
        if e <= 2 * l {
            let mut w = r[s..].to_vec();
            w.extend_from_slice(&r[..e - l]);
            return w;
        }
        r.repeat(e / l + 1)[s..e].to_vec()
    }
}

/// get_candidates()'s sample positions: round(j * last / (n_samples - 1)) for each j,
/// with Python's float division and round-half-to-even
fn sample_positions(last: usize, n_samples: usize) -> Vec<usize> {
    if n_samples > 1 {
        let mut p: Vec<usize> = (0..n_samples)
            .map(|j| ((j * last) as f64 / (n_samples - 1) as f64).round_ties_even() as usize)
            .collect();
        p.sort_unstable();
        p.dedup();
        p
    } else {
        vec![0]
    }
}

pub struct Exact {
    refs: Refs,
    k: usize,
    n_samples: usize,
    step: usize,
    pool: SetPool,
    index: KmerMap<u32>,                // canonical k-mer -> set id (see build_kmer_index)
    pos_index: PosIndex, // forward k-mer -> (ref, position) (build_position_index)
}

impl Exact {
    pub fn new(seqs: Vec<Vec<u8>>, k: usize, n_samples: usize, circular: bool, step: usize) -> Self {
        let refs = Refs { seqs, circular };

        // build_kmer_index(): canonical k-mer -> the set of sequences holding it
        let (packed, other) = gather_kmers(refs.seqs.len(), k, true, |i| {
            let seq = &refs.seqs[i];
            let l = seq.len();
            let (circ, n_kmers): (Vec<u8>, usize) = if circular {
                let rep = seq.repeat(k / l + 2);
                (rep[..(l + k - 1).min(rep.len())].to_vec(), l)
            } else {
                (seq.clone(), (l + 1).saturating_sub(k))
            };
            Some((circ, (0..n_kmers).collect()))
        });
        let mut pool = SetPool::new();
        let mut keys = Vec::new();
        let mut buf = Vec::new();
        for (km, run) in key_runs(&packed) {
            keys.push((km, pool.intern_sorted(run_seqs(run, &mut buf))));
        }
        drop(packed);
        let mut sets: FxHashMap<Box<[u8]>, Vec<u32>> = FxHashMap::default();
        for (km, i, _) in other {
            let v = sets.entry(km).or_default();
            if v.last() != Some(&i) {
                v.push(i);
            }
        }
        let bytes = sets.into_iter().map(|(km, s)| (km, pool.intern(s))).collect();
        let index = KmerMap::build(keys, bytes);

        // build_position_index(): forward k-mer -> (sequence, position), every `step` positions
        let (packed, other) = gather_kmers(refs.seqs.len(), k, false, |i| {
            let seq = &refs.seqs[i];
            let l = seq.len();
            if l < k {
                return None;
            }
            Some(if circular {
                let mut c = seq.clone();
                c.extend_from_slice(&seq[..k - 1]);
                (c, (0..l).step_by(step).collect())
            } else {
                let last = l - k;
                let mut p: Vec<usize> = (0..=last).step_by(step).collect();
                if *p.last().unwrap() != last {
                    p.push(last);
                }
                (seq.clone(), p)
            })
        });
        let pos_index = PosIndex::build(&packed, other);

        Exact {
            refs,
            k,
            n_samples,
            step,
            pool,
            index,
            pos_index,
        }
    }

    /// get_candidates(): refs holding every sampled canonical k-mer (None: none do)
    fn candidates(&self, s: &[u8]) -> Option<Vec<u32>> {
        let k = self.k;
        let last = s.len() - k;
        let mut cands: Option<Vec<u32>> = None;
        for p in sample_positions(last, self.n_samples) {
            let hit = match self.index.get_key(canon_key(&s[p..p + k])) {
                Some(&id) => &self.pool.sets[id as usize],
                None => return None,
            };
            cands = Some(match cands {
                None => hit.clone(),
                Some(c) => c.into_iter().filter(|x| hit.binary_search(x).is_ok()).collect(),
            });
            if cands.as_ref().unwrap().is_empty() {
                return None;
            }
        }
        cands
    }

    /// ExactMatcher.hits(): sorted indices of refs `seq` (uppercased) exactly matches
    pub fn hits(&self, s: &[u8], min_frac: f64) -> Vec<u32> {
        let l = s.len();
        if l == 0 {
            return Vec::new();
        }
        let k = self.k;
        let n_refs = self.refs.seqs.len() as u32;
        let cands: Vec<u32> = if l >= k {
            match self.candidates(s) {
                Some(c) => c,
                None => return Vec::new(),
            }
        } else {
            (0..n_refs).collect()
        };
        let cands: Vec<u32> = cands
            .into_iter()
            .filter(|&i| {
                let rl = self.refs.seqs[i as usize].len();
                (l as f64) >= min_frac * (rl as f64) && (self.refs.circular || l <= rl)
            })
            .collect();
        if cands.is_empty() {
            return Vec::new();
        }
        let s_rc = revcomp(s);

        // short reads: plain substring checks
        if (l + 1) < k + self.step {
            // (L - k + 1 < step)
            return cands
                .into_iter()
                .filter(|&i| self.refs.contains(i as usize, s) || self.refs.contains(i as usize, &s_rc))
                .collect();
        }

        let is_cand: FxHashSet<u32> = cands.iter().copied().collect();
        let mut hits: FxHashSet<u32> = FxHashSet::default();
        for strand in [s, s_rc.as_slice()] {
            let mut tried: FxHashSet<(u32, i64)> = FxHashSet::default();
            for t in 0..self.step {
                let entries = match self.pos_index.get(&strand[t..t + k]) {
                    Some(e) if !e.is_empty() => e,
                    _ => continue,
                };
                for &(i, pos) in entries {
                    if hits.contains(&i) || !is_cand.contains(&i) {
                        continue;
                    }
                    let mut start = pos as i64 - t as i64;
                    if self.refs.circular {
                        start = pymod(start, self.refs.seqs[i as usize].len() as i64);
                    } else if start < 0 {
                        continue;
                    }
                    if !tried.insert((i, start)) {
                        continue;
                    }
                    if self.refs.matches_at(i as usize, strand, start) {
                        hits.insert(i);
                    }
                }
            }
            if hits.len() == cands.len() {
                break;
            }
        }
        let mut out: Vec<u32> = hits.into_iter().collect();
        out.sort_unstable();
        out
    }
}

pub struct Edit {
    refs: Refs,
    k: usize,
    max_edits: usize,
    pool: SetPool,
    dense: PosIndex,                 // forward k-mer -> (ref, position), every position
    kmer_refs: KmerMap<u32>,         // forward k-mer -> set id of refs holding it
}

impl Edit {
    pub fn new(seqs: Vec<Vec<u8>>, k: usize, circular: bool, max_edits: usize) -> Self {
        let refs = Refs { seqs, circular };

        // build_dense_position_index(): forward k-mer -> (sequence, position), and the set of
        // sequences holding it
        let (packed, other) = gather_kmers(refs.seqs.len(), k, false, |i| {
            let seq = &refs.seqs[i];
            let l = seq.len();
            if l < k {
                return None;
            }
            Some(if circular {
                let mut c = seq.clone();
                c.extend_from_slice(&seq[..k - 1]);
                (c, (0..l).collect())
            } else {
                (seq.clone(), (0..l - k + 1).collect())
            })
        });
        let mut pool = SetPool::new();
        let mut keys = Vec::new();
        let mut buf = Vec::new();
        for (km, run) in key_runs(&packed) {
            keys.push((km, pool.intern_sorted(run_seqs(run, &mut buf))));
        }
        let mut sets: FxHashMap<Box<[u8]>, Vec<u32>> = FxHashMap::default();
        for (km, i, _) in &other {
            let v = sets.entry(km.clone()).or_default();
            if v.last() != Some(i) {
                v.push(*i);
            }
        }
        let bytes = sets.into_iter().map(|(km, s)| (km, pool.intern(s))).collect();
        let kmer_refs = KmerMap::build(keys, bytes);
        let dense = PosIndex::build(&packed, other);

        Edit {
            refs,
            k,
            max_edits,
            pool,
            dense,
            kmer_refs,
        }
    }

    fn eligible(&self, i: usize, n: usize, min_frac: f64) -> bool {
        let l = self.refs.seqs[i].len();
        (n as f64) >= min_frac * (l as f64) && (self.refs.circular || n <= l + self.max_edits)
    }

    /// EditMatcher.distances(): (ref index, edit distance) for refs within max_edits
    pub fn distances(&self, s: &[u8], min_frac: f64) -> Vec<(u32, usize)> {
        let n = s.len();
        if n == 0 {
            return Vec::new();
        }
        let (k, e) = (self.k, self.max_edits);
        let n_slots = n / k;
        if n_slots < e + 1 {
            return self.full_alignment(s, min_frac);
        }
        let needed = n_slots - e;

        let mut best: FxHashMap<u32, usize> = FxHashMap::default();
        let s_rc = revcomp(s);
        for strand in [s, s_rc.as_slice()] {
            let kmers: Vec<&[u8]> = (0..n_slots).map(|j| &strand[j * k..(j + 1) * k]).collect();
            let mut counts: FxHashMap<u32, usize> = FxHashMap::default();
            for km in &kmers {
                if let Some(&id) = self.kmer_refs.get(km) {
                    for &i in &self.pool.sets[id as usize] {
                        *counts.entry(i).or_insert(0) += 1;
                    }
                }
            }
            let mut passing: Vec<u32> = counts
                .iter()
                .filter(|&(&i, &c)| c >= needed && self.eligible(i as usize, n, min_frac))
                .map(|(&i, _)| i)
                .collect();
            if passing.is_empty() {
                continue;
            }
            passing.sort_unstable();

            let mut diagonals: FxHashMap<u32, Vec<i64>> = FxHashMap::default();
            for (j, km) in kmers.iter().enumerate().take(e + 1) {
                if let Some(entries) = self.dense.get(km) {
                    for &(i, p) in entries {
                        let mut diag = p as i64 - (j * k) as i64;
                        if self.refs.circular {
                            diag = pymod(diag, self.refs.seqs[i as usize].len() as i64);
                        }
                        diagonals.entry(i).or_default().push(diag);
                    }
                }
            }
            for v in diagonals.values_mut() {
                v.sort_unstable();
                v.dedup();
            }

            let none = Vec::new();
            for i in passing {
                let d = self.best_distance(i as usize, strand, diagonals.get(&i).unwrap_or(&none));
                if let Some(d) = d {
                    let cur = best.entry(i).or_insert(usize::MAX);
                    if d < *cur {
                        *cur = d;
                    }
                }
            }
        }
        best.into_iter().collect()
    }

    /// EditMatcher._best_distance()
    fn best_distance(&self, i: usize, read: &[u8], diagonals: &[i64]) -> Option<usize> {
        let l = self.refs.seqs[i].len() as i64;
        let (n, e) = (read.len(), self.max_edits as i64);

        for &diag in diagonals {
            if (self.refs.circular || 0 <= diag) && self.refs.matches_at(i, read, diag) {
                return Some(0);
            }
        }

        let mut starts: Vec<i64> = Vec::new();
        for &diag in diagonals {
            for st in (diag - e)..=(diag + e) {
                if self.refs.circular {
                    starts.push(pymod(st, l));
                } else if 0 <= st && st < l {
                    starts.push(st);
                }
            }
        }
        starts.sort_unstable();
        starts.dedup();

        let (mut best, mut limit): (Option<usize>, usize) = (None, self.max_edits);
        for st in starts {
            let window = self.refs.window(i, st, n + limit);
            if let Some(d) = edlib::prefix_distance(read, &window, limit)
                && best.is_none_or(|b| d < b)
            {
                best = Some(d);
                if d <= 1 {
                    // 0 was already ruled out above
                    break;
                }
                limit = d - 1;
            }
        }
        best
    }

    /// EditMatcher._distances_full_alignment()
    fn full_alignment(&self, s: &[u8], min_frac: f64) -> Vec<(u32, usize)> {
        let (n, e) = (s.len(), self.max_edits);
        let s_rc = revcomp(s);
        let mut best = Vec::new();
        for i in 0..self.refs.seqs.len() {
            if !self.eligible(i, n, min_frac) {
                continue;
            }
            let l = self.refs.seqs[i].len();
            let target = if self.refs.circular {
                self.refs.window(i, 0, l + n + e - 1)
            } else {
                self.refs.seqs[i].clone()
            };
            let d = [s, s_rc.as_slice()]
                .iter()
                .filter_map(|r| edlib::distance(r, &target, e))
                .min();
            if let Some(d) = d {
                best.push((i as u32, d));
            }
        }
        best
    }
}

pub enum Matcher {
    Exact(Exact),
    Edit(Edit),
}

impl Matcher {
    /// {seq index: edit distance} as (index, distance) pairs (`seq` is uppercased here)
    pub fn distances(&self, seq: &[u8], min_frac: f64) -> Vec<(u32, usize)> {
        let s = seq.to_ascii_uppercase();
        match self {
            Matcher::Exact(m) => m.hits(&s, min_frac).into_iter().map(|i| (i, 0)).collect(),
            Matcher::Edit(m) => m.distances(&s, min_frac),
        }
    }
}

/// collapse_to_refs(): ref -> (lowest distance among its seqs, its seqs at that distance)
fn collapse(seq_dists: &[(u32, usize)], seq_to_ref: &[u32]) -> FxHashMap<u32, (usize, Vec<u32>)> {
    let mut out: FxHashMap<u32, (usize, Vec<u32>)> = FxHashMap::default();
    for &(si, d) in seq_dists {
        let ri = seq_to_ref[si as usize];
        match out.get_mut(&ri) {
            None => {
                out.insert(ri, (d, vec![si]));
            }
            Some(cur) => {
                if d < cur.0 {
                    *cur = (d, vec![si]);
                } else if d == cur.0 {
                    cur.1.push(si);
                }
            }
        }
    }
    out
}

/// summarize_distances(): (best refs, sorted; best distance; next-best distance)
fn summarize(dists: &[(u32, usize)]) -> Option<(Vec<u32>, usize, Option<usize>)> {
    let best = dists.iter().map(|&(_, d)| d).min()?;
    let mut hits: Vec<u32> = dists.iter().filter(|&&(_, d)| d == best).map(|&(i, _)| i).collect();
    hits.sort_unstable();
    let next = dists.iter().map(|&(_, d)| d).filter(|&d| d != best).min();
    Some((hits, best, next))
}

/// one read's (or pair's) result, as process_batch() builds it (minus the read id, lengths
/// and records, which Python adds back)
pub struct Assigned {
    pub hits: Vec<u32>,
    pub hit_seqs: Vec<u32>,
    pub best: usize,
    pub next_best: Option<usize>,
}

pub fn assign_single(m: &Matcher, seq_to_ref: &[u32], seq: &[u8], min_frac: f64) -> Option<Assigned> {
    let collapsed = collapse(&m.distances(seq, min_frac), seq_to_ref);
    let per_ref: Vec<(u32, usize)> = collapsed.iter().map(|(&ri, v)| (ri, v.0)).collect();
    let (hits, best, next_best) = summarize(&per_ref)?;
    let mut hit_seqs = Vec::new();
    for ri in &hits {
        let mut s = collapsed[ri].1.clone();
        s.sort_unstable();
        hit_seqs.extend(s);
    }
    Some(Assigned {
        hits,
        hit_seqs,
        best,
        next_best,
    })
}

pub fn assign_pair(
    m: &Matcher,
    seq_to_ref: &[u32],
    r1: &[u8],
    r2: &[u8],
    min_frac: f64,
) -> Option<Assigned> {
    let c1 = collapse(&m.distances(r1, min_frac), seq_to_ref);
    if c1.is_empty() {
        return None;
    }
    let c2 = collapse(&m.distances(r2, min_frac), seq_to_ref);
    if c2.is_empty() {
        return None;
    }
    let shared: Vec<(u32, usize)> = c1
        .iter()
        .filter_map(|(ri, a)| c2.get(ri).map(|b| (*ri, a.0 + b.0)))
        .collect();
    let (hits, best, next_best) = summarize(&shared)?;
    let mut hit_seqs = Vec::new();
    for ri in &hits {
        let mut s: Vec<u32> = c1[ri].1.iter().chain(c2[ri].1.iter()).copied().collect();
        s.sort_unstable();
        s.dedup();
        hit_seqs.extend(s);
    }
    Some(Assigned {
        hits,
        hit_seqs,
        best,
        next_best,
    })
}

/// runs a batch (reads, or pairs as consecutive mates) across threads, keeping input order
pub fn assign_batch(
    m: &Matcher,
    seq_to_ref: &[u32],
    seqs: &[&[u8]],
    paired: bool,
    min_frac: f64,
) -> Vec<Option<Assigned>> {
    if paired {
        seqs.par_chunks(2)
            .map(|p| assign_pair(m, seq_to_ref, p[0], p[1], min_frac))
            .collect()
    } else {
        seqs.par_iter()
            .map(|s| assign_single(m, seq_to_ref, s, min_frac))
            .collect()
    }
}

/// A read-hits table written sorted by key, ties in the order the rows came (as
/// SortedRowWriter in match_reads.py): up to `buffer_rows` rows are held in memory, beyond that
/// sorted chunks go to a temporary directory beside the table and are merged at the end, so
/// memory stays bounded however many rows there are. The directory is removed when this is
/// dropped, finished or not; the table is written only by finish().
pub struct SortedRows {
    path: std::path::PathBuf,
    header: String,
    buffer_rows: usize,
    rows: Vec<((i64, i64), String)>,
    chunks: Vec<std::path::PathBuf>,
    tmp_dir: Option<std::path::PathBuf>,
}

impl SortedRows {
    pub fn new(path: &str, header: String, buffer_rows: usize) -> Self {
        SortedRows {
            path: path.into(),
            header,
            buffer_rows: buffer_rows.max(1),
            rows: Vec::new(),
            chunks: Vec::new(),
            tmp_dir: None,
        }
    }

    pub fn add(&mut self, key: (i64, i64), line: String) -> std::io::Result<()> {
        self.rows.push((key, line));
        if self.rows.len() >= self.buffer_rows {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> std::io::Result<()> {
        use std::io::Write;
        if self.tmp_dir.is_none() {
            let parent = match self.path.parent() {
                Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
                _ => std::path::PathBuf::from("."),
            };
            let mut n = 0;
            let dir = loop {
                let d = parent.join(format!(".sorting-{}-{n}", std::process::id()));
                match std::fs::create_dir(&d) {
                    Ok(()) => break d,
                    Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => n += 1,
                    Err(e) => return Err(e),
                }
            };
            self.tmp_dir = Some(dir);
        }
        self.rows.sort_by_key(|r| r.0); // stable: ties stay in the order added
        let chunk = self.tmp_dir.as_ref().unwrap().join(format!("chunk-{}.tsv", self.chunks.len()));
        let mut f = std::io::BufWriter::new(std::fs::File::create(&chunk)?);
        for ((a, b), line) in self.rows.drain(..) {
            writeln!(f, "{a}\t{b}\t{line}")?;
        }
        f.flush()?;
        self.chunks.push(chunk);
        Ok(())
    }

    pub fn finish(mut self) -> std::io::Result<()> {
        use std::io::{BufRead, Write};
        let mut out = std::io::BufWriter::new(std::fs::File::create(&self.path)?);
        writeln!(out, "{}", self.header)?;
        if self.chunks.is_empty() {
            self.rows.sort_by_key(|r| r.0);
            for (_, line) in &self.rows {
                writeln!(out, "{line}")?;
            }
        } else {
            if !self.rows.is_empty() {
                self.spill()?;
            }
            // a k-way merge; on equal keys the earlier chunk goes first, so ties stay in input order
            let mut readers = Vec::new();
            for c in &self.chunks {
                readers.push(std::io::BufReader::new(std::fs::File::open(c)?).lines());
            }
            let parse = |l: String| -> std::io::Result<((i64, i64), String)> {
                let mut parts = l.splitn(3, '\t');
                let bad = || std::io::Error::new(std::io::ErrorKind::InvalidData, "bad sort chunk");
                let a = parts.next().and_then(|x| x.parse().ok()).ok_or_else(bad)?;
                let b = parts.next().and_then(|x| x.parse().ok()).ok_or_else(bad)?;
                Ok(((a, b), parts.next().ok_or_else(bad)?.to_string()))
            };
            let mut heads: Vec<Option<((i64, i64), String)>> = Vec::new();
            let mut heap = std::collections::BinaryHeap::new();
            for (ci, r) in readers.iter_mut().enumerate() {
                let head = r.next().transpose()?.map(parse).transpose()?;
                if let Some((key, _)) = &head {
                    heap.push(std::cmp::Reverse((*key, ci)));
                }
                heads.push(head);
            }
            while let Some(std::cmp::Reverse((_, ci))) = heap.pop() {
                let (_, line) = heads[ci].take().unwrap();
                writeln!(out, "{line}")?;
                let head = readers[ci].next().transpose()?.map(parse).transpose()?;
                if let Some((key, _)) = &head {
                    heap.push(std::cmp::Reverse((*key, ci)));
                }
                heads[ci] = head;
            }
        }
        out.flush()
    }
}

impl Drop for SortedRows {
    fn drop(&mut self) {
        if let Some(d) = &self.tmp_dir {
            let _ = std::fs::remove_dir_all(d);
        }
    }
}
