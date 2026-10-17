//! Placing found sequences in reads and counting them (the rules are set out in the
//! docstrings of Matcher and count_in_reads in bit/modules/find_overrep.py).

use crate::edlib;
use crate::kmers::{BATCH, Prog, for_each_kmer, revcomp};
use crate::reads::SeqReader;
use crate::trim::{Cut, Trimmer};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io;

/// a found sequence's best placement in a read
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub seq: usize,
    pub start: Option<usize>, // None: the read is within the sequence
    pub is_rc: bool,
    pub edits: usize,
}

struct Pattern {
    seq: usize,
    is_rc: bool,
    bases: Vec<u8>,
    edits: usize,
}

pub struct Matcher {
    patterns: Vec<Pattern>,
    lengths: Vec<usize>,
    max_diff: f64,
    min_read_len: usize,
    q: usize,
    lut: Vec<Vec<u64>>, // per 64 patterns, per q-mer code: which patterns hold it
}

pub fn allowed_edits(len: usize, max_diff: f64) -> usize {
    (max_diff * len as f64).floor() as usize
}

impl Matcher {
    pub fn new(seqs: &[Vec<u8>], max_diff: f64, min_read_len: usize) -> Self {
        let mut patterns = Vec::with_capacity(2 * seqs.len());
        for (i, s) in seqs.iter().enumerate() {
            let e = allowed_edits(s.len(), max_diff);
            patterns.push(Pattern {
                seq: i,
                is_rc: false,
                bases: s.clone(),
                edits: e,
            });
            patterns.push(Pattern {
                seq: i,
                is_rc: true,
                bases: revcomp(s),
                edits: e,
            });
        }
        // q such that a sequence within its allowed edits must share an exact q-mer with
        // the read (pigeonhole), so the prefilter can't cause a miss
        let mut q = seqs
            .iter()
            .map(|s| s.len() / (allowed_edits(s.len(), max_diff) + 1))
            .min()
            .unwrap_or(11);
        q = q.min(min_read_len / (allowed_edits(min_read_len, max_diff) + 1));
        let q = q.clamp(4, 11);
        let words = patterns.len().div_ceil(64).max(1);
        let mut lut = vec![vec![0u64; 1 << (2 * q)]; words];
        for (p, pat) in patterns.iter().enumerate() {
            for_each_kmer(&pat.bases, q, |_, fwd, _| {
                lut[p / 64][fwd as usize] |= 1u64 << (p % 64)
            });
        }
        Matcher {
            patterns,
            lengths: seqs.iter().map(|s| s.len()).collect(),
            max_diff,
            min_read_len,
            q,
            lut,
        }
    }

    /// patterns sharing a q-mer with the read, in order
    fn candidates(&self, read: &[u8]) -> Vec<usize> {
        let mut masks = vec![0u64; self.lut.len()];
        for_each_kmer(read, self.q, |_, fwd, _| {
            for (w, table) in self.lut.iter().enumerate() {
                masks[w] |= table[fwd as usize];
            }
        });
        let mut out = Vec::new();
        for (w, &m) in masks.iter().enumerate() {
            let mut m = m;
            while m != 0 {
                out.push(w * 64 + m.trailing_zeros() as usize);
                m &= m - 1;
            }
        }
        out
    }

    /// the best placement per found sequence (see Matcher.hits in find_overrep.py), from
    /// the prefilter's candidates, or from every pattern if `all`
    pub fn hits(&self, read: &[u8], all: bool) -> Vec<Hit> {
        let cands: Vec<usize> = if all {
            (0..self.patterns.len()).collect()
        } else {
            self.candidates(read)
        };
        let mut found: Vec<Hit> = Vec::new();
        for p in cands {
            let pat = &self.patterns[p];
            let hit = if read.len() + pat.edits >= pat.bases.len() {
                match edlib::locate(&pat.bases, read, pat.edits) {
                    Some((d, start)) => Hit {
                        seq: pat.seq,
                        start: Some(start),
                        is_rc: pat.is_rc,
                        edits: d,
                    },
                    None => continue,
                }
            } else if read.len() >= self.min_read_len {
                let er = allowed_edits(read.len(), self.max_diff);
                match edlib::distance(read, &pat.bases, er) {
                    Some(d) => Hit {
                        seq: pat.seq,
                        start: None,
                        is_rc: pat.is_rc,
                        edits: d,
                    },
                    None => continue,
                }
            } else {
                continue;
            };
            match found.iter_mut().find(|h| h.seq == hit.seq) {
                Some(h) => {
                    if hit.edits < h.edits {
                        *h = hit;
                    }
                }
                None => found.push(hit),
            }
        }
        // a read holding one found sequence in full isn't also counted as a piece of another
        if found.iter().any(|h| h.start.is_some()) {
            found.retain(|h| h.start.is_some());
        } else if found.len() > 1 {
            // and one only within found sequences counts for the shortest
            let shortest = found.iter().map(|h| self.lengths[h.seq]).min().unwrap();
            found.retain(|h| self.lengths[h.seq] == shortest);
        }
        found
    }
}

#[derive(Default)]
pub struct Tally {
    pub reads: u64,
    pub reads_with: Vec<u64>,
    pub full: Vec<u64>,
    pub rc: Vec<u64>,
    pub edits: Vec<FxHashMap<usize, u64>>,
    pub starts: Vec<FxHashMap<usize, u64>>,
    pub adapters: FxHashMap<String, u64>,
    pub adapter_only: u64,
    pub explained: u64,
}

impl Tally {
    fn new(n: usize) -> Self {
        Tally {
            reads_with: vec![0; n],
            full: vec![0; n],
            rc: vec![0; n],
            edits: vec![FxHashMap::default(); n],
            starts: vec![FxHashMap::default(); n],
            ..Default::default()
        }
    }
}

#[derive(Default)]
pub struct PairTally {
    pub pairs: u64,
    pub pairs_with: Vec<u64>, // either read holds it
    pub pairs_both: Vec<u64>, // both reads hold it
    pub explained: u64,
    pub explained_both: u64,
    pub adapter_only: u64,
}

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

struct CountState {
    tallies: Vec<Tally>,
    pairs: PairTally,
    hits: FxHashMap<Vec<u8>, Vec<Hit>>, // trimmed read -> its hits (only a cache)
    names: FxHashMap<Cut, String>,
}

fn read_batches(readers: &mut [SeqReader]) -> io::Result<Vec<Vec<Vec<u8>>>> {
    readers.iter_mut().map(|r| r.batch(BATCH)).collect()
}

/// tallies one batch of reads per file (read pairs in step, for two)
fn tally_batches(
    batches: &[Vec<Vec<u8>>],
    st: &mut CountState,
    matcher: &Matcher,
    trimmer: &Trimmer,
    min_insert_len: usize,
) {
    let mut per_mate: Vec<Vec<(Vec<usize>, bool)>> = Vec::new();
    for (fi, batch) in batches.iter().enumerate() {
        // trimming, then placing each distinct trimmed read not seen before, in parallel
        let looked: Vec<(Vec<u8>, Option<Cut>)> = batch
            .par_iter()
            .map(|read| {
                let (keep, cut) = trimmer.trim(read);
                (read[..keep].to_vec(), cut)
            })
            .collect();
        let mut todo: Vec<&Vec<u8>> = looked
            .iter()
            .map(|(t, _)| t)
            .filter(|t| !st.hits.contains_key(*t))
            .collect();
        todo.sort_unstable();
        todo.dedup();
        let placed: Vec<(Vec<u8>, Vec<Hit>)> = todo
            .par_iter()
            .map(|t| ((*t).clone(), matcher.hits(t, false)))
            .collect();
        st.hits.extend(placed);

        let t = &mut st.tallies[fi];
        let mut results = Vec::with_capacity(batch.len());
        for (trimmed, cut) in &looked {
            let found = &st.hits[trimmed];
            t.reads += 1;
            if let Some(c) = cut {
                let name = st
                    .names
                    .entry(c.clone())
                    .or_insert_with(|| trimmer.cut_name(c))
                    .clone();
                *t.adapters.entry(name).or_insert(0) += 1;
            }
            let short = cut.is_some() && trimmed.len() < min_insert_len;
            if short {
                t.adapter_only += 1;
            }
            for h in found {
                t.reads_with[h.seq] += 1;
                *t.edits[h.seq].entry(h.edits).or_insert(0) += 1;
                if let Some(s) = h.start {
                    t.full[h.seq] += 1;
                    *t.starts[h.seq].entry(s).or_insert(0) += 1;
                }
                if h.is_rc {
                    t.rc[h.seq] += 1;
                }
            }
            if !found.is_empty() {
                t.explained += 1;
            }
            results.push((found.iter().map(|h| h.seq).collect(), short));
        }
        per_mate.push(results);
    }
    // bounding memory on diverse libraries (it's only a cache)
    if st.hits.len() > 2_000_000 {
        st.hits.clear();
    }
    if per_mate.len() == 2 {
        let pairs = &mut st.pairs;
        for ((a, sa), (b, sb)) in per_mate[0].iter().zip(&per_mate[1]) {
            pairs.pairs += 1;
            let mut either: Vec<usize> = a.iter().chain(b).copied().collect();
            either.sort_unstable();
            either.dedup();
            for i in &either {
                pairs.pairs_with[*i] += 1;
                if a.contains(i) && b.contains(i) {
                    pairs.pairs_both[*i] += 1;
                }
            }
            if !either.is_empty() {
                pairs.explained += 1;
            }
            if !a.is_empty() && !b.is_empty() {
                pairs.explained_both += 1;
            }
            if *sa && *sb {
                pairs.adapter_only += 1;
            }
        }
    }
}

/// counts the found sequences in reads (one file, or two for pairs); the next batch is read
/// while the current one is tallied
pub fn count_in_reads(
    paths: &[String],
    matcher: &Matcher,
    trimmer: &Trimmer,
    min_insert_len: usize,
    progress: Prog,
) -> io::Result<(Vec<Tally>, Option<PairTally>)> {
    let n = matcher.lengths.len();
    let mut st = CountState {
        tallies: paths.iter().map(|_| Tally::new(n)).collect(),
        pairs: PairTally {
            pairs_with: vec![0; n],
            pairs_both: vec![0; n],
            ..Default::default()
        },
        hits: FxHashMap::default(),
        names: FxHashMap::default(),
    };
    let mut readers = paths
        .iter()
        .map(|p| SeqReader::open(p, progress.cloned()))
        .collect::<io::Result<Vec<_>>>()?;
    let mut cur = read_batches(&mut readers)?;
    while !cur.iter().all(|b| b.is_empty()) {
        if cur.len() == 2 && cur[0].len() != cur[1].len() {
            return Err(bad("The read files don't have the same number of reads."));
        }
        let (next, ()) = rayon::join(
            || read_batches(&mut readers),
            || tally_batches(&cur, &mut st, matcher, trimmer, min_insert_len),
        );
        cur = next?;
    }
    Ok((st.tallies, (paths.len() == 2).then_some(st.pairs)))
}
