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

use crate::readio::{GzOut, Rec, cut_record, gz_member, spawn_reader};
use crate::reads::Progress;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io;
use std::sync::Arc;
use std::sync::mpsc::{SyncSender, sync_channel};

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
#[pyo3(signature = (r1, r2=None, k=13, min_shared=16, tolerance=20))]
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

// Whole files (`bit filter-foldbacks`): reading, finding and cutting folds, and writing the
// gzipped output all happen here, so the work spreads over the threads. Records are read in
// batches on threads of their own and written out gzipped on several threads (readio.rs); each
// batch is cut into chunks scanned on the threads of the pool, and a writer thread appends the
// chunks' output in input order. Reads are written as they were in the input (header, comments
// and qualities), cut ones to the stretch kept.

const SHORT_BATCH: usize = 16_384;
const SHORT_CHUNK: usize = 128;
const LONG_BATCH: usize = 1_024;
const LONG_CHUNK: usize = 8;
// the unpaired middles kept for the summary's quartiles (long reads)
const MAX_MIDDLES: usize = 1_000_000;

/// what a pass found, in the order of bit's Counts.FIELDS
#[derive(Default)]
struct Tally([u64; 12]);

// Tally's fields
const IN: usize = 0;
const FOLDBACK_PAIRS: usize = 1;
const PAIRS_FOLD_IN_READ: usize = 2;
const PAIRS_FOLD_BETWEEN: usize = 3;
const READS_CUT: usize = 4;
const FOLDS_CUT: usize = 5;
const MATES_AS_SINGLES: usize = 6;
const DROPPED_SHORT: usize = 7;
const PAIRS_OUT: usize = 8;
const READS_OUT: usize = 9;
const BASES_IN: usize = 10;
const BASES_OUT: usize = 11;

/// one chunk's output: a gzip member for each output file (empty if nothing went to it), what
/// it found, and the unpaired middles of its cut long reads
struct ChunkOut {
    text: Vec<Vec<u8>>,
    tally: Tally,
    middles: Vec<i32>,
}

/// a record onto `out`, whole or cut to `start..end` of its bases (`scratch`: two buffers)
fn put(out: &mut Vec<u8>, raw: &[u8], len: usize, start: usize, end: usize, fastq: bool, scratch: &mut [Vec<u8>; 2]) {
    if start == 0 && end == len {
        out.extend_from_slice(raw);
    } else {
        let [rec, buf] = scratch;
        rec.clear();
        rec.extend_from_slice(raw);
        cut_record(rec, start, end, fastq, buf);
        out.extend_from_slice(rec);
    }
}

struct Settings {
    k: usize,
    min_shared: usize,
    min_length: usize,
    threshold: u64,
    fastq: bool,
}

/// short reads or pairs: outputs are (reads) single-end, or (r1, r2, fold-back mates) paired
fn short_chunk(chunk: &[Rec], st: &Settings) -> ChunkOut {
    let paired = chunk.first().is_some_and(|r| r.seq2.is_some());
    let mut text: Vec<Vec<u8>> = vec![Vec::new(); if paired { 3 } else { 1 }];
    let mut t = Tally::default();
    let mut scratch = [Vec::new(), Vec::new()];
    for rec in chunk {
        let (a, raw_a) = (&rec.seq1, &rec.raw1);
        let (s1, e1, f1) = trim(a, st.k, st.min_shared, 20);
        t.0[IN] += 1;
        t.0[BASES_IN] += a.len() as u64;
        let (Some(b), Some(raw_b)) = (&rec.seq2, &rec.raw2) else {
            t.0[READS_CUT] += (f1 > 0) as u64;
            t.0[FOLDS_CUT] += f1 as u64;
            if f1 == 0 || e1 - s1 >= st.min_length {
                put(&mut text[0], raw_a, a.len(), s1, e1, st.fastq, &mut scratch);
                t.0[READS_OUT] += 1;
                t.0[BASES_OUT] += (e1 - s1) as u64;
            } else {
                t.0[DROPPED_SHORT] += 1;
            }
            continue;
        };
        let (s2, e2, f2) = trim(b, st.k, st.min_shared, 20);
        let (same, _) = shared(a, b, st.k);
        t.0[BASES_IN] += b.len() as u64;
        t.0[READS_CUT] += (f1 > 0) as u64 + (f2 > 0) as u64;
        t.0[FOLDS_CUT] += f1 as u64 + f2 as u64;
        let in_read = f1 > 0 || f2 > 0;
        let between = same >= st.min_shared;
        if !in_read && !between {
            text[0].extend_from_slice(raw_a);
            text[1].extend_from_slice(raw_b);
            t.0[PAIRS_OUT] += 1;
            t.0[BASES_OUT] += (a.len() + b.len()) as u64;
            continue;
        }
        t.0[FOLDBACK_PAIRS] += 1;
        t.0[PAIRS_FOLD_IN_READ] += in_read as u64;
        t.0[PAIRS_FOLD_BETWEEN] += (between && !in_read) as u64;
        for (raw, len, s, e) in [(raw_a, a.len(), s1, e1), (raw_b, b.len(), s2, e2)] {
            if e - s >= st.min_length {
                put(&mut text[2], raw, len, s, e, st.fastq, &mut scratch);
                t.0[MATES_AS_SINGLES] += 1;
                t.0[BASES_OUT] += (e - s) as u64;
            } else {
                t.0[DROPPED_SHORT] += 1;
            }
        }
    }
    ChunkOut { text: text.iter().map(|x| gz_member(x)).collect(), tally: t, middles: Vec::new() }
}

/// long reads: each cut at its folds as shifia --pta does (see long_trim)
fn long_chunk(chunk: &[Rec], st: &Settings) -> ChunkOut {
    let mut text = Vec::new();
    let mut t = Tally::default();
    let mut middles = Vec::new();
    let mut scratch = [Vec::new(), Vec::new()];
    for rec in chunk {
        let a = &rec.seq1;
        let (start, end, folds, _, middle) = long_trim(a, st.k, st.threshold, 0, 5, 0.25, 50, 4);
        let (start, end) = (start as usize, end as usize);
        t.0[IN] += 1;
        t.0[BASES_IN] += a.len() as u64;
        if folds > 0 {
            t.0[READS_CUT] += 1;
            t.0[FOLDS_CUT] += folds as u64;
            middles.push(middle);
            if end - start < st.min_length {
                t.0[DROPPED_SHORT] += 1;
                continue;
            }
        }
        put(&mut text, &rec.raw1, a.len(), start, end, st.fastq, &mut scratch);
        t.0[READS_OUT] += 1;
        t.0[BASES_OUT] += (end - start) as u64;
    }
    ChunkOut { text: vec![gz_member(&text)], tally: t, middles }
}

/// appends each batch's chunks to the outputs, in order, on a thread of its own. The first two
/// outputs always exist (an empty one as a valid empty gzip file); a third (paired: the
/// fold-back mates) is made only once something goes to it. Returns whether each got any reads.
#[allow(clippy::type_complexity)]
fn spawn_writer(
    paths: Vec<String>,
) -> io::Result<(SyncSender<Vec<ChunkOut>>, std::thread::JoinHandle<io::Result<Vec<bool>>>)> {
    let mut outs: Vec<Option<GzOut>> = Vec::new();
    for (i, p) in paths.iter().enumerate() {
        outs.push(if i < 2 { Some(GzOut::create(p)?) } else { None });
    }
    let (tx, rx) = sync_channel::<Vec<ChunkOut>>(2);
    let handle = std::thread::spawn(move || {
        let mut wrote = vec![false; paths.len()];
        for chunks in rx {
            for c in chunks {
                for (i, text) in c.text.iter().enumerate() {
                    if text.is_empty() {
                        continue;
                    }
                    if outs[i].is_none() {
                        outs[i] = Some(GzOut::create(&paths[i])?);
                    }
                    outs[i].as_mut().unwrap().push_member(text)?;
                    wrote[i] = true;
                }
            }
        }
        for o in outs.into_iter().flatten() {
            o.finish()?;
        }
        Ok(wrote)
    });
    Ok((tx, handle))
}

fn run_filter(
    read_1: String,
    read_2: Option<String>,
    outputs: Vec<String>,
    mut st: Settings,
    hifi: bool,
    progress: Option<Arc<Progress>>,
) -> io::Result<(Tally, Vec<bool>, Vec<i32>)> {
    let (batch, chunk) = if hifi { (LONG_BATCH, LONG_CHUNK) } else { (SHORT_BATCH, SHORT_CHUNK) };
    let (batches, fastq, spent) = spawn_reader(read_1, read_2, batch, progress)?;
    st.fastq = fastq;
    let (to_writer, writer) = spawn_writer(outputs)?;
    let mut tally = Tally::default();
    let mut middles = Vec::new();
    let mut failed = None;
    for b in batches {
        let b = match b {
            Ok(b) => b,
            Err(e) => {
                failed = Some(e);
                break;
            }
        };
        let outs: Vec<ChunkOut> = b
            .par_chunks(chunk)
            .map(|c| if hifi { long_chunk(c, &st) } else { short_chunk(c, &st) })
            .collect();
        let _ = spent.send(b);
        for o in &outs {
            for (sum, x) in tally.0.iter_mut().zip(o.tally.0) {
                *sum += x;
            }
            let room = MAX_MIDDLES.saturating_sub(middles.len());
            middles.extend(o.middles.iter().take(room));
        }
        if to_writer.send(outs).is_err() {
            break; // the writer failed: its error is below
        }
    }
    drop(to_writer);
    let wrote = writer.join().expect("writer thread")?;
    if let Some(e) = failed {
        return Err(e);
    }
    Ok((tally, wrote, middles))
}

/// Filters a read file, or a pair of them, to the outputs given: single-end or long reads
/// [reads]; paired [read 1, read 2, fold-back mates] (gzipped). Returns what was found (bit's
/// Counts.FIELDS, in order), whether each output got any reads, and the unpaired middles of the
/// cut long reads (the first million). Short reads: k-mers of `k`, `min_shared` to call a fold.
/// Long reads (`hifi`): anchors of `k` sampled at about 1 in `density`. Cut pieces under
/// `min_length` are dropped.
#[pyfunction]
#[pyo3(signature = (read_1, read_2, outputs, k, min_shared, min_length, hifi=false, density=4,
                    progress=None))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn filter_foldback_files(
    py: Python<'_>,
    read_1: String,
    read_2: Option<String>,
    outputs: Vec<String>,
    k: usize,
    min_shared: usize,
    min_length: usize,
    hifi: bool,
    density: u64,
    progress: Option<PyRef<crate::Progress>>,
) -> PyResult<(Vec<u64>, Vec<bool>, Vec<i32>)> {
    if !(5..=31).contains(&k) {
        return Err(pyo3::exceptions::PyValueError::new_err("k must be between 5 and 31"));
    }
    let want = if read_2.is_some() { 3 } else { 1 };
    if outputs.len() != want || (hifi && read_2.is_some()) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "outputs: [reads] for single-end or long reads, [read 1, read 2, fold-back mates] for pairs",
        ));
    }
    let st = Settings { k, min_shared, min_length, threshold: crate::anchors::threshold(density), fastq: false };
    let progress = crate::arc(&progress);
    let (tally, wrote, middles) = py
        .detach(|| run_filter(read_1, read_2, outputs, st, hifi, progress))
        .map_err(crate::to_py_err)?;
    Ok((tally.0.to_vec(), wrote, middles))
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
