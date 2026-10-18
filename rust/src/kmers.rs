//! K-mer counting passes over reads (see bit/modules/find_overrep.py): counts of every
//! k-mer, per read file, with the 2- and 3-mers of the background model and the valid k-mer
//! windows at each read position; and reads holding chosen seed k-mers.

use crate::reads::{Progress, SeqReader};
use crate::trim::Trimmer;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io;
use std::sync::{Arc, Mutex};

pub type Prog<'a> = Option<&'a Arc<Progress>>;

pub const BATCH: usize = 200_000;

/// reads holding a seed k-mer, each as (read, seed position), per seed
pub type Contexts = Vec<Vec<(Vec<u8>, usize)>>;

static CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

/// calls f(position, forward code, reverse-complement code) for each k-mer of `seq` holding
/// only A/C/G/T, in order
#[inline]
pub fn for_each_kmer<F: FnMut(usize, u64, u64)>(seq: &[u8], k: usize, mut f: F) {
    if k == 0 || seq.len() < k {
        return;
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k - 1);
    let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = CODE[b as usize];
        if c > 3 {
            run = 0;
            continue;
        }
        fwd = ((fwd << 2) | c as u64) & mask;
        rev = (rev >> 2) | (((3 - c) as u64) << shift);
        run += 1;
        if run >= k {
            f(i + 1 - k, fwd, rev);
        }
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// reads per chunk handed to a thread within a batch
const CHUNK: usize = 4096;

fn trim_all(batch: &mut [Vec<u8>], trimmer: Option<&Trimmer>) {
    if let Some(t) = trimmer {
        batch.par_iter_mut().for_each(|r| {
            let (keep, _) = t.trim(r);
            r.truncate(keep);
        });
    }
}

/// calls f(batch) for batches of reads from one file in order, cut where adapter
/// read-through starts if a trimmer is given; the next batch is read (and trimmed) while f
/// works on the current one. Stops early if f returns false.
fn file_batches<F: FnMut(&[Vec<u8>]) -> bool + Send>(
    path: &str,
    trimmer: Option<&Trimmer>,
    progress: Prog,
    max_reads: u64,
    mut f: F,
) -> io::Result<()> {
    let mut reader = SeqReader::open(path, progress.cloned())?;
    // reads taken so far; at most max_reads (0: no limit) are
    let mut taken: u64 = 0;
    let mut next_batch = |taken: &mut u64| -> io::Result<Vec<Vec<u8>>> {
        let want = if max_reads == 0 {
            BATCH
        } else {
            (max_reads - *taken).min(BATCH as u64) as usize
        };
        if want == 0 {
            return Ok(Vec::new());
        }
        let mut b = reader.batch(want)?;
        *taken += b.len() as u64;
        trim_all(&mut b, trimmer);
        Ok(b)
    };
    let mut cur = next_batch(&mut taken)?;
    while !cur.is_empty() {
        let (next, go) = rayon::join(|| next_batch(&mut taken), || f(&cur));
        if !go {
            break;
        }
        cur = next?;
    }
    Ok(())
}

/// as file_batches, for each file at the same time: f(file index, batch)
fn all_batches<F: Fn(usize, &[Vec<u8>]) + Sync>(
    paths: &[String],
    trimmer: Option<&Trimmer>,
    progress: Prog,
    max_reads: u64,
    f: F,
) -> io::Result<()> {
    paths.par_iter().enumerate().try_for_each(|(fi, path)| {
        file_batches(path, trimmer, progress, max_reads, |batch| {
            f(fi, batch);
            true
        })
    })
}

/// one accumulator per thread (each only ever locked by its own thread), so counting needs
/// no sharing until the end
fn per_thread<T: Default>() -> Vec<Mutex<T>> {
    (0..rayon::current_num_threads() + 1)
        .map(|_| Mutex::new(T::default()))
        .collect()
}

fn mine<T>(slots: &[Mutex<T>]) -> std::sync::MutexGuard<'_, T> {
    let i = rayon::current_thread_index().unwrap_or(slots.len() - 1);
    slots[i].lock().unwrap()
}

pub struct Counts {
    pub by_file: Vec<Vec<u64>>,
    pub tri: Vec<u64>,
    pub di: Vec<u64>,
    pub windows_at: Vec<u64>,
    pub n_reads: u64,
}

#[derive(Default)]
struct CountAcc {
    by_file: Vec<Option<Vec<u64>>>,
    tri: Vec<u64>,
    di: Vec<u64>,
    windows_at: Vec<u64>,
    n_reads: u64,
}

pub fn count(
    paths: &[String],
    k: usize,
    trimmer: Option<&Trimmer>,
    progress: Prog,
    max_reads: u64,
) -> io::Result<Counts> {
    let size = 1usize << (2 * k);
    let slots: Vec<Mutex<CountAcc>> = per_thread();
    all_batches(paths, trimmer, progress, max_reads, |fi, batch| {
        batch.par_chunks(CHUNK).for_each(|chunk| {
            let mut guard = mine(&slots);
            let a = &mut *guard;
            if a.by_file.is_empty() {
                a.by_file = vec![None; paths.len()];
                a.tri = vec![0; 64];
                a.di = vec![0; 16];
            }
            let counts = a.by_file[fi].get_or_insert_with(|| vec![0u64; size]);
            for read in chunk {
                a.n_reads += 1;
                for_each_kmer(read, k, |p, fwd, _| {
                    counts[fwd as usize] += 1;
                    if p >= a.windows_at.len() {
                        a.windows_at.resize(p + 1, 0);
                    }
                    a.windows_at[p] += 1;
                });
                for_each_kmer(read, 3, |_, fwd, _| a.tri[fwd as usize] += 1);
                for_each_kmer(read, 2, |_, fwd, _| a.di[fwd as usize] += 1);
            }
        });
    })?;
    let mut c = Counts {
        by_file: vec![vec![0u64; size]; paths.len()],
        tri: vec![0; 64],
        di: vec![0; 16],
        windows_at: Vec::new(),
        n_reads: 0,
    };
    for slot in slots {
        let a = slot.into_inner().unwrap();
        for (fi, counts) in a.by_file.into_iter().enumerate() {
            if let Some(counts) = counts {
                c.by_file[fi]
                    .par_iter_mut()
                    .zip(counts.par_iter())
                    .for_each(|(x, y)| *x += y);
            }
        }
        for (x, y) in c.tri.iter_mut().zip(&a.tri) {
            *x += y;
        }
        for (x, y) in c.di.iter_mut().zip(&a.di) {
            *x += y;
        }
        if a.windows_at.len() > c.windows_at.len() {
            c.windows_at.resize(a.windows_at.len(), 0);
        }
        for (x, y) in c.windows_at.iter_mut().zip(&a.windows_at) {
            *x += y;
        }
        c.n_reads += a.n_reads;
    }
    Ok(c)
}

/// up to max_contexts reads holding each seed k-mer (either strand), oriented so the seed
/// reads forward, as (read, seed position), taken in file order (so the same every run);
/// stops reading once all are full
pub fn contexts(
    paths: &[String],
    k: usize,
    trimmer: Option<&Trimmer>,
    seeds: &[Vec<u8>],
    max_contexts: usize,
    progress: Prog,
    max_reads: u64,
) -> io::Result<Contexts> {
    let mut lut: FxHashMap<u64, usize> = FxHashMap::default();
    let mut fwd_code = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let mut code = None;
        for_each_kmer(seed, k, |p, fwd, rev| {
            if p == 0 {
                code = Some((fwd, rev));
            }
        });
        let (fwd, rev) = code.unwrap_or((u64::MAX, u64::MAX));
        fwd_code.push(fwd);
        lut.entry(fwd).or_insert(i);
        lut.entry(rev).or_insert(i);
    }
    let mut out: Contexts = vec![Vec::new(); seeds.len()];
    let mut open = seeds.len();
    let mut full = vec![false; seeds.len()];
    for path in paths {
        if open == 0 {
            break;
        }
        file_batches(path, trimmer, progress, max_reads, |batch| {
            // finding the seeds in each read in parallel, then taking them in read order
            let found: Vec<Vec<(usize, usize, bool)>> = batch
                .par_iter()
                .map(|read| {
                    let mut hits = Vec::new();
                    for_each_kmer(read, k, |p, fwd, _| {
                        if let Some(&w) = lut.get(&fwd) {
                            hits.push((w, p, fwd == fwd_code[w]));
                        }
                    });
                    hits
                })
                .collect();
            for (read, hits) in batch.iter().zip(found) {
                for (w, p, is_fwd) in hits {
                    if full[w] {
                        continue;
                    }
                    if is_fwd {
                        out[w].push((read.clone(), p));
                    } else {
                        out[w].push((revcomp(read), read.len() - p - k));
                    }
                    if out[w].len() >= max_contexts {
                        full[w] = true;
                        open -= 1;
                    }
                }
            }
            open > 0
        })?;
    }
    Ok(out)
}
