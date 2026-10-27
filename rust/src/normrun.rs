//! The two passes of `bit normalize-reads` over the read files (see
//! bit/modules/normalize_reads.py for the approach): counting every read's anchors into a
//! sketch, then judging and writing the reads with those counts known. Reads are read on a
//! thread of their own, judged a batch at a time on the threads set with set_threads, and kept
//! ones written as they were in the input (names, comments, qualities and all), gzipped in
//! parallel.

use crate::anchors::{anchors_into, threshold};
use crate::reads::{Progress, SeqReader};
use crate::sketch::Sketch;
use flate2::Compression;
use flate2::write::GzEncoder;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::fs::File;
use std::hash::Hasher;
use std::io::{self, BufWriter, Write};
use std::sync::Arc;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{Receiver, sync_channel};

// reads (or pairs) judged at a time (each is judged on its own: batches change nothing but speed)
const BATCH: usize = 8192;
// an anchor in under this share of the reads around it (within 2k bases) is taken for a
// recurrent sequencing error, not sequence of its own (see judge_read())
const ERROR_SHARE: f64 = 0.01;
// output is gzipped this many bytes at a time, each its own gzip member (a file of several
// members is still one gzip file to every reader), so it can be done on several threads
const GZ_CHUNK: usize = 4 << 20;
const GZ_PENDING: usize = 8 * GZ_CHUNK;
// the deepest depth kept apart in the depth histograms (deeper goes in the last bin)
const MAX_DEPTH_BIN: usize = 1 << 16;

fn to_py_err(e: io::Error) -> PyErr {
    if e.kind() == io::ErrorKind::InvalidData {
        PyValueError::new_err(e.to_string())
    } else {
        PyIOError::new_err(e.to_string())
    }
}

fn bad(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// a read, or a pair: sequences, and with `raw` the records as in the files
#[derive(Default)]
struct Rec {
    seq1: Vec<u8>,
    raw1: Vec<u8>,
    seq2: Option<Vec<u8>>,
    raw2: Option<Vec<u8>>,
}

/// a record's name: its header's first word, less the marker and any /1 or /2
fn base_name(raw: &[u8]) -> &[u8] {
    let line = raw.get(1..).unwrap_or(&[]);
    let end = line
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(line.len());
    let n = &line[..end];
    if n.len() >= 2 && n[n.len() - 2] == b'/' && matches!(n[n.len() - 1], b'1' | b'2') {
        &n[..n.len() - 2]
    } else {
        n
    }
}

fn name_str(raw: &[u8]) -> String {
    let line = raw.get(1..).unwrap_or(&[]);
    let end = line
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(line.len());
    String::from_utf8_lossy(&line[..end]).into_owned()
}

/// batches of reads (or pairs, the two files read in step, their names checked to match),
/// read on a thread of their own, a batch ahead
fn spawn_reader(
    path1: String,
    path2: Option<String>,
    batch_size: usize,
    progress: Option<Arc<Progress>>,
) -> io::Result<(Receiver<io::Result<Vec<Rec>>>, bool)> {
    let mut r1 = SeqReader::open(&path1, progress.clone())?;
    let mut r2 = match &path2 {
        Some(p) => Some(SeqReader::open(p, progress.clone())?),
        None => None,
    };
    let fastq = r1.is_fastq();
    let (tx, rx) = sync_channel(2);
    std::thread::spawn(move || {
        let result = (|| -> io::Result<()> {
            loop {
                let mut batch = Vec::with_capacity(batch_size);
                while batch.len() < batch_size {
                    let mut rec = Rec::default();
                    if !r1.next_record(&mut rec.seq1, &mut rec.raw1)? {
                        break;
                    }
                    if let Some(r2) = r2.as_mut() {
                        let (mut s, mut raw) = (Vec::new(), Vec::new());
                        if !r2.next_record(&mut s, &mut raw)? {
                            return Err(bad(format!(
                                "Read-2 file ran out before read-1 file (at read '{}').",
                                name_str(&rec.raw1)
                            )));
                        }
                        if base_name(&rec.raw1) != base_name(&raw) {
                            return Err(bad(format!(
                                "Read names don't match between read files: '{}' vs '{}'. \
                                 The read files need to be in the same order.",
                                name_str(&rec.raw1),
                                name_str(&raw)
                            )));
                        }
                        rec.seq2 = Some(s);
                        rec.raw2 = Some(raw);
                    }
                    batch.push(rec);
                }
                let done = batch.len() < batch_size;
                if done
                    && let Some(r2) = r2.as_mut()
                    && r2.next_seq(&mut Vec::new())?
                {
                    return Err(bad("Read-1 file ran out before read-2 file.".to_string()));
                }
                if let Some(p) = &progress {
                    p.reads.fetch_add(batch.len() as u64, Ordering::Relaxed);
                }
                if tx.send(Ok(batch)).is_err() || done {
                    return Ok(());
                }
            }
        })();
        if let Err(e) = result {
            let _ = tx.send(Err(e));
        }
    });
    Ok((rx, fastq))
}

/// a gzipped output, compressed GZ_CHUNK at a time on several threads (see GZ_CHUNK)
struct GzOut {
    file: BufWriter<File>,
    pending: Vec<u8>,
    wrote: bool,
}

impl GzOut {
    fn create(path: &str) -> io::Result<Self> {
        Ok(GzOut {
            file: BufWriter::with_capacity(1 << 20, File::create(path)?),
            pending: Vec::with_capacity(GZ_PENDING + GZ_CHUNK),
            wrote: false,
        })
    }

    fn push(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.pending.extend_from_slice(bytes);
        if self.pending.len() >= GZ_PENDING {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let members: Vec<io::Result<Vec<u8>>> = self
            .pending
            .par_chunks(GZ_CHUNK)
            .map(|chunk| {
                let mut gz = GzEncoder::new(Vec::with_capacity(chunk.len() / 3), Compression::new(1));
                gz.write_all(chunk)?;
                gz.finish()
            })
            .collect();
        for m in members {
            self.file.write_all(&m?)?;
        }
        self.pending.clear();
        self.wrote = true;
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        self.flush()?;
        if !self.wrote {
            // an empty output is still a valid (empty) gzip file
            self.file.write_all(&GzEncoder::new(Vec::new(), Compression::new(1)).finish()?)?;
        }
        self.file.flush()
    }
}

/// where kept reads go: one file, or two for pairs, and an optional depth table
struct Outputs {
    r1: GzOut,
    r2: Option<GzOut>,
    depths: Option<GzOut>,
}

impl Outputs {
    fn create(out1: &str, out2: Option<&str>, depth_table: Option<&str>) -> io::Result<Self> {
        let mut depths = match depth_table {
            Some(p) => Some(GzOut::create(p)?),
            None => None,
        };
        if let Some(d) = depths.as_mut() {
            d.push(b"read\tdepth\n")?;
        }
        Ok(Outputs {
            r1: GzOut::create(out1)?,
            r2: match out2 {
                Some(p) => Some(GzOut::create(p)?),
                None => None,
            },
            depths,
        })
    }

    fn write(&mut self, rec: &Rec, depth: u32) -> io::Result<()> {
        self.r1.push(&rec.raw1)?;
        if let (Some(o), Some(raw)) = (self.r2.as_mut(), rec.raw2.as_ref()) {
            o.push(raw)?;
        }
        if let Some(d) = self.depths.as_mut() {
            d.push(base_name(&rec.raw1))?;
            d.push(format!("\t{depth}\n").as_bytes())?;
        }
        Ok(())
    }

    fn finish(self) -> io::Result<()> {
        self.r1.finish()?;
        if let Some(o) = self.r2 {
            o.finish()?;
        }
        if let Some(d) = self.depths {
            d.finish()?;
        }
        Ok(())
    }
}

/// running totals, for the summary
#[derive(Default)]
struct Totals {
    reads_in: u64,
    bases_in: u64,
    reads_out: u64,
    bases_out: u64,
    no_anchors: u64,
    kept_no_solid: u64,
    anchor_occurrences: u64,
    error_occurrences: u64,
    depth_in: Vec<u64>,  // reads (or pairs) by typical depth
    depth_out: Vec<u64>, // kept ones by what their depth comes to after normalization
}

impl Totals {
    fn new() -> Self {
        Totals {
            depth_in: vec![0; MAX_DEPTH_BIN],
            depth_out: vec![0; MAX_DEPTH_BIN],
            ..Default::default()
        }
    }

    fn median(hist: &[u64]) -> u64 {
        let n: u64 = hist.iter().sum();
        if n == 0 {
            return 0;
        }
        let mut seen = 0;
        for (d, c) in hist.iter().enumerate() {
            seen += c;
            if 2 * seen >= n {
                return d as u64;
            }
        }
        0
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("reads_or_pairs_in", self.reads_in)?;
        d.set_item("bases_in", self.bases_in)?;
        d.set_item("reads_or_pairs_out", self.reads_out)?;
        d.set_item("bases_out", self.bases_out)?;
        d.set_item("no_anchors", self.no_anchors)?;
        d.set_item("kept_without_solid_kmers", self.kept_no_solid)?;
        d.set_item("kmer_occurrences", self.anchor_occurrences)?;
        d.set_item("kmer_occurrences_judged_errors", self.error_occurrences)?;
        d.set_item("median_depth_in", Self::median(&self.depth_in))?;
        d.set_item("median_depth_out", Self::median(&self.depth_out))?;
        Ok(d)
    }
}

fn bin(d: f64) -> usize {
    (d.round().max(0.0) as usize).min(MAX_DEPTH_BIN - 1)
}

fn bases(rec: &Rec) -> u64 {
    (rec.seq1.len() + rec.seq2.as_ref().map_or(0, |s| s.len())) as u64
}

/// a read's (or pair's) distinct anchors, read 1's and read 2's together
fn distinct_anchors(rec: &Rec, k: usize, t: u64, along: &mut Vec<(u32, u64)>, out: &mut Vec<u64>) {
    out.clear();
    anchors_into(&rec.seq1, k, t, along);
    out.extend(along.iter().map(|a| a.1));
    if let Some(s2) = &rec.seq2 {
        anchors_into(s2, k, t, along);
        out.extend(along.iter().map(|a| a.1));
    }
    out.sort_unstable();
    out.dedup();
}

/// The read counts of every anchor in the reads, from the counting pass, and the judging
/// and writing pass that uses them.
#[pyclass(frozen)]
pub struct AnchorCounts {
    sketch: Arc<Sketch>,
    k: usize,
    density: u64,
    reads: u64,
    bases: u64,
}

/// the counting pass: every read's (or pair's) distinct anchors into a sketch of about
/// `sketch_bytes` (see sketch.rs)
#[pyfunction]
#[pyo3(signature = (read_1, read_2, k, density, sketch_bytes, progress=None))]
pub fn count_anchors(
    py: Python<'_>,
    read_1: String,
    read_2: Option<String>,
    k: usize,
    density: u64,
    sketch_bytes: usize,
    progress: Option<PyRef<crate::Progress>>,
) -> PyResult<AnchorCounts> {
    if !(1..=32).contains(&k) {
        return Err(PyValueError::new_err("k must be between 1 and 32"));
    }
    let pr = crate::arc(&progress);
    py.detach(|| {
        let sketch = Sketch::new(sketch_bytes);
        let t = threshold(density);
        let (rx, _) = spawn_reader(read_1, read_2, BATCH, pr)?;
        let (mut reads, mut nbases) = (0u64, 0u64);
        for batch in rx {
            let batch = batch?;
            reads += batch.len() as u64;
            nbases += batch.iter().map(bases).sum::<u64>();
            batch.par_iter().for_each_init(
                || (Vec::new(), Vec::new()),
                |(along, hs), rec| {
                    distinct_anchors(rec, k, t, along, hs);
                    for &h in hs.iter() {
                        sketch.add(h);
                    }
                },
            );
        }
        Ok(AnchorCounts {
            sketch: Arc::new(sketch),
            k,
            density,
            reads,
            bases: nbases,
        })
    })
    .map_err(to_py_err)
}

/// what judging a read (or pair) with known counts comes to
struct Call {
    keep_prob: f64,
    typical: u32,
    occurrences: u64,
    errors: u64,
    has_anchors: bool,
    has_solid: bool,
}

#[derive(Clone, Copy)]
struct Rules {
    k: usize,
    t: u64,
    target: f64,
    novel_fraction: f64,
    minor_ratio: f64,
    minor_reads: u32,
    minor_anchors: usize,
    minor_cap: f64,
    min_reads: u32,
}

/// the median of `v` (sorted in place), 0 for none
fn median(v: &mut [u32]) -> u32 {
    if v.is_empty() {
        return 0;
    }
    let mid = v.len() / 2;
    *v.select_nth_unstable(mid).1
}

/// the upper quartile of `v` (reordered in place), 0 for none
fn upper_quartile(v: &mut [u32]) -> u32 {
    if v.is_empty() {
        return 0;
    }
    let at = (3 * v.len()) / 4;
    *v.select_nth_unstable(at.min(v.len() - 1)).1
}

/// a read's anchors as (position, read count), in read order; each is solid when it's in at
/// least min_reads reads and at least ERROR_SHARE of the upper quartile of those within 2k
/// bases of it (an error's neighbors are the same sequence without it). The upper quartile
/// rather than the median: two errors close together put about half the anchors around
/// either one on an error, which brought the median down to the errors' own counts, so at
/// 1,000x (where an error recurs in 2-3 reads) reads with two or three errors were kept at 3-6
/// times the rate of clean ones (Oct 2026). Adds (hash, count, solid) for each to `out`
fn judge_read(seq: &[u8], r: &Rules, sketch: &Sketch, along: &mut Vec<(u32, u64)>,
              out: &mut Vec<(u64, u32, bool)>) {
    anchors_into(seq, r.k, r.t, along);
    let counts: Vec<u32> = along.iter().map(|a| sketch.get(a.1)).collect();
    // only an anchor under ERROR_SHARE of the read's deepest one can be under that share of
    // those around it, so most skip the closer look
    let deepest = counts.iter().copied().max().unwrap_or(0) as f64;
    let reach = 2 * r.k as u32;
    let mut lo = 0;
    let mut near = Vec::new();
    for (i, &(p, h)) in along.iter().enumerate() {
        while along[lo].0 + reach < p {
            lo += 1;
        }
        let c = counts[i];
        if (c as f64) >= ERROR_SHARE * deepest {
            out.push((h, c, c >= r.min_reads));
            continue;
        }
        near.clear();
        near.extend(
            along[lo..]
                .iter()
                .zip(&counts[lo..])
                .take_while(|(a, _)| a.0 <= p + reach)
                .map(|(_, c)| *c),
        );
        let local = upper_quartile(&mut near) as f64;
        out.push((h, c, c >= r.min_reads && (c as f64) >= ERROR_SHARE * local));
    }
}

/// How likely a read (or pair) is to be kept, from its anchors' read counts in the whole input
/// (both mates' together). Anchors that aren't solid (see judge_read) are errors and count for
/// nothing. The read's depth is that of its solid anchors at the novel_fraction quantile (the
/// 10% least covered, by default), not their median, so a read is judged by its least-covered
/// part; it's kept with probability target / that depth (always, at or under the target).
/// With minor_anchors of its solid anchors on a minor branch (in minor_reads+ reads, at
/// minor_ratio to half the read's typical depth), it's kept with at least probability
/// minor_cap / the branch's depth. A read with anchors but none solid is in no more than one
/// or two reads' worth of anything, and is kept.
fn call(rec: &Rec, r: &Rules, sketch: &Sketch, along: &mut Vec<(u32, u64)>,
        all: &mut Vec<(u64, u32, bool)>) -> Call {
    all.clear();
    judge_read(&rec.seq1, r, sketch, along, all);
    if let Some(s2) = &rec.seq2 {
        judge_read(s2, r, sketch, along, all);
    }
    let occurrences = all.len() as u64;
    let errors = all.iter().filter(|a| !a.2).count() as u64;
    if all.is_empty() {
        return Call { keep_prob: 0.0, typical: 0, occurrences, errors, has_anchors: false, has_solid: false };
    }
    // distinct anchors (a stretch both mates hold counts once), solid if any copy is
    all.sort_unstable_by_key(|a| (a.0, !a.2));
    all.dedup_by_key(|a| a.0);
    let mut solid: Vec<u32> = all.iter().filter(|a| a.2).map(|a| a.1).collect();
    if solid.is_empty() {
        return Call { keep_prob: 1.0, typical: 0, occurrences, errors, has_anchors: true, has_solid: false };
    }
    solid.sort_unstable();
    let typical = solid[solid.len() / 2];
    let at = ((r.novel_fraction * solid.len() as f64).ceil() as usize).clamp(1, solid.len()) - 1;
    let low = solid[at].max(1) as f64;
    let mut p = (r.target / low).min(1.0);
    if r.minor_anchors > 0 {
        let lo_d = (r.minor_ratio * typical as f64).max(r.minor_reads as f64);
        let hi_d = 0.5 * typical as f64;
        let mut branch: Vec<u32> = solid
            .iter()
            .copied()
            .filter(|&c| c as f64 >= lo_d && c as f64 <= hi_d)
            .collect();
        if branch.len() >= r.minor_anchors {
            let d = median(&mut branch).max(1) as f64;
            p = p.max((r.minor_cap / d).min(1.0));
        }
    }
    Call { keep_prob: p, typical, occurrences, errors, has_anchors: true, has_solid: true }
}

/// a number in [0, 1) from the read's own sequence and `seed`, the same whatever the order
/// or threads
fn coin(rec: &Rec, seed: u64) -> f64 {
    let mut h = rustc_hash::FxHasher::default();
    h.write_u64(seed);
    h.write(&rec.seq1);
    if let Some(s2) = &rec.seq2 {
        h.write_u8(0);
        h.write(s2);
    }
    let mut x = h.finish() ^ seed.rotate_left(17);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^= x >> 31;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

#[pymethods]
impl AnchorCounts {
    /// The judging and writing pass (see call()). Writes kept reads to `out_1` (and their
    /// mates to `out_2`), as they were in the input and in input order, and with
    /// `depth_table` each kept read's name and typical depth; returns the totals for the
    /// summary.
    #[pyo3(signature = (read_1, read_2, out_1, out_2, target_depth, novel_fraction=0.1,
                        minor_ratio=0.05, minor_reads=3, minor_anchors=3, min_reads=2, seed=1,
                        depth_table=None, progress=None))]
    #[allow(clippy::too_many_arguments)]
    fn normalize<'py>(
        &self,
        py: Python<'py>,
        read_1: String,
        read_2: Option<String>,
        out_1: String,
        out_2: Option<String>,
        target_depth: u32,
        novel_fraction: f64,
        minor_ratio: f64,
        minor_reads: u32,
        minor_anchors: usize,
        min_reads: u32,
        seed: u64,
        depth_table: Option<String>,
        progress: Option<PyRef<crate::Progress>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        if target_depth == 0 {
            return Err(PyValueError::new_err("target_depth must be 1 or greater"));
        }
        let rules = Rules {
            k: self.k,
            t: threshold(self.density),
            target: target_depth as f64,
            novel_fraction,
            minor_ratio,
            minor_reads,
            minor_anchors,
            minor_cap: (target_depth / 2).max(10).min(target_depth) as f64,
            min_reads,
        };
        let pr = crate::arc(&progress);
        let sketch = self.sketch.clone();
        let totals = py
            .detach(|| -> io::Result<Totals> {
                let (rx, _) = spawn_reader(read_1, read_2, BATCH, pr)?;
                let mut outs = Outputs::create(&out_1, out_2.as_deref(), depth_table.as_deref())?;
                let mut tot = Totals::new();
                for batch in rx {
                    let batch = batch?;
                    let calls: Vec<Call> = batch
                        .par_iter()
                        .map_init(
                            || (Vec::new(), Vec::new()),
                            |(along, all), rec| call(rec, &rules, &sketch, along, all),
                        )
                        .collect();
                    for (rec, c) in batch.iter().zip(&calls) {
                        tot.reads_in += 1;
                        tot.bases_in += bases(rec);
                        tot.anchor_occurrences += c.occurrences;
                        tot.error_occurrences += c.errors;
                        if !c.has_anchors {
                            tot.no_anchors += 1;
                            continue;
                        }
                        tot.depth_in[bin(c.typical as f64)] += 1;
                        if c.keep_prob >= 1.0 || coin(rec, seed) < c.keep_prob {
                            tot.reads_out += 1;
                            tot.bases_out += bases(rec);
                            tot.kept_no_solid += !c.has_solid as u64;
                            tot.depth_out[bin(c.typical as f64 * c.keep_prob)] += 1;
                            outs.write(rec, c.typical)?;
                        }
                    }
                }
                outs.finish()?;
                Ok(tot)
            })
            .map_err(to_py_err)?;
        let d = totals.to_dict(py)?;
        d.set_item("sketch_bytes", self.sketch.bytes())?;
        d.set_item("sketch_fill", self.sketch.fill())?;
        Ok(d)
    }

    #[getter]
    fn reads(&self) -> u64 {
        self.reads
    }

    #[getter]
    fn bases(&self) -> u64 {
        self.bases
    }

    /// (how many reads hold each of these sequences' anchors), for checking: the sketch's
    /// count for every anchor of each sequence, in order
    fn counts_of(&self, seq: &str) -> Vec<u32> {
        let mut along = Vec::new();
        anchors_into(seq.as_bytes(), self.k, threshold(self.density), &mut along);
        along.iter().map(|a| self.sketch.get(a.1)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_names() {
        assert_eq!(base_name(b"@r1/1 extra\n"), b"r1");
        assert_eq!(base_name(b"@r1/2\n"), b"r1");
        assert_eq!(base_name(b">contig_3 len=4\n"), b"contig_3");
        assert_eq!(base_name(b"@SRR1.1 1/1\n"), b"SRR1.1");
    }

    #[test]
    fn coin_is_order_free_and_spread() {
        let a = Rec { seq1: b"ACGTACGTTT".to_vec(), ..Default::default() };
        assert_eq!(coin(&a, 1), coin(&a, 1));
        assert_ne!(coin(&a, 1), coin(&a, 2));
        let n = 10_000;
        let below = (0..n)
            .filter(|i| {
                let r = Rec { seq1: format!("ACGT{i}").into_bytes(), ..Default::default() };
                coin(&r, 1) < 0.25
            })
            .count();
        assert!((2_200..2_800).contains(&below));
    }
}
