//! A small BAM reader (records in file order, as pysam's fetch(until_eof=True) gives them) and
//! the per-read passes of `bit mapped-read-stats` and `bit cov-stats` (see get_mapped_reads_pids
//! and get_per_contig_pid_stats in bit/modules/mapped_read_stats.py, which this follows
//! step for step), with cov-stats' per-base coverage counted in the same pass. BGZF blocks are decompressed in parallel.
//!
//! Anything unexpected (not BGZF, a malformed record, an NM tag that isn't an integer, a mapped
//! read without a CIGAR or reference) is returned as an error, and the caller falls back to the
//! Python (pysam) path, so results (or errors) stay exactly as they were.

use flate2::{Decompress, FlushDecompress};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufReader, Read};

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// reads BGZF: a thread of its own reads the file a batch of blocks at a time and inflates
/// each batch in parallel, a couple of batches ahead of what's being read, so inflating
/// and whatever is done with the data overlap rather than take turns
struct Bgzf {
    batches: std::sync::mpsc::Receiver<io::Result<Vec<Vec<u8>>>>,
    // inflated blocks not read yet, and how far into the first the reading is (blocks are
    // handed over whole and copied out a slice at a time: going through a queue of bytes
    // copied each one twice, one by one)
    out: VecDeque<Vec<u8>>,
    at: usize,
    done: bool,
}

/// the next raw block's compressed deflate data, expected size and CRC, or None at the end
fn next_block(raw: &mut impl Read) -> io::Result<Option<(Vec<u8>, usize, u32)>> {
    let mut head = [0u8; 12];
    let mut got = 0;
    while got < 12 {
        let n = raw.read(&mut head[got..])?;
        if n == 0 {
            return if got == 0 { Ok(None) } else { Err(bad("truncated BGZF block")) };
        }
        got += n;
    }
    if head[0] != 31 || head[1] != 139 || head[2] != 8 || head[3] & 4 == 0 {
        return Err(bad("not BGZF"));
    }
    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    raw.read_exact(&mut extra)?;
    let mut bsize = None;
    let mut i = 0;
    while i + 4 <= xlen {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 && i + 6 <= xlen {
            bsize = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize);
        }
        i += 4 + slen;
    }
    let bsize = bsize.ok_or_else(|| bad("not BGZF"))?;
    let cdata_len = (bsize + 1)
        .checked_sub(xlen + 20)
        .ok_or_else(|| bad("bad BGZF block size"))?;
    let mut cdata = vec![0u8; cdata_len];
    raw.read_exact(&mut cdata)?;
    let mut tail = [0u8; 8];
    raw.read_exact(&mut tail)?;
    let crc = u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]);
    let isize = u32::from_le_bytes([tail[4], tail[5], tail[6], tail[7]]) as usize;
    Ok(Some((cdata, isize, crc)))
}

/// the next batch of blocks, inflated (the empty ones left out), and whether it's the last
fn next_batch(raw: &mut impl Read) -> io::Result<(Vec<Vec<u8>>, bool)> {
    let mut blocks = Vec::with_capacity(Bgzf::BATCH);
    let mut last = false;
    while blocks.len() < Bgzf::BATCH {
        match next_block(raw)? {
            Some(b) => blocks.push(b),
            None => {
                last = true;
                break;
            }
        }
    }
    let inflated: Vec<io::Result<Vec<u8>>> = blocks
        .par_iter()
        .map(|(cdata, isize, crc)| {
            let mut out = vec![0u8; *isize];
            let mut d = Decompress::new(false);
            d.decompress(cdata, &mut out, FlushDecompress::Finish)
                .map_err(|_| bad("bad BGZF deflate data"))?;
            if d.total_out() as usize != *isize {
                return Err(bad("bad BGZF block"));
            }
            // (htslib checks each block's CRC, and pysam then errors)
            let mut c = flate2::Crc::new();
            c.update(&out);
            if c.sum() != *crc {
                return Err(bad("BGZF CRC mismatch"));
            }
            Ok(out)
        })
        .collect();
    let mut out = Vec::with_capacity(inflated.len());
    for b in inflated {
        let b = b?;
        if !b.is_empty() {
            out.push(b);
        }
    }
    Ok((out, last))
}

impl Bgzf {
    const BATCH: usize = 256;

    fn open(path: &str) -> io::Result<Self> {
        let mut raw = BufReader::with_capacity(1 << 20, File::open(path)?);
        let (tx, batches) = std::sync::mpsc::sync_channel(2);
        std::thread::spawn(move || {
            loop {
                let batch = next_batch(&mut raw);
                let stop = !matches!(batch, Ok((_, false)));
                // (an error from send: the reader is gone, so there's no one to read for)
                if tx.send(batch.map(|(blocks, _)| blocks)).is_err() || stop {
                    break;
                }
            }
        });
        Ok(Bgzf {
            batches,
            out: VecDeque::new(),
            at: 0,
            done: false,
        })
    }
}

impl Read for Bgzf {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.out.is_empty() && !self.done {
            match self.batches.recv() {
                Ok(batch) => self.out.extend(batch?),
                // the reading thread is done: the file's end, or an error already passed on
                Err(_) => self.done = true,
            }
        }
        let Some(block) = self.out.front() else { return Ok(0) };
        let n = buf.len().min(block.len() - self.at);
        buf[..n].copy_from_slice(&block[self.at..self.at + n]);
        self.at += n;
        if self.at == block.len() {
            self.out.pop_front();
            self.at = 0;
        }
        Ok(n)
    }
}

fn read_i32(r: &mut impl Read) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

pub struct Record {
    pub flag: u16,
    pub ref_id: i32,
    pub pos: i32,
    pub name: Vec<u8>,
    pub cigar: Vec<u32>,
    pub nm: Option<i64>,
}

pub struct Bam {
    inner: Bgzf,
    pub ref_names: Vec<String>,
    pub ref_lens: Vec<usize>,
    buf: Vec<u8>,
}

/// size in bytes of one value of an aux type
fn aux_size(t: u8) -> Option<usize> {
    match t {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        _ => None,
    }
}

struct Aux<'a> {
    nm: Option<Result<i64, ()>>,
    // the first CG tag: Some(its values) if it's B:I/B:i, Some(None) if of another type
    cg: Option<Option<Vec<u32>>>,
    _d: std::marker::PhantomData<&'a ()>,
}

/// finds NM (the first, as htslib's bam_aux_get does) and CG (a long CIGAR) in the aux data
fn scan_aux(mut a: &[u8]) -> io::Result<Aux<'_>> {
    let mut out = Aux {
        nm: None,
        cg: None,
        _d: std::marker::PhantomData,
    };
    while !a.is_empty() {
        if a.len() < 3 {
            return Err(bad("bad aux data"));
        }
        let (tag, t) = ([a[0], a[1]], a[2]);
        a = &a[3..];
        match t {
            b'Z' | b'H' => {
                let end = memchr::memchr(0, a).ok_or_else(|| bad("bad aux data"))?;
                if &tag == b"NM" && out.nm.is_none() {
                    out.nm = Some(Err(()));
                }
                if &tag == b"CG" && out.cg.is_none() {
                    out.cg = Some(None);
                }
                a = &a[end + 1..];
            }
            b'B' => {
                if a.len() < 5 {
                    return Err(bad("bad aux data"));
                }
                let sub = a[0];
                let n = u32::from_le_bytes([a[1], a[2], a[3], a[4]]) as usize;
                let size = aux_size(sub).ok_or_else(|| bad("bad aux data"))?;
                let body = a.get(5..5 + n * size).ok_or_else(|| bad("bad aux data"))?;
                if &tag == b"NM" && out.nm.is_none() {
                    out.nm = Some(Err(()));
                }
                if &tag == b"CG" && out.cg.is_none() {
                    out.cg = Some((sub == b'I' || sub == b'i').then(|| {
                        body.chunks_exact(4)
                            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect()
                    }));
                }
                a = &a[5 + n * size..];
            }
            _ => {
                let size = aux_size(t).ok_or_else(|| bad("bad aux data"))?;
                let v = a.get(..size).ok_or_else(|| bad("bad aux data"))?;
                if &tag == b"CG" && out.cg.is_none() {
                    out.cg = Some(None);
                }
                if &tag == b"NM" && out.nm.is_none() {
                    out.nm = Some(match t {
                        b'c' => Ok(v[0] as i8 as i64),
                        b'C' => Ok(v[0] as i64),
                        b's' => Ok(i16::from_le_bytes([v[0], v[1]]) as i64),
                        b'S' => Ok(u16::from_le_bytes([v[0], v[1]]) as i64),
                        b'i' => Ok(i32::from_le_bytes([v[0], v[1], v[2], v[3]]) as i64),
                        b'I' => Ok(u32::from_le_bytes([v[0], v[1], v[2], v[3]]) as i64),
                        _ => Err(()), // A (a character) or f (a float): not an integer
                    });
                }
                a = &a[size..];
            }
        }
    }
    Ok(out)
}

impl Bam {
    pub fn open(path: &str) -> io::Result<Self> {
        let mut inner = Bgzf::open(path)?;
        let mut magic = [0u8; 4];
        inner.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            return Err(bad("not a BAM file"));
        }
        let l_text = read_i32(&mut inner)?;
        if l_text < 0 {
            return Err(bad("bad BAM header"));
        }
        io::copy(&mut (&mut inner).take(l_text as u64), &mut io::sink())?;
        let n_ref = read_i32(&mut inner)?;
        let mut ref_names = Vec::new();
        let mut ref_lens = Vec::new();
        for _ in 0..n_ref.max(0) {
            let l_name = read_i32(&mut inner)?;
            if l_name < 1 {
                return Err(bad("bad BAM header"));
            }
            let mut name = vec![0u8; l_name as usize];
            inner.read_exact(&mut name)?;
            name.pop(); // the NUL
            if !name.is_ascii() {
                return Err(bad("non-ASCII reference name"));
            }
            ref_names.push(String::from_utf8(name).unwrap());
            let l_ref = read_i32(&mut inner)?;
            if l_ref < 0 {
                return Err(bad("bad BAM header"));
            }
            ref_lens.push(l_ref as usize);
        }
        Ok(Bam {
            inner,
            ref_names,
            ref_lens,
            buf: Vec::new(),
        })
    }

    /// the next record (None at the end)
    pub fn next(&mut self) -> io::Result<Option<Record>> {
        let mut size = [0u8; 4];
        let mut got = 0;
        while got < 4 {
            let n = self.inner.read(&mut size[got..])?;
            if n == 0 {
                return if got == 0 { Ok(None) } else { Err(bad("truncated BAM")) };
            }
            got += n;
        }
        let size = i32::from_le_bytes(size);
        if size < 32 {
            return Err(bad("bad BAM record"));
        }
        self.buf.resize(size as usize, 0);
        self.inner.read_exact(&mut self.buf)?;
        let b = &self.buf;
        let i32_at = |o: usize| i32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]]);
        let u16_at = |o: usize| u16::from_le_bytes([b[o], b[o + 1]]);
        let ref_id = i32_at(0);
        let mate_ref_id = i32_at(20);
        let pos = i32_at(4);
        let l_name = b[8] as usize;
        let n_cigar = u16_at(12) as usize;
        let flag = u16_at(14);
        let l_seq = i32_at(16);
        if l_seq < 0 || l_name < 1 {
            return Err(bad("bad BAM record"));
        }
        let l_seq = l_seq as usize;
        let name_end = 32 + l_name;
        let cigar_end = name_end + 4 * n_cigar;
        let aux_start = cigar_end + l_seq.div_ceil(2) + l_seq;
        if aux_start > b.len() || b[name_end - 1] != 0 {
            return Err(bad("bad BAM record"));
        }
        // (htslib may pad the name with extra NULs; pysam's query_name stops at the first)
        let name_bytes = &b[32..name_end - 1];
        let name = name_bytes[..memchr::memchr(0, name_bytes).unwrap_or(name_bytes.len())].to_vec();
        let mut cigar: Vec<u32> = b[name_end..cigar_end]
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        let aux = scan_aux(&b[aux_start..])?;
        // a long CIGAR held in CG, as htslib swaps in when reading (bam_tag2cigar)
        if !cigar.is_empty()
            && ref_id >= 0
            && pos >= 0
            && cigar[0] & 0xf == 4
            && (cigar[0] >> 4) as usize == l_seq
            && let Some(Some(cg)) = aux.cg
            && cg.len() >= cigar.len()
            && cg.len() < (1 << 29)
        {
            cigar = cg;
        }
        // records htslib rejects (pysam then raises): reference or mate ids out of range,
        // unknown CIGAR operations, and a CIGAR whose query length differs from the sequence's
        let n_ref = self.ref_names.len() as i32;
        if !(-1..n_ref).contains(&ref_id) || !(-1..n_ref).contains(&mate_ref_id) {
            return Err(bad("reference id out of range"));
        }
        if cigar.iter().any(|&c| c & 0xf > 8) {
            return Err(bad("unknown CIGAR operation"));
        }
        if !cigar.is_empty() && l_seq > 0 && flag & 0x4 == 0 {
            let qlen: u64 = cigar
                .iter()
                .filter(|&&c| matches!(c & 0xf, 0 | 1 | 4 | 7 | 8))
                .map(|&c| (c >> 4) as u64)
                .sum();
            if qlen != l_seq as u64 {
                return Err(bad("CIGAR and query sequence lengths differ"));
            }
        }
        let nm = match aux.nm {
            None => None,
            Some(Ok(v)) => Some(v),
            Some(Err(())) => return Err(bad("NM tag that isn't an integer")),
        };
        Ok(Some(Record {
            flag,
            ref_id,
            pos,
            name,
            cigar,
            nm,
        }))
    }
}

/// parse_cigar_breakdown()
#[derive(Clone, Copy)]
pub struct Breakdown {
    pub aligned_cols: i64,
    pub query_aligned: i64,
    pub matches_mismatches: i64,
    pub gap_events: i64,
    pub soft_clipped: i64,
    pub hard_clipped: i64,
    pub clipped: i64,
    pub full_read_length: i64,
}

fn breakdown(cigar: &[u32]) -> Breakdown {
    let mut counts = [0i64; 16];
    let (mut ins, mut del) = (0i64, 0i64);
    for &c in cigar {
        let (op, len) = ((c & 0xf) as usize, (c >> 4) as i64);
        counts[op] += len;
        if op == 1 {
            ins += 1;
        } else if op == 2 {
            del += 1;
        }
    }
    let mm = counts[0] + counts[7] + counts[8];
    Breakdown {
        aligned_cols: mm + counts[1] + counts[2],
        query_aligned: mm + counts[1],
        matches_mismatches: mm,
        gap_events: ins + del,
        soft_clipped: counts[4],
        hard_clipped: counts[5],
        clipped: counts[4] + counts[5],
        full_read_length: mm + counts[1] + counts[4] + counts[5],
    }
}

/// compute_pids(): (gap-aware, gap-compressed), with Python's int / int -> float
fn pids(b: &Breakdown, nm: i64) -> (f64, f64) {
    let matches = b.aligned_cols - nm;
    let ga = if b.aligned_cols != 0 {
        matches as f64 / b.aligned_cols as f64 * 100.0
    } else {
        0.0
    };
    let gc_denom = b.matches_mismatches + b.gap_events;
    let gc = if gc_denom != 0 {
        matches as f64 / gc_denom as f64 * 100.0
    } else {
        0.0
    };
    (ga, gc)
}

/// PidStats, minus the rounded histogram: counts per distinct value instead (Python rounds
/// them into its histogram, so its own round() is used). The total is summed in file order,
/// as Python does, so it's the same float.
#[derive(Default)]
pub struct Acc {
    pub count: u64,
    pub total: f64,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub values: FxHashMap<u64, u64>, // value's bits -> count
    pub order: Vec<f64>,             // distinct values in first-seen order
}

impl Acc {
    fn add(&mut self, v: f64) {
        self.count += 1;
        self.total += v;
        if self.min.is_none_or(|m| v < m) {
            self.min = Some(v);
        }
        if self.max.is_none_or(|m| v > m) {
            self.max = Some(v);
        }
        let e = self.values.entry(v.to_bits()).or_insert(0);
        if *e == 0 {
            self.order.push(v);
        }
        *e += 1;
    }
}

/// ClipAggregate (length_stats are of integers, so held as such)
#[derive(Default)]
pub struct Clips {
    pub n_reads: u64,
    pub total_aligned: i64,
    pub total_soft: i64,
    pub total_hard: i64,
    pub clip_frac_sum: f64,
    pub len_count: u64,
    pub len_total: f64,
    pub len_min: Option<i64>,
    pub len_max: Option<i64>,
    pub len_values: FxHashMap<i64, u64>,
    pub len_order: Vec<i64>,
}

impl Clips {
    fn add(&mut self, b: &Breakdown) {
        self.n_reads += 1;
        self.total_aligned += b.query_aligned;
        self.total_soft += b.soft_clipped;
        self.total_hard += b.hard_clipped;
        let l = b.full_read_length;
        self.len_count += 1;
        self.len_total += l as f64;
        if self.len_min.is_none_or(|m| l < m) {
            self.len_min = Some(l);
        }
        if self.len_max.is_none_or(|m| l > m) {
            self.len_max = Some(l);
        }
        let e = self.len_values.entry(l).or_insert(0);
        if *e == 0 {
            self.len_order.push(l);
        }
        *e += 1;
        if l != 0 {
            self.clip_frac_sum += b.clipped as f64 / l as f64;
        }
    }
}

/// what's kept per read for mapped-read-stats' per-read output
pub struct ReadRow {
    pub read_id: String,
    pub pid: f64,
    pub gc_pid: f64,
    pub query_aligned: i64,
    pub full_read_length: i64,
    pub soft: i64,
    pub hard: i64,
    pub clipped_frac: f64,
    pub ref_start: i64,
    pub ref_end: i64,
}

/// whether a mapped alignment's gap-compressed percent identity is at least min_pid (false
/// without an NM tag or CIGAR, as it can't be judged)
fn passes_min_pid(r: &Record, min_pid: f64) -> bool {
    match r.nm {
        Some(nm) if !r.cigar.is_empty() => pids(&breakdown(&r.cigar), nm).1 >= min_pid,
        _ => false,
    }
}

/// for a read that's counted (mapped, primary unless include_non_primary, with NM): its
/// breakdown, its pids, whether it's non-primary
#[allow(clippy::type_complexity)]
fn counted(bam: &Bam, r: &Record, include_non_primary: bool) -> io::Result<Option<(Breakdown, (f64, f64), bool)>> {
    if r.flag & 0x4 != 0 {
        return Ok(None);
    }
    let non_primary = r.flag & 0x100 != 0 || r.flag & 0x800 != 0;
    if !include_non_primary && non_primary {
        return Ok(None);
    }
    let nm = match r.nm {
        Some(nm) => nm,
        None => return Ok(None),
    };
    if r.cigar.is_empty() || r.ref_id < 0 || r.ref_id as usize >= bam.ref_names.len() {
        // pysam's path errors on these
        return Err(bad("mapped read without a CIGAR or reference"));
    }
    let b = breakdown(&r.cigar);
    let p = pids(&b, nm);
    Ok(Some((b, p, non_primary)))
}

/// the reads counted: (record, its breakdown, its pids, whether it's non-primary)
fn each_counted(
    path: &str,
    include_non_primary: bool,
    mut f: impl FnMut(&Bam, &Record, &Breakdown, (f64, f64), bool) -> io::Result<()>,
) -> io::Result<Bam> {
    let mut bam = Bam::open(path)?;
    while let Some(r) = bam.next()? {
        if let Some((b, p, non_primary)) = counted(&bam, &r, include_non_primary)? {
            f(&bam, &r, &b, p, non_primary)?;
        }
    }
    Ok(bam)
}

/// get_mapped_reads_pids(): (rows grouped by reference in first-seen order, or None),
/// gap-aware stats, gap-compressed stats, clip aggregate
#[allow(clippy::type_complexity)]
pub fn mapped_reads_pids(
    path: &str,
    include_non_primary: bool,
    store_rows: bool,
) -> io::Result<(Option<Vec<(String, Vec<ReadRow>)>>, Acc, Acc, Clips)> {
    let (mut pid, mut gc, mut clips) = (Acc::default(), Acc::default(), Clips::default());
    let mut groups: Vec<(String, Vec<ReadRow>)> = Vec::new();
    let mut group_of: FxHashMap<String, usize> = FxHashMap::default();
    each_counted(path, include_non_primary, |bam, r, b, (p, g), non_primary| {
        pid.add(p);
        gc.add(g);
        if !non_primary {
            clips.add(b);
        }
        if store_rows {
            if !r.name.is_ascii() {
                return Err(bad("non-ASCII read name"));
            }
            let mut read_id = String::from_utf8(r.name.clone()).unwrap();
            if r.flag & 0x1 != 0 {
                if r.flag & 0x40 != 0 {
                    read_id.push_str("/1");
                } else if r.flag & 0x80 != 0 {
                    read_id.push_str("/2");
                }
            }
            let clipped_frac = if b.full_read_length != 0 {
                b.clipped as f64 / b.full_read_length as f64
            } else {
                0.0
            };
            // pysam's reference_end (htslib's bam_endpos): a span of 0 counts as 1
            let rlen: i64 = r
                .cigar
                .iter()
                .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
                .map(|&c| (c >> 4) as i64)
                .sum();
            let row = ReadRow {
                read_id,
                pid: p,
                gc_pid: g,
                query_aligned: b.query_aligned,
                full_read_length: b.full_read_length,
                soft: b.soft_clipped,
                hard: b.hard_clipped,
                clipped_frac,
                ref_start: r.pos as i64 + 1,
                ref_end: r.pos as i64 + if rlen == 0 { 1 } else { rlen },
            };
            // (grouped by name, as Python's dict is, which matters if a header repeats one)
            let name = &bam.ref_names[r.ref_id as usize];
            let gi = match group_of.get(name) {
                Some(&gi) => gi,
                None => {
                    groups.push((name.clone(), Vec::new()));
                    group_of.insert(name.clone(), groups.len() - 1);
                    groups.len() - 1
                }
            };
            groups[gi].1.push(row);
        }
        Ok(())
    })?;
    Ok((store_rows.then_some(groups), pid, gc, clips))
}

/// get_per_contig_pid_stats(): per reference (first-seen order), stats of one pid flavor
pub fn per_contig_pid_stats(
    path: &str,
    include_non_primary: bool,
    gap_compressed: bool,
) -> io::Result<Vec<(String, Acc)>> {
    let mut contig_pids = ContigPids::default();
    each_counted(path, include_non_primary, |bam, r, _, (ga, gc), _| {
        contig_pids.add(&bam.ref_names[r.ref_id as usize], if gap_compressed { gc } else { ga });
        Ok(())
    })?;
    Ok(contig_pids.accs)
}

/// the per-contig percent-ID stats of cov-stats, gathered a read at a time
#[derive(Default)]
struct ContigPids {
    // keyed by name (as Python's dict is, which matters if a header repeats one)
    accs: Vec<(String, Acc)>,
    idx: FxHashMap<String, usize>,
}

impl ContigPids {
    fn add(&mut self, name: &str, v: f64) {
        let i = match self.idx.get(name) {
            Some(&i) => i,
            None => {
                self.accs.push((name.to_string(), Acc::default()));
                self.idx.insert(name.to_string(), self.accs.len() - 1);
                self.accs.len() - 1
            }
        };
        self.accs[i].1.add(v);
    }
}

/// a contig's coverage: (total coverage count, bases detected at all, bases detected at 10x,
/// depth -> bases)
pub type CovTotals = (u128, u128, u128, FxHashMap<i64, u128>);

/// what the error says when a bam's alignments aren't grouped by reference (the caller sorts
/// the bam and comes back)
pub const NOT_GROUPED: &str = "alignments not grouped by reference";

/// Per-base depth of each wanted reference, from the alignments as they come: one reference
/// at a time is held as an array of depth changes (+1 where an alignment starts, -1 where it
/// ends), then summed along into depths and totalled when the next reference starts. So the
/// alignments of a reference have to be together (as in a sorted bam), in any order.
struct Coverage {
    slot: Vec<Option<usize>>, // per bam reference: its place in `wanted`
    seen: Vec<bool>,
    cur: usize,
    diff: Vec<i32>,
    small: Vec<u64>, // bases at each depth under SMALL, of the reference being totalled
    out: Vec<Option<CovTotals>>,
}

impl Coverage {
    const SMALL: usize = 1 << 16;

    fn new(bam: &Bam, wanted: &[String]) -> Self {
        let idx: FxHashMap<&str, usize> = wanted.iter().enumerate().map(|(i, w)| (w.as_str(), i)).collect();
        Coverage {
            slot: bam.ref_names.iter().map(|n| idx.get(n.as_str()).copied()).collect(),
            seen: vec![false; bam.ref_names.len()],
            cur: usize::MAX,
            diff: Vec::new(),
            small: vec![0; Self::SMALL],
            out: vec![None; wanted.len()],
        }
    }

    /// an alignment over [start, end) of a reference
    fn add(&mut self, bam: &Bam, ref_id: usize, start: usize, end: usize) -> io::Result<()> {
        if self.slot[ref_id].is_none() {
            return Ok(());
        }
        let len = bam.ref_lens[ref_id];
        if ref_id != self.cur {
            self.close(bam);
            if self.seen[ref_id] {
                return Err(bad(NOT_GROUPED));
            }
            self.seen[ref_id] = true;
            self.cur = ref_id;
            self.diff.clear();
            self.diff.resize(len + 1, 0);
        }
        let (start, end) = (start.min(len), end.min(len));
        if start < end {
            self.diff[start] += 1;
            self.diff[end] -= 1;
        }
        Ok(())
    }

    /// totals the reference held, if one is
    fn close(&mut self, bam: &Bam) {
        if self.cur == usize::MAX {
            return;
        }
        let len = bam.ref_lens[self.cur];
        let (mut depth, mut total, mut top) = (0i64, 0u64, 0usize);
        let mut big: FxHashMap<i64, u128> = FxHashMap::default();
        for &d in &self.diff[..len] {
            depth += d as i64;
            total += depth as u64;
            if (depth as usize) < Self::SMALL {
                self.small[depth as usize] += 1;
                top = top.max(depth as usize);
            } else {
                *big.entry(depth).or_insert(0) += 1;
            }
        }
        let under_10: u64 = self.small[..10].iter().sum();
        let (detected, detected_10x) = (len as u64 - self.small[0], len as u64 - under_10);
        let mut hist = big;
        for (depth, n) in self.small[..=top].iter_mut().enumerate() {
            if *n > 0 {
                hist.insert(depth as i64, *n as u128);
                *n = 0;
            }
        }
        self.merge(self.cur, total as u128, detected as u128, detected_10x as u128, hist);
        self.cur = usize::MAX;
    }

    fn merge(&mut self, ref_id: usize, total: u128, detected: u128, detected_10x: u128, hist: FxHashMap<i64, u128>) {
        let e = self.out[self.slot[ref_id].unwrap()].get_or_insert_with(|| (0, 0, 0, FxHashMap::default()));
        e.0 += total;
        e.1 += detected;
        e.2 += detected_10x;
        for (depth, n) in hist {
            *e.3.entry(depth).or_insert(0) += n;
        }
    }

    /// the totals, with the references that had no alignments in as all depth 0
    fn finish(mut self, bam: &Bam) -> Vec<Option<CovTotals>> {
        self.close(bam);
        for ref_id in 0..self.slot.len() {
            let len = bam.ref_lens[ref_id];
            if self.slot[ref_id].is_some() && !self.seen[ref_id] && len > 0 {
                self.merge(ref_id, 0, 0, 0, [(0i64, len as u128)].into_iter().collect());
            }
        }
        self.out
    }
}

/// cov-stats' one pass over a bam: per wanted contig its coverage totals (None if it isn't
/// in the bam), and, if asked for, the per-contig percent-ID stats (as per_contig_pid_stats()
/// gives them).
///
/// Coverage is counted as `mosdepth -x` does: every alignment not flagged with any of
/// `exclude` adds 1 from its start to its end on the reference, whatever is between
/// (deletions and skips included), and mates that overlap both count.
///
/// With `min_pid` above 0, alignments under that gap-compressed percent identity (and mapped
/// ones without an NM tag or CIGAR, which can't be judged) add to neither.
#[allow(clippy::type_complexity)]
pub fn coverage_and_pids(
    path: &str,
    wanted: &[String],
    exclude: u16,
    pid: Option<(bool, bool)>, // (include_non_primary, gap_compressed)
    min_pid: f64,
    progress: Option<&std::sync::Arc<crate::reads::Progress>>,
) -> io::Result<(Vec<Option<CovTotals>>, Option<Vec<(String, Acc)>>)> {
    let mut bam = Bam::open(path)?;
    let mut cov = Coverage::new(&bam, wanted);
    let mut contig_pids = ContigPids::default();
    let mut n = 0u64;
    while let Some(r) = bam.next()? {
        n += 1;
        if n.is_multiple_of(65536)
            && let Some(p) = progress
        {
            p.reads.store(n, std::sync::atomic::Ordering::Relaxed);
        }
        if min_pid > 0.0 && r.flag & 0x4 == 0 && !passes_min_pid(&r, min_pid) {
            continue;
        }
        if r.flag & exclude == 0 && r.ref_id >= 0 && r.pos >= 0 {
            // htslib's bam_endpos: a span of 0 counts as 1
            let rlen: usize = r
                .cigar
                .iter()
                .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
                .map(|&c| (c >> 4) as usize)
                .sum();
            let start = r.pos as usize;
            cov.add(&bam, r.ref_id as usize, start, start + rlen.max(1))?;
        }
        if let Some((include_non_primary, gap_compressed)) = pid
            && let Some((_, (ga, gc), _)) = counted(&bam, &r, include_non_primary)?
        {
            contig_pids.add(&bam.ref_names[r.ref_id as usize], if gap_compressed { gc } else { ga });
        }
    }
    if let Some(p) = progress {
        p.reads.store(n, std::sync::atomic::Ordering::Relaxed);
    }
    Ok((cov.finish(&bam), pid.map(|_| contig_pids.accs)))
}
