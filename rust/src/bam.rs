//! A small BAM reader (records in file order, as pysam's fetch(until_eof=True) gives them) and
//! the per-read passes of `bit mapped-read-stats` and `bit cov-stats` (see get_mapped_reads_pids
//! and get_per_contig_pid_stats in bit/modules/mapped_read_stats.py, which this follows
//! step for step). BGZF blocks are decompressed in parallel.
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

/// reads BGZF, decompressing a batch of blocks at a time in parallel
struct Bgzf {
    raw: BufReader<File>,
    out: VecDeque<u8>,
    done: bool,
}

impl Bgzf {
    const BATCH: usize = 256;

    fn open(path: &str) -> io::Result<Self> {
        Ok(Bgzf {
            raw: BufReader::with_capacity(1 << 20, File::open(path)?),
            out: VecDeque::new(),
            done: false,
        })
    }

    /// the next raw block's compressed deflate data and expected size, or None at the end
    fn next_block(&mut self) -> io::Result<Option<(Vec<u8>, usize, u32)>> {
        let mut head = [0u8; 12];
        let mut got = 0;
        while got < 12 {
            let n = self.raw.read(&mut head[got..])?;
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
        self.raw.read_exact(&mut extra)?;
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
        self.raw.read_exact(&mut cdata)?;
        let mut tail = [0u8; 8];
        self.raw.read_exact(&mut tail)?;
        let crc = u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]);
        let isize = u32::from_le_bytes([tail[4], tail[5], tail[6], tail[7]]) as usize;
        Ok(Some((cdata, isize, crc)))
    }

    fn fill(&mut self) -> io::Result<()> {
        let mut blocks = Vec::with_capacity(Self::BATCH);
        while blocks.len() < Self::BATCH {
            match self.next_block()? {
                Some(b) => blocks.push(b),
                None => {
                    self.done = true;
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
        for b in inflated {
            self.out.extend(b?);
        }
        Ok(())
    }
}

impl Read for Bgzf {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.out.is_empty() && !self.done {
            self.fill()?;
        }
        let n = buf.len().min(self.out.len());
        for (i, b) in self.out.drain(..n).enumerate() {
            buf[i] = b;
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
            read_i32(&mut inner)?;
        }
        Ok(Bam {
            inner,
            ref_names,
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

/// the reads counted (mapped, primary unless include_non_primary, with NM): (record, its
/// breakdown, its pids, whether it's non-primary)
fn each_counted(
    path: &str,
    include_non_primary: bool,
    mut f: impl FnMut(&Bam, &Record, &Breakdown, (f64, f64), bool) -> io::Result<()>,
) -> io::Result<Bam> {
    let mut bam = Bam::open(path)?;
    while let Some(r) = bam.next()? {
        if r.flag & 0x4 != 0 {
            continue;
        }
        let non_primary = r.flag & 0x100 != 0 || r.flag & 0x800 != 0;
        if !include_non_primary && non_primary {
            continue;
        }
        let nm = match r.nm {
            Some(nm) => nm,
            None => continue,
        };
        if r.cigar.is_empty() || r.ref_id < 0 || r.ref_id as usize >= bam.ref_names.len() {
            // pysam's path errors on these
            return Err(bad("mapped read without a CIGAR or reference"));
        }
        let b = breakdown(&r.cigar);
        let p = pids(&b, nm);
        f(&bam, &r, &b, p, non_primary)?;
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
    // keyed by name (as Python's dict is, which matters if a header repeats one)
    let mut accs: Vec<(String, Acc)> = Vec::new();
    let mut idx: FxHashMap<String, usize> = FxHashMap::default();
    each_counted(path, include_non_primary, |bam, r, _, (ga, gc), _| {
        let name = &bam.ref_names[r.ref_id as usize];
        let i = match idx.get(name) {
            Some(&i) => i,
            None => {
                accs.push((name.clone(), Acc::default()));
                idx.insert(name.clone(), accs.len() - 1);
                accs.len() - 1
            }
        };
        accs[i].1.add(if gap_compressed { gc } else { ga });
        Ok(())
    })?;
    Ok(accs)
}

/// the per-base.bed.gz pass of cov-stats (parse_bed_file): per contig wanted, (total coverage
/// count, bases detected at all, bases detected at 10x, depth -> bases). Lines for other
/// contigs are parsed (as Python does) but not kept. Anything Python would read differently
/// (non-UTF-8, a lone carriage return, a number int() takes but that isn't plain digits) is
/// an error, so the caller falls back to Python.
pub type BedStats = (u128, u128, u128, FxHashMap<i64, u128>);

pub fn per_base_bed(
    path: &str,
    wanted: &[String],
    progress: Option<&std::sync::Arc<crate::reads::Progress>>,
) -> io::Result<Vec<Option<BedStats>>> {
    use std::io::BufRead;
    let idx: FxHashMap<&[u8], usize> = wanted.iter().enumerate().map(|(i, w)| (w.as_bytes(), i)).collect();
    let mut out: Vec<Option<BedStats>> = vec![None; wanted.len()];

    let file = File::open(path)?;
    let gz = flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file));
    let mut rdr = BufReader::with_capacity(1 << 20, gz);

    // int() of a field: optional surrounding ASCII whitespace, optional sign, digits
    fn int(f: &[u8]) -> io::Result<i64> {
        let t = f.trim_ascii();
        let (neg, d) = match t.first() {
            Some(b'-') => (true, &t[1..]),
            Some(b'+') => (false, &t[1..]),
            _ => (false, t),
        };
        if d.is_empty() || d.len() > 18 || !d.iter().all(|c| c.is_ascii_digit()) {
            return Err(bad("unusual number"));
        }
        let v: i64 = std::str::from_utf8(d).unwrap().parse().map_err(|_| bad("unusual number"))?;
        Ok(if neg { -v } else { v })
    }

    let mut line = Vec::new();
    let mut n_lines = 0u64;
    loop {
        line.clear();
        let n = rdr.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        n_lines += 1;
        if n_lines.is_multiple_of(65536)
            && let Some(p) = progress
        {
            p.reads.store(n_lines, std::sync::atomic::Ordering::Relaxed);
        }
        // Python's text mode turns "\r\n" into "\n" and splits on a lone "\r"; only the
        // first is handled here
        let mut l: &[u8] = &line;
        if l.ends_with(b"\r\n") {
            l = &l[..l.len() - 2];
        } else if l.ends_with(b"\n") {
            l = &l[..l.len() - 1];
        }
        if memchr::memchr(b'\r', l).is_some() || !l.is_ascii() {
            return Err(bad("unusual line"));
        }
        let mut fields = l.splitn(5, |&c| c == b'\t');
        let header = fields.next().unwrap();
        let (a, b, c) = match (fields.next(), fields.next(), fields.next()) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => return Err(bad("short line")),
        };
        // (a 5th+ field only matters to Python as part of fields[3] if there's no 5th tab;
        // with one, fields[3] is just the 4th field, as here)
        let span = int(b)? - int(a)?;
        let n_reads = int(c)?;
        if let Some(&i) = idx.get(header) {
            if span < 0 || n_reads < 0 {
                return Err(bad("negative span or depth"));
            }
            let e = out[i].get_or_insert_with(|| (0, 0, 0, FxHashMap::default()));
            let (span_u, reads_u) = (span as u128, n_reads as u128);
            e.0 += reads_u * span_u;
            if n_reads > 0 {
                e.1 += span_u;
            }
            if n_reads >= 10 {
                e.2 += span_u;
            }
            *e.3.entry(n_reads).or_insert(0) += span_u;
        }
    }
    if let Some(p) = progress {
        p.reads.store(n_lines, std::sync::atomic::Ordering::Relaxed);
    }
    Ok(out)
}
