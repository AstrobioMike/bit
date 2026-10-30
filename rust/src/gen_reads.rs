//! The per-read loops of `bit gen-reads` (gen_one_file_paired and gen_one_file_single in
//! bit/modules/gen_reads.py), followed step for step, drawing from the same Mersenne Twister
//! state as the Python random.Random they'd otherwise use (its state is passed in and back
//! out via getstate()/setstate()), with the same derivations of random(), randint() and
//! gauss(). So with the same seed, reads, headers and coordinates are identical. The Python
//! side checks this generator against random.Random at runtime before using it.

use std::f64::consts::TAU;
use std::fmt::Write;

const N: usize = 624;
const M: usize = 397;

/// CPython's random.Random core (Modules/_randommodule.c)
pub struct Mt {
    pub state: [u32; N],
    pub index: usize,
    pub gauss_next: Option<f64>,
}

impl Mt {
    pub fn genrand_u32(&mut self) -> u32 {
        const MAG01: [u32; 2] = [0, 0x9908_b0df];
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        let mt = &mut self.state;
        if self.index >= N {
            let mut kk = 0;
            while kk < N - M {
                let y = (mt[kk] & UPPER) | (mt[kk + 1] & LOWER);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
                kk += 1;
            }
            while kk < N - 1 {
                let y = (mt[kk] & UPPER) | (mt[kk + 1] & LOWER);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
                kk += 1;
            }
            let y = (mt[N - 1] & UPPER) | (mt[0] & LOWER);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
            self.index = 0;
        }
        let mut y = mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// random()
    pub fn random(&mut self) -> f64 {
        let a = (self.genrand_u32() >> 5) as f64;
        let b = (self.genrand_u32() >> 6) as f64;
        (a * 67108864.0 + b) * (1.0 / 9007199254740992.0)
    }

    /// getrandbits(k), for k <= 64 (words least significant first, as CPython fills them)
    fn getrandbits(&mut self, k: u32) -> u64 {
        if k == 0 {
            return 0;
        }
        if k <= 32 {
            return (self.genrand_u32() >> (32 - k)) as u64;
        }
        let lo = self.genrand_u32() as u64;
        let hi = (self.genrand_u32() >> (64 - k)) as u64;
        lo | (hi << 32)
    }

    /// _randbelow_with_getrandbits(n)
    fn randbelow(&mut self, n: u64) -> u64 {
        let k = 64 - n.leading_zeros();
        let mut r = self.getrandbits(k);
        while r >= n {
            r = self.getrandbits(k);
        }
        r
    }

    /// randint(a, b) (randrange(a, b + 1)); None where Python raises (an empty range)
    pub fn randint(&mut self, a: i64, b: i64) -> Option<i64> {
        let width = b + 1 - a;
        if width <= 0 {
            return None;
        }
        Some(a + self.randbelow(width as u64) as i64)
    }

    /// gauss(mu, sigma)
    pub fn gauss(&mut self, mu: f64, sigma: f64) -> f64 {
        let z = match self.gauss_next.take() {
            Some(z) => z,
            None => {
                let x2pi = self.random() * TAU;
                let g2rad = (-2.0 * (1.0 - self.random()).ln()).sqrt();
                self.gauss_next = Some(x2pi.sin() * g2rad);
                x2pi.cos() * g2rad
            }
        };
        mu + z * sigma
    }
}

/// Python's round() of a float to an int
fn pyround(x: f64) -> i64 {
    x.round_ties_even() as i64
}

/// reverse_complement() (ACGTacgt complemented, others kept, reversed)
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
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
            o => o,
        })
        .collect()
}

/// extract_subsequence(): (subsequence, start); None where Python raises
fn extract(
    rng: &mut Mt,
    seq: &[u8],
    subseq_len: usize,
    circular: bool,
    include_ns: bool,
) -> Option<(Vec<u8>, i64)> {
    let l = seq.len();
    let mut last: Option<(Vec<u8>, i64)> = None;
    for _ in 0..100 {
        let (sub, start) = if circular {
            let start = rng.randint(0, l as i64 - 1)? as usize;
            let end = start + subseq_len;
            if end <= l {
                (seq[start..end].to_vec(), start)
            } else {
                let mut s = seq[start..].to_vec();
                s.extend_from_slice(&seq[..(end - l).min(l)]);
                (s, start)
            }
        } else if subseq_len >= l {
            let start = rng.randint(0, l as i64 - 1)? as usize;
            (seq[start..].to_vec(), start)
        } else {
            let start = rng.randint(0, (l - subseq_len) as i64)? as usize;
            (seq[start..start + subseq_len].to_vec(), start)
        };
        let has_n = sub.iter().any(|&b| b == b'N' || b == b'n');
        if include_ns || !has_n {
            return Some((sub, start as i64));
        }
        last = Some((sub, start as i64));
    }
    last
}

pub struct Params<'a> {
    pub paired: bool,
    pub long: bool,
    pub read_length: i64,
    pub mu: f64,     // fragment size (paired) or read length (long), for gauss
    pub sigma: f64,  // its spread
    pub min_len: i64, // min/max fragment (paired) or read (long) length
    pub max_len: i64,
    pub circular: bool,
    pub include_ns: bool,
    pub per_read_tsv: bool,
    pub seq_id: &'a str,
    pub source_fasta: &'a str,
}

#[derive(Default)]
pub struct Chunk {
    pub r1: String,
    pub r2: String,
    pub tsv: String,
    pub spans: Vec<(i64, i64)>, // merged, sorted
    pub realized_bases: u64,
}

fn header(out: &mut String, p: &Params, id: &str, start: i64, end: i64, strand: char, wrapped: bool) {
    out.push('@');
    out.push_str(id);
    if !p.per_read_tsv {
        let _ = write!(out, " contig={};start={};end={};strand={}", p.seq_id, start, end, strand);
        if wrapped {
            out.push_str(";wrapped=true");
        }
    }
    out.push('\n');
}

#[allow(clippy::too_many_arguments)]
fn record(out: &mut String, p: &Params, id: &str, start: i64, end: i64, strand: char, wrapped: bool, seq: &[u8]) {
    header(out, p, id, start, end, strand, wrapped);
    // (sequences are ASCII here, checked by the caller)
    out.push_str(std::str::from_utf8(seq).unwrap());
    out.push_str("\n+\n");
    for _ in 0..seq.len() {
        out.push('I');
    }
    out.push('\n');
}

fn tsv_row(out: &mut String, p: &Params, id: &str, start: i64, end: i64, strand: char, wrapped: bool) {
    let _ = writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        id,
        p.source_fasta,
        p.seq_id,
        start,
        end,
        strand,
        if wrapped { "true" } else { "false" }
    );
}

/// _record_span(): a span, split in two if it wraps the origin
fn add_span(spans: &mut Vec<(i64, i64)>, start: i64, end: i64, l: i64) {
    if end <= l {
        spans.push((start, end));
    } else {
        spans.push((start, l));
        spans.push((0, end - l));
    }
}

fn merge(mut spans: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    spans.retain(|&(s, e)| s < e);
    spans.sort_unstable();
    let mut out: Vec<(i64, i64)> = Vec::with_capacity(spans.len());
    for (s, e) in spans {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// `n` reads (or pairs) from one contig's sequence `seq` (already uppercased for paired, as
/// Python does), numbered from read_count + 1. None where Python would raise (the caller then
/// runs Python's loop for these reads from the same state).
pub fn gen_chunk(rng: &mut Mt, seq: &[u8], n: u64, read_count: u64, p: &Params) -> Option<Chunk> {
    let l = seq.len() as i64;
    let mut c = Chunk::default();
    let mut spans = Vec::new();
    for i in 0..n {
        let id = format!("r{}", read_count + i + 1);
        if p.paired {
            let frag_len = pyround(rng.gauss(p.mu, p.sigma)).max(p.min_len).min(p.max_len).min(l);
            if frag_len < 0 {
                return None;
            }
            let (mut fragment, start) = extract(rng, seq, frag_len as usize, p.circular, p.include_ns)?;
            let strand = if rng.random() < 0.5 {
                '+'
            } else {
                fragment = revcomp(&fragment);
                '-'
            };
            let fl = fragment.len();
            let rl = p.read_length;
            if rl < 0 {
                return None;
            }
            // fragment[:rl] and fragment[-rl:] (with Python's slicing: [-0:] is everything)
            let fwd_len = if rl >= 0 { (rl as usize).min(fl) } else { fl.saturating_sub((-rl) as usize) };
            let fwd = &fragment[..fwd_len];
            let rev_from = if rl == 0 {
                0
            } else if rl > 0 {
                fl.saturating_sub(rl as usize)
            } else {
                ((-rl) as usize).min(fl)
            };
            let rev: Vec<u8> = fragment[rev_from..]
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    o => o,
                })
                .collect();
            // compute_mate_coords()
            let (mut r1s, mut r1e) = (start, start + fwd.len() as i64);
            let (mut r2e, mut r2s) = (start + fl as i64, start + fl as i64 - rev.len() as i64);
            let wrapped = start + fl as i64 > l;
            if strand == '-' {
                std::mem::swap(&mut r1s, &mut r2s);
                std::mem::swap(&mut r1e, &mut r2e);
            }
            let s2 = if strand == '+' { '-' } else { '+' };
            let (id1, id2) = (format!("{id}/1"), format!("{id}/2"));
            record(&mut c.r1, p, &id1, r1s, r1e, strand, wrapped, fwd);
            record(&mut c.r2, p, &id2, r2s, r2e, s2, wrapped, &rev);
            if p.per_read_tsv {
                tsv_row(&mut c.tsv, p, &id1, r1s, r1e, strand, wrapped);
                tsv_row(&mut c.tsv, p, &id2, r2s, r2e, s2, wrapped);
            }
            add_span(&mut spans, r1s, r1e, l);
            add_span(&mut spans, r2s, r2e, l);
            c.realized_bases += (fwd.len() + rev.len()) as u64;
        } else {
            let read_len = if p.long {
                pyround(rng.gauss(p.mu, p.sigma)).max(p.min_len).min(p.max_len).min(l)
            } else {
                p.read_length.min(l)
            };
            if read_len < 0 {
                return None;
            }
            let (mut read, start) = extract(rng, seq, read_len as usize, p.circular, p.include_ns)?;
            let read_end = start + read.len() as i64;
            let wrapped = read_end > l;
            let strand = if rng.random() < 0.5 {
                '+'
            } else {
                read = revcomp(&read);
                '-'
            };
            record(&mut c.r1, p, &id, start, read_end, strand, wrapped, &read);
            if p.per_read_tsv {
                tsv_row(&mut c.tsv, p, &id, start, read_end, strand, wrapped);
            }
            add_span(&mut spans, start, read_end, l);
            c.realized_bases += read.len() as u64;
        }
        if spans.len() > 1 << 20 {
            spans = merge(spans);
        }
    }
    c.spans = merge(spans);
    Some(c)
}
