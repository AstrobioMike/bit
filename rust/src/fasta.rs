//! FASTA passes for `bit fasta` subcommands and `bit summarize-assembly`, each following the
//! parser its Python version uses, so results are identical:
//!   - Biopython's SeqIO "fasta" (FastaIterator): count, calc-gc, filter-by-length
//!   - bit's own read_fasta() (seqs.py): extract-by-primers
//!   - plain line handling: remove-wraps
//!   - summarize-assembly: as pyfastx read assemblies (see assembly_summary)
//!
//! Input is plain text or gzipped (told by its first two bytes), its lines ending in "\n", "\r\n" or
//! a lone "\r" (as Python's text mode reads them). Header lines (">...") may hold any UTF-8 text; any
//! other line that isn't ASCII is an error (it's no sequence file).

use crate::edlib;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

pub fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// ASCII characters Python's str.split()/strip() count as whitespace
fn py_ws(b: u8) -> bool {
    matches!(b, 0x09..=0x0d | 0x1c..=0x1f | 0x20)
}

fn rstrip_py(s: &[u8]) -> &[u8] {
    let mut e = s.len();
    while e > 0 && py_ws(s[e - 1]) {
        e -= 1;
    }
    &s[..e]
}

fn strip_py(s: &[u8]) -> &[u8] {
    let s = rstrip_py(s);
    let mut b = 0;
    while b < s.len() && py_ws(s[b]) {
        b += 1;
    }
    &s[b..]
}

/// title.split(None, 1)[0], or "" if there's no word
fn first_word(title: &[u8]) -> &[u8] {
    let mut b = 0;
    while b < title.len() && py_ws(title[b]) {
        b += 1;
    }
    let mut e = b;
    while e < title.len() && !py_ws(title[e]) {
        e += 1;
    }
    &title[b..e]
}

/// a file's bytes, gunzipped if it's gzipped
pub fn open_text(path: &str) -> io::Result<Box<dyn Read + Send>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    let start = io::Cursor::new(magic[..n].to_vec()).chain(file);
    Ok(if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, start)))
    } else {
        Box::new(start)
    })
}

/// lines as Python's text mode gives them ("\r\n" and a lone "\r" become "\n"), each with its
/// "\n" if it had one
pub struct Lines {
    rdr: BufReader<Box<dyn Read + Send>>,
    buf: Vec<u8>,
    // the further lines of a chunk holding lone "\r"s
    queued: VecDeque<Vec<u8>>,
}

impl Lines {
    pub fn open(path: &str) -> io::Result<Self> {
        Ok(Lines { rdr: BufReader::with_capacity(1 << 20, open_text(path)?), buf: Vec::new(), queued: VecDeque::new() })
    }

    pub fn next_line(&mut self) -> io::Result<Option<&[u8]>> {
        if let Some(line) = self.queued.pop_front() {
            check_text(&line)?;
            self.buf = line;
            return Ok(Some(&self.buf));
        }
        self.buf.clear();
        let n = self.rdr.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(None);
        }
        if self.buf.ends_with(b"\r\n") {
            let l = self.buf.len();
            self.buf.truncate(l - 2);
            self.buf.push(b'\n');
        }
        if memchr::memchr(b'\r', &self.buf).is_some() {
            // old Mac line endings: each "\r" ends a line
            let chunk = std::mem::take(&mut self.buf);
            let pieces: Vec<&[u8]> = chunk.split(|&b| b == b'\r').collect();
            let last = pieces.len() - 1;
            for (i, piece) in pieces.into_iter().enumerate() {
                if i < last {
                    let mut line = piece.to_vec();
                    line.push(b'\n');
                    self.queued.push_back(line);
                } else if !piece.is_empty() {
                    self.queued.push_back(piece.to_vec());
                }
            }
            return self.next_line();
        }
        check_text(&self.buf)?;
        Ok(Some(&self.buf))
    }
}

/// a header line may hold any UTF-8 text (e.g. ">seq1 β-proteobacterium"); other lines need to be ASCII
fn check_text(line: &[u8]) -> io::Result<()> {
    if line.is_ascii() {
        return Ok(());
    }
    if line.iter().find(|&&b| !py_ws(b)) == Some(&b'>') {
        return match std::str::from_utf8(line) {
            Ok(_) => Ok(()),
            Err(_) => Err(bad("has a header line that isn't UTF-8 text")),
        };
    }
    Err(bad("holds a line that isn't ASCII text outside the headers, so it isn't a plain-text sequence file"))
}

fn line_body(l: &[u8]) -> &[u8] {
    l.strip_suffix(b"\n").unwrap_or(l)
}

/// Biopython's FastaIterator: (title, sequence) for each record
pub struct BioFasta {
    lines: Lines,
    pending: Option<Vec<u8>>,
}

impl BioFasta {
    pub fn open(path: &str) -> io::Result<Self> {
        let mut lines = Lines::open(path)?;
        let pending = match lines.next_line()? {
            None => None,
            Some(l) if l.first() == Some(&b'>') => Some(l.to_vec()),
            // Biopython raises on text before the first record
            Some(_) => return Err(bad("doesn't start with '>'")),
        };
        Ok(BioFasta { lines, pending })
    }

    /// (title: the header line minus '>', right-stripped; sequence: the following lines
    /// with spaces, tabs, carriage returns and newlines removed)
    pub fn next_record(&mut self) -> io::Result<Option<(Vec<u8>, Vec<u8>)>> {
        let header = match self.pending.take() {
            Some(h) => h,
            None => return Ok(None),
        };
        let title = rstrip_py(&header[1..]).to_vec();
        let mut seq = Vec::new();
        while let Some(l) = self.lines.next_line()? {
            if l[0] == b'>' {
                self.pending = Some(l.to_vec());
                break;
            }
            seq.extend(l.iter().filter(|&&b| !matches!(b, b' ' | b'\t' | b'\r' | b'\n')));
        }
        Ok(Some((title, seq)))
    }
}

/// per record: (id/name, length, count of 'G' and 'C')
pub fn records_summary(path: &str) -> io::Result<Vec<(String, usize, usize)>> {
    let mut f = BioFasta::open(path)?;
    let mut out = Vec::new();
    while let Some((title, seq)) = f.next_record()? {
        let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
        out.push((String::from_utf8(first_word(&title).to_vec()).unwrap(), seq.len(), gc));
    }
    Ok(out)
}

/// per record: (name, length, count of 'G' and 'C', G+C count of each window), windows
/// starting every `step` bases while a whole window fits (as range(0, len - window + 1, step))
#[allow(clippy::type_complexity)]
pub fn gc_windows(path: &str, window: usize, step: usize) -> io::Result<Vec<(String, usize, usize, Vec<u64>)>> {
    let mut f = BioFasta::open(path)?;
    let mut out = Vec::new();
    while let Some((title, seq)) = f.next_record()? {
        let mut prefix = Vec::with_capacity(seq.len() + 1);
        prefix.push(0u64);
        let mut c = 0u64;
        for &b in &seq {
            if b == b'G' || b == b'C' {
                c += 1;
            }
            prefix.push(c);
        }
        let windows = if seq.len() >= window {
            (0..=seq.len() - window)
                .step_by(step)
                .map(|i| prefix[i + window] - prefix[i])
                .collect()
        } else {
            Vec::new()
        };
        out.push((
            String::from_utf8(first_word(&title).to_vec()).unwrap(),
            seq.len(),
            c as usize,
            windows,
        ));
    }
    Ok(out)
}

/// filter_fasta_by_length(): writes kept records as ">{title}\n{seq}\n"; returns (seqs,
/// seqs kept, bases, bases kept). Nothing is written if the input calls for the Python path.
pub fn filter_by_length(
    path: &str,
    out_path: &str,
    min_length: Option<i64>,
    max_length: i64,
) -> io::Result<(u64, u64, u64, u64)> {
    // a first pass checks the whole input can be read here before anything is written
    {
        let mut f = BioFasta::open(path)?;
        while f.next_record()?.is_some() {}
    }
    let mut f = BioFasta::open(path)?;
    let mut out = BufWriter::with_capacity(1 << 20, File::create(out_path)?);
    let (mut n, mut kept, mut bases, mut bases_kept) = (0u64, 0u64, 0u64, 0u64);
    while let Some((title, seq)) = f.next_record()? {
        let l = seq.len() as i64;
        n += 1;
        bases += l as u64;
        if min_length.is_none_or(|m| l >= m) && l <= max_length {
            kept += 1;
            bases_kept += l as u64;
            out.write_all(b">")?;
            out.write_all(&title)?;
            out.write_all(b"\n")?;
            out.write_all(&seq)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok((n, kept, bases, bases_kept))
}

/// checks a whole file is readable here: ASCII outside its (UTF-8) headers
pub fn check_plain(path: &str) -> io::Result<()> {
    let mut lines = Lines::open(path)?;
    while lines.next_line()?.is_some() {}
    Ok(())
}

/// remove_wraps(): header lines as they are, other lines joined up to the next header
pub fn remove_wraps(path: &str, out: &mut impl Write) -> io::Result<()> {
    let mut lines = Lines::open(path)?;
    let mut out = BufWriter::with_capacity(1 << 20, out);
    let mut n: &[u8] = b"";
    while let Some(l) = lines.next_line()? {
        if l.first() == Some(&b'>') {
            out.write_all(n)?;
            out.write_all(l)?;
            n = b"";
        } else {
            out.write_all(line_body(l))?;
            n = b"\n";
        }
    }
    out.write_all(n)?;
    out.flush()
}

/// bit's read_fasta(): (header, sequence) records. Lines are stripped; a stripped line
/// starting with '>' is a header (minus the '>'); a record is only kept if its header is
/// non-empty (with the sequence lines before the first header going nowhere)
pub fn read_fasta(path: &str) -> io::Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let mut lines = Lines::open(path)?;
    let mut out = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut seq = Vec::new();
    while let Some(l) = lines.next_line()? {
        let l = strip_py(l);
        if l.first() == Some(&b'>') {
            if let Some(h) = header.take()
                && !h.is_empty()
            {
                out.push((h, std::mem::take(&mut seq)));
            }
            header = Some(l[1..].to_vec());
            seq.clear();
        } else {
            seq.extend_from_slice(l);
        }
    }
    if let Some(h) = header
        && !h.is_empty()
    {
        out.push((h, seq));
    }
    Ok(out)
}

/// _find_primer_sites() of extract_seqs.py: every site of `primer` within `masked` (non-ACGT
/// already masked) allowing up to `k` edits, found by searching, masking the sites found, and
/// searching again until none are left within the limit, as (start, end exclusive)
fn primer_sites(primer: &[u8], masked: &[u8], k: usize, eq: &[(u8, u8)]) -> io::Result<Vec<(usize, usize)>> {
    let mut sites = Vec::new();
    if masked.is_empty() {
        return Ok(sites);
    }
    let mut target = masked.to_vec();
    loop {
        let locs = edlib::all_locations(primer, &target, k, eq).map_err(|_| bad("edlib error"))?;
        if locs.is_empty() {
            break;
        }
        for &(start, end_incl) in &locs {
            let end = end_incl + 1;
            sites.push((start, end));
            let stop = end.min(target.len());
            for b in &mut target[start.min(stop)..stop] {
                *b = b'X';
            }
        }
    }
    Ok(sites)
}

/// each primer's sites in one sequence (uppercased, non-ACGT positions masked first; see primer_sites)
pub fn seq_primer_sites(seq: &[u8], primers: &[Vec<u8>], k: usize, eq: &[(u8, u8)]) -> io::Result<Vec<Vec<(usize, usize)>>> {
    if primers.iter().any(|p| p.len() <= k) {
        return Err(bad("primers need to be longer than the allowed differences"));
    }
    let masked: Vec<u8> = seq
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .map(|b| if matches!(b, b'A' | b'C' | b'G' | b'T') { b } else { b'X' })
        .collect();
    primers.iter().map(|p| primer_sites(p, &masked, k, eq)).collect()
}

/// For extract-by-primers: each record of `path` (read_fasta), uppercased, with each primer's
/// sites (see primer_sites; non-ACGT target positions masked first). Returns (header,
/// uppercased sequence, [[(start, end exclusive), ...] per primer]) for records with at least
/// two sites in all, since amplicons need two. Primers must be longer than k (checked by the
/// caller), or the search could repeat forever.
#[allow(clippy::type_complexity)]
pub fn primer_hits(
    path: &str,
    primers: &[Vec<u8>],
    k: usize,
    eq: &[(u8, u8)],
) -> io::Result<Vec<(String, String, Vec<Vec<(usize, usize)>>)>> {
    if primers.iter().any(|p| p.len() <= k) {
        return Err(bad("primers need to be longer than the allowed differences"));
    }
    let records = read_fasta(path)?;
    let found: Vec<Option<(String, String, Vec<Vec<(usize, usize)>>)>> = records
        .into_par_iter()
        .map(|(h, s)| {
            let s = s.to_ascii_uppercase();
            let masked: Vec<u8> = s
                .iter()
                .map(|&b| if matches!(b, b'A' | b'C' | b'G' | b'T') { b } else { b'X' })
                .collect();
            // (the primers are searched in parallel too, which is what helps for a few long records)
            let sites: Vec<Vec<(usize, usize)>> = primers
                .par_iter()
                .map(|p| primer_sites(p, &masked, k, eq))
                .collect::<io::Result<_>>()?;
            if sites.iter().map(|l| l.len()).sum::<usize>() < 2 {
                return Ok(None);
            }
            Ok(Some((String::from_utf8(h).unwrap(), String::from_utf8(s).unwrap(), sites)))
        })
        .collect::<io::Result<_>>()?;
    Ok(found.into_iter().flatten().collect())
}

/// For summarize-assembly: (number of records, total length, count of each byte value in the
/// sequences, each record's length). A record is a header line ('>') and the lines after it, every
/// character of which counts except spaces and tabs (dropped, as Biopython drops them for `bit fasta
/// count` and calc-gc; blank lines are skipped, and a record with no sequence has length 0). Plain
/// or gzipped, any line endings. Text before the first header is an error.
pub fn assembly_summary(path: &str) -> io::Result<(u64, u64, Vec<u64>, Vec<u64>)> {
    let mut lines = Lines::open(path)?;
    let mut comp = vec![0u64; 256];
    let mut lengths: Vec<u64> = Vec::new();
    while let Some(line) = lines.next_line()? {
        let body = line_body(line);
        if body.first() == Some(&b'>') {
            lengths.push(0);
            continue;
        }
        let mut bases = body.iter().filter(|&&b| b != b' ' && b != b'\t').peekable();
        if bases.peek().is_none() {
            continue;
        }
        let Some(length) = lengths.last_mut() else {
            return Err(bad("isn't a fasta file (it doesn't start with a '>' header line)"));
        };
        for &b in bases {
            comp[b as usize] += 1;
            *length += 1;
        }
    }
    let size = lengths.iter().sum();
    Ok((lengths.len() as u64, size, comp, lengths))
}
