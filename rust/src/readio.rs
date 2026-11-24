//! Reading reads in and writing them back out, for commands that pass a read file (or a pair of
//! them) through and write out the reads kept, as they were in the input (names, comments,
//! qualities and all; `bit filter-foldbacks`, `bit match-reads`): batches of records read on
//! threads of their own (each mate's file on its own), records cut to a stretch, and gzipped
//! output compressed on several threads.

use crate::reads::{Progress, SeqReader};
use flate2::Compression;
use flate2::write::GzEncoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::Arc;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{Receiver, Sender, channel, sync_channel};

// output is gzipped this many bytes at a time, each its own gzip member (a file of several
// members is still one gzip file to every reader), so it can be done on several threads
const GZ_CHUNK: usize = 1 << 20;
const GZ_PENDING: usize = 8 * GZ_CHUNK;
// gz_members() cuts text this small, so even one batch's text spreads over many threads
const GZ_PIECE: usize = 256 << 10;
// outputs are compressed at level 4, as BBTools does by default (bbnorm's output, Nov 2026):
// level 1 was quicker but gave files up to twice the size; 6 (gzip's default) gave files 3-7%
// smaller than 4 for 1.2-1.5x the time
const GZ_LEVEL: u32 = 4;

pub fn bad(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// a read, or a pair: sequences (uppercased), and the records as in the files
#[derive(Default)]
pub struct Rec {
    pub seq1: Vec<u8>,
    pub raw1: Vec<u8>,
    pub seq2: Option<Vec<u8>>,
    pub raw2: Option<Vec<u8>>,
}

/// a record's name: its header's first word, less the marker and any /1 or /2
pub fn base_name(raw: &[u8]) -> &[u8] {
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

pub fn name_str(raw: &[u8]) -> String {
    let line = raw.get(1..).unwrap_or(&[]);
    let end = line
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(line.len());
    String::from_utf8_lossy(&line[..end]).into_owned()
}

/// batches of reads (or pairs, their names checked to match), read on a thread of their own, a
/// batch ahead. With pairs, read 2s are read on another thread still, so the two files are
/// decompressed side by side (decompressing gzip is the slowest part of reading). Spent batches
/// sent back on the returned sender are refilled in place, so the reads' buffers are reused
/// instead of freed and allocated again for every read (that churn, and freeing each batch on
/// the thread that hands out the next, held up both passes of normalize-reads, Nov 2026).
/// Also returns whether the input is fastq.
#[allow(clippy::type_complexity)]
pub fn spawn_reader(
    path1: String,
    path2: Option<String>,
    batch_size: usize,
    progress: Option<Arc<Progress>>,
) -> io::Result<(Receiver<io::Result<Vec<Rec>>>, bool, Sender<Vec<Rec>>)> {
    let mut r1 = SeqReader::open(&path1, progress.clone())?;
    let mates = match &path2 {
        Some(p) => Some(spawn_mate_reader(SeqReader::open(p, progress.clone())?, batch_size)),
        None => None,
    };
    let fastq = r1.is_fastq();
    let (tx, rx) = sync_channel(2);
    let (spent_tx, spent) = channel::<Vec<Rec>>();
    std::thread::spawn(move || {
        let result = (|| -> io::Result<()> {
            loop {
                let mut batch = spent.try_recv().unwrap_or_else(|_| Vec::with_capacity(batch_size));
                let mut n = 0;
                while n < batch_size {
                    if n == batch.len() {
                        batch.push(Rec::default());
                    }
                    let rec = &mut batch[n];
                    if !r1.next_record(&mut rec.seq1, &mut rec.raw1)? {
                        break;
                    }
                    n += 1;
                }
                batch.truncate(n);
                let done = n < batch_size;
                if let Some((mates, mates_back)) = &mates {
                    // the mate reader sends batches of batch_size until its last, shorter one
                    let mut twos = mates.recv().unwrap_or_else(|_| Ok(Vec::new()))?;
                    if twos.len() < batch.len() {
                        return Err(bad(format!(
                            "Read-2 file ran out before read-1 file (at read '{}').",
                            name_str(&batch[twos.len()].raw1)
                        )));
                    }
                    if twos.len() > batch.len() {
                        return Err(bad("Read-1 file ran out before read-2 file.".to_string()));
                    }
                    for (rec, (s, raw)) in batch.iter_mut().zip(twos.iter_mut()) {
                        if base_name(&rec.raw1) != base_name(raw) {
                            return Err(bad(format!(
                                "Read names don't match between read files: '{}' vs '{}'. \
                                 The read files need to be in the same order.",
                                name_str(&rec.raw1),
                                name_str(raw)
                            )));
                        }
                        // the read 2 into the pair, the pair's old buffers back to the mate reader
                        std::mem::swap(rec.seq2.get_or_insert_with(Vec::new), s);
                        std::mem::swap(rec.raw2.get_or_insert_with(Vec::new), raw);
                    }
                    let _ = mates_back.send(twos);
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
    Ok((rx, fastq, spent_tx))
}

/// read 2s, as (sequence, record) batches of batch_size, the last one shorter (empty, if the
/// reads come out even), read on a thread of their own; batches sent back are refilled in place
#[allow(clippy::type_complexity)]
fn spawn_mate_reader(
    mut r2: SeqReader,
    batch_size: usize,
) -> (Receiver<io::Result<Vec<(Vec<u8>, Vec<u8>)>>>, Sender<Vec<(Vec<u8>, Vec<u8>)>>) {
    let (tx, rx) = sync_channel(2);
    let (spent_tx, spent) = channel::<Vec<(Vec<u8>, Vec<u8>)>>();
    std::thread::spawn(move || {
        let result = (|| -> io::Result<()> {
            loop {
                let mut batch = spent.try_recv().unwrap_or_else(|_| Vec::with_capacity(batch_size));
                let mut n = 0;
                while n < batch_size {
                    if n == batch.len() {
                        batch.push((Vec::new(), Vec::new()));
                    }
                    let (s, raw) = &mut batch[n];
                    if !r2.next_record(s, raw)? {
                        break;
                    }
                    n += 1;
                }
                batch.truncate(n);
                let done = n < batch_size;
                if tx.send(Ok(batch)).is_err() || done {
                    return Ok(());
                }
            }
        })();
        if let Err(e) = result {
            let _ = tx.send(Err(e));
        }
    });
    (rx, spent_tx)
}

/// a record cut to its bases from `start` to `end` (and their quality scores), sequence on one
/// line, keeping its header and the bases' case
pub fn cut_record(raw: &mut Vec<u8>, start: usize, end: usize, fastq: bool, scratch: &mut Vec<u8>) {
    // the part of `line` (whose first base is base `at` of the read) between start and end
    fn take(line: &[u8], at: usize, start: usize, end: usize) -> &[u8] {
        let from = start.saturating_sub(at).min(line.len());
        let to = end.saturating_sub(at).min(line.len());
        &line[from..to.max(from)]
    }
    scratch.clear();
    let mut lines = raw.split(|&b| b == b'\n');
    scratch.extend_from_slice(lines.next().unwrap_or(&[]));
    scratch.push(b'\n');
    let mut n = 0;
    let mut qual_from = None;
    for line in lines.by_ref() {
        if fastq && line.first() == Some(&b'+') {
            qual_from = Some(());
            break;
        }
        let line = line.trim_ascii();
        scratch.extend_from_slice(take(line, n, start, end));
        n += line.len();
    }
    scratch.push(b'\n');
    if qual_from.is_some() {
        scratch.extend_from_slice(b"+\n");
        let mut q = 0;
        for line in lines {
            let line = line.trim_ascii();
            scratch.extend_from_slice(take(line, q, start, end));
            q += line.len();
        }
        scratch.push(b'\n');
    }
    std::mem::swap(raw, scratch);
}

/// `text` as one gzip member at the outputs' level, for compressing on the thread that made the
/// text and handing to GzOut::push_member (quicker than leaving it to the output's own thread
/// where text comes from many threads at once)
pub fn gz_member(text: &[u8]) -> Vec<u8> {
    if text.is_empty() {
        return Vec::new();
    }
    let mut gz = GzEncoder::new(Vec::with_capacity(text.len() / 3), Compression::new(GZ_LEVEL));
    gz.write_all(text).expect("writing to memory");
    gz.finish().expect("writing to memory")
}

/// `text` as gzip members of up to GZ_PIECE bytes each, compressed across the threads, joined
/// (for push_member; text can be cut anywhere, as joined members read back as one stream)
pub fn gz_members(text: &[u8]) -> Vec<u8> {
    if text.len() <= GZ_PIECE {
        return gz_member(text);
    }
    let parts: Vec<Vec<u8>> = text.par_chunks(GZ_PIECE).map(gz_member).collect();
    parts.concat()
}

/// a gzipped output, compressed GZ_CHUNK at a time on several threads (see GZ_CHUNK)
pub struct GzOut {
    file: BufWriter<File>,
    pending: Vec<u8>,
    // compressed once this much is pending
    flush_at: usize,
    wrote: bool,
}

impl GzOut {
    pub fn create(path: &str) -> io::Result<Self> {
        Self::with_buffer(path, GZ_PENDING)
    }

    /// one that holds at most about `flush_at` bytes before compressing them (for when many
    /// outputs are open at once, e.g., one per reference)
    pub fn with_buffer(path: &str, flush_at: usize) -> io::Result<Self> {
        Ok(GzOut {
            file: BufWriter::with_capacity(1 << 16, File::create(path)?),
            pending: Vec::new(),
            flush_at,
            wrote: false,
        })
    }

    pub fn push(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.pending.extend_from_slice(bytes);
        if self.pending.len() >= self.flush_at {
            self.flush()?;
        }
        Ok(())
    }

    /// appends an already compressed gzip member (see gz_member), after any text pending
    pub fn push_member(&mut self, member: &[u8]) -> io::Result<()> {
        if member.is_empty() {
            return Ok(());
        }
        self.flush()?;
        self.file.write_all(member)?;
        self.wrote = true;
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
                let mut gz = GzEncoder::new(Vec::with_capacity(chunk.len() / 3), Compression::new(GZ_LEVEL));
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

    pub fn finish(mut self) -> io::Result<()> {
        self.flush()?;
        if !self.wrote {
            // an empty output is still a valid (empty) gzip file
            self.file.write_all(&GzEncoder::new(Vec::new(), Compression::new(GZ_LEVEL)).finish()?)?;
        }
        self.file.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_cut() {
        let mut scratch = Vec::new();
        let mut r = b"@r1 x\nACgtACGT\n+r1\nABCDEFGH\n".to_vec();
        cut_record(&mut r, 0, 5, true, &mut scratch);
        assert_eq!(r, b"@r1 x\nACgtA\n+\nABCDE\n");
        // multi-line records come out on one line each
        let mut r = b"@r1\nACG\nTAC\n+\nABCD\nEF\n".to_vec();
        cut_record(&mut r, 0, 4, true, &mut scratch);
        assert_eq!(r, b"@r1\nACGT\n+\nABCD\n");
        let mut r = b">r1\nACGT\nACGT\n".to_vec();
        cut_record(&mut r, 0, 6, false, &mut scratch);
        assert_eq!(r, b">r1\nACGTAC\n");
        // a stretch from inside the read, across line breaks
        let mut r = b"@r1\nACG\nTAC\n+\nABCD\nEF\n".to_vec();
        cut_record(&mut r, 2, 5, true, &mut scratch);
        assert_eq!(r, b"@r1\nGTA\n+\nCDE\n");
        let mut r = b">r1\nACGT\nACGT\n".to_vec();
        cut_record(&mut r, 5, 8, false, &mut scratch);
        assert_eq!(r, b">r1\nCGT\n");
    }

    #[test]
    fn base_names() {
        assert_eq!(base_name(b"@r1/1 extra\n"), b"r1");
        assert_eq!(base_name(b"@r1/2\n"), b"r1");
        assert_eq!(base_name(b">contig_3 len=4\n"), b"contig_3");
        assert_eq!(base_name(b"@SRR1.1 1/1\n"), b"SRR1.1");
    }
}
