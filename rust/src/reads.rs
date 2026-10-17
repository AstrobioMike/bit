//! Read input: sequences from FASTA or FASTQ, gzipped or not (names and qualities aren't
//! needed). Sequences come back uppercased. Adapted from shifia's reader.

use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// how far a pass over read files has got, for progress bars: bytes of the files on disk
/// read so far, and reads
#[derive(Default)]
pub struct Progress {
    pub bytes: AtomicU64,
    pub reads: AtomicU64,
}

/// a file that counts the bytes read from it into a Progress
struct Counting {
    inner: File,
    progress: Option<Arc<Progress>>,
}

impl Read for Counting {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if let Some(p) = &self.progress {
            p.bytes.fetch_add(n as u64, Ordering::Relaxed);
        }
        Ok(n)
    }
}

/// decompresses on a thread of its own, a chunk ahead, so inflating and parsing (and
/// everything done with the reads) overlap rather than take turns
struct Inflating {
    chunks: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    current: Vec<u8>,
    at: usize,
}

impl Inflating {
    const CHUNK: usize = 4 << 20;

    fn new<R: Read + Send + 'static>(mut source: R) -> Self {
        let (tx, rx) = std::sync::mpsc::sync_channel(4);
        std::thread::spawn(move || {
            loop {
                let mut buf = vec![0u8; Self::CHUNK];
                let mut filled = 0;
                let result = loop {
                    match source.read(&mut buf[filled..]) {
                        Ok(0) => break Ok(()),
                        Ok(n) => {
                            filled += n;
                            if filled == buf.len() {
                                break Ok(());
                            }
                        }
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                        Err(e) => break Err(e),
                    }
                };
                buf.truncate(filled);
                let last = filled < Self::CHUNK;
                let failed = result.is_err();
                let msg = result.map(|_| buf);
                // the reader is gone (dropped early): stop
                if tx.send(msg).is_err() || last || failed {
                    break;
                }
            }
        });
        Inflating {
            chunks: rx,
            current: Vec::new(),
            at: 0,
        }
    }
}

impl Read for Inflating {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.at == self.current.len() {
            match self.chunks.recv() {
                Ok(Ok(chunk)) => {
                    if chunk.is_empty() {
                        return Ok(0);
                    }
                    self.current = chunk;
                    self.at = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0), // all sent
            }
        }
        let n = buf.len().min(self.current.len() - self.at);
        buf[..n].copy_from_slice(&self.current[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

#[derive(PartialEq)]
enum Format {
    Fasta,
    Fastq,
    Empty,
}

pub struct SeqReader {
    input: Box<dyn BufRead + Send>,
    format: Format,
    line: Vec<u8>,
    pending_header: bool, // fasta: a header line was just read
    path: String,
    progress: Option<Arc<Progress>>,
}

fn bad(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl SeqReader {
    pub fn open(path: &str, progress: Option<Arc<Progress>>) -> io::Result<Self> {
        let file = Counting {
            inner: File::open(path)?,
            progress: progress.clone(),
        };
        let mut raw = BufReader::with_capacity(1 << 20, file);
        let gz = raw.fill_buf()?.starts_with(&[0x1f, 0x8b]);
        let mut input: Box<dyn BufRead + Send> = if gz {
            Box::new(BufReader::with_capacity(
                1 << 20,
                Inflating::new(MultiGzDecoder::new(raw)),
            ))
        } else {
            Box::new(raw)
        };
        // past any leading blank lines
        let format = loop {
            let head = input.fill_buf()?;
            match head.first() {
                None => break Format::Empty,
                Some(b'\n') | Some(b'\r') | Some(b' ') | Some(b'\t') => input.consume(1),
                Some(b'>') => break Format::Fasta,
                Some(b'@') => break Format::Fastq,
                Some(_) => {
                    return Err(bad(format!(
                        "'{path}' doesn't look like fasta or fastq format."
                    )));
                }
            }
        };
        let mut r = SeqReader {
            input,
            format,
            line: Vec::new(),
            pending_header: false,
            path: path.to_string(),
            progress,
        };
        if r.format == Format::Fasta {
            r.read_line()?;
            r.pending_header = true;
        }
        Ok(r)
    }

    fn read_line(&mut self) -> io::Result<usize> {
        self.line.clear();
        self.input.read_until(b'\n', &mut self.line)
    }

    /// the next sequence into `seq` (cleared first); false at the end of the file
    pub fn next_seq(&mut self, seq: &mut Vec<u8>) -> io::Result<bool> {
        seq.clear();
        let got = match self.format {
            Format::Empty => false,
            Format::Fasta => self.next_fasta(seq)?,
            Format::Fastq => self.next_fastq(seq)?,
        };
        seq.make_ascii_uppercase();
        Ok(got)
    }

    fn next_fasta(&mut self, seq: &mut Vec<u8>) -> io::Result<bool> {
        if !self.pending_header {
            return Ok(false);
        }
        self.pending_header = false;
        loop {
            if self.read_line()? == 0 {
                break;
            }
            if self.line.first() == Some(&b'>') {
                self.pending_header = true;
                break;
            }
            seq.extend_from_slice(self.line.trim_ascii());
        }
        Ok(true)
    }

    fn next_fastq(&mut self, seq: &mut Vec<u8>) -> io::Result<bool> {
        loop {
            if self.read_line()? == 0 {
                return Ok(false);
            }
            if !self.line.trim_ascii().is_empty() {
                break;
            }
        }
        if self.line.first() != Some(&b'@') {
            let near = String::from_utf8_lossy(&self.line[..self.line.len().min(60)])
                .trim()
                .to_string();
            return Err(bad(format!(
                "Unexpected fastq format in '{}' near: {near}",
                self.path
            )));
        }
        // sequence line(s) up to the '+' line, then as many quality characters as bases
        loop {
            if self.read_line()? == 0 {
                return Err(bad(format!(
                    "A fastq record in '{}' ends early.",
                    self.path
                )));
            }
            if self.line.first() == Some(&b'+') {
                break;
            }
            seq.extend_from_slice(self.line.trim_ascii());
        }
        let mut quals = 0;
        while quals < seq.len() {
            if self.read_line()? == 0 {
                return Err(bad(format!(
                    "A fastq record in '{}' ends early.",
                    self.path
                )));
            }
            quals += self.line.trim_ascii().len();
        }
        Ok(true)
    }

    /// up to n sequences; fewer only at the end of the file
    pub fn batch(&mut self, n: usize) -> io::Result<Vec<Vec<u8>>> {
        let mut out = Vec::with_capacity(n);
        let mut seq = Vec::new();
        while out.len() < n && self.next_seq(&mut seq)? {
            out.push(std::mem::take(&mut seq));
        }
        if let Some(p) = &self.progress {
            p.reads.fetch_add(out.len() as u64, Ordering::Relaxed);
        }
        Ok(out)
    }
}
