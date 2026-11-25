//! bit's compiled core: the per-read loops of `bit find-overrep-seqs` (reading, adapter trimming,
//! k-mer counting, gathering reads around seeds, and placing found sequences in reads
//! allowing for edits), of `bit filter-foldbacks` (finding fold-back reads and pairs, short or long),
//! of `bit match-reads` (assigning reads to references), and the BAM and per-base-coverage
//! passes of `bit mapped-read-stats` and `bit cov-stats`, and FASTA passes of `bit fasta` subcommands
//! and `bit summarize-assembly`, and the per-read loops of `bit gen-reads`.
//! Everything else stays in Python (bit/modules/find_overrep.py,
//! bit/modules/filter_foldbacks.py), which calls into this as `bit._native`.

mod anchors;
mod bam;
mod edlib;
mod fasta;
mod foldback;
mod gen_reads;
mod kmers;
mod match_reads;
mod matcher;
mod readio;
mod reads;
mod trim;

use numpy::PyArray1;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

fn to_py_err(e: io::Error) -> PyErr {
    if e.kind() == io::ErrorKind::InvalidData {
        PyValueError::new_err(e.to_string())
    } else {
        PyIOError::new_err(e.to_string())
    }
}

/// How far a pass over the read files has got (bytes of the files on disk read so far, and
/// reads), for progress bars; pass one to a function as `progress` and read it from another
/// thread while that runs.
#[pyclass(frozen)]
struct Progress {
    inner: Arc<reads::Progress>,
}

#[pymethods]
impl Progress {
    #[new]
    fn new() -> Self {
        Progress {
            inner: Arc::new(reads::Progress::default()),
        }
    }

    #[getter]
    fn bytes(&self) -> u64 {
        self.inner.bytes.load(Ordering::Relaxed)
    }

    #[getter]
    fn reads(&self) -> u64 {
        self.inner.reads.load(Ordering::Relaxed)
    }
}

fn arc(progress: &Option<PyRef<Progress>>) -> Option<Arc<reads::Progress>> {
    progress.as_ref().map(|p| p.inner.clone())
}

fn trimmer(adapters: Option<Vec<(String, String)>>) -> Option<trim::Trimmer> {
    adapters.map(|a| trim::Trimmer::new(&a))
}

/// sets how many threads the compiled core uses (once per process; later calls change
/// nothing and return false)
#[pyfunction]
fn set_threads(n: usize) -> bool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n.max(1))
        .build_global()
        .is_ok()
}

/// (trimmed read, what cut it or None): the read cut where adapter read-through starts
#[pyfunction]
fn trim_read(seq: &str, adapters: Vec<(String, String)>) -> (String, Option<String>) {
    let t = trim::Trimmer::new(&adapters);
    let (keep, cut) = t.trim(seq.as_bytes());
    (seq[..keep].to_string(), cut.map(|c| t.cut_name(&c)))
}

/// (forward k-mer counts per file, 3-mer counts, 2-mer counts, valid k-mer windows at each
/// read position, reads) over all reads, trimmed with `adapters` if given
#[pyfunction]
#[pyo3(signature = (paths, k, adapters=None, progress=None, max_reads=0))]
#[allow(clippy::type_complexity)]
fn count_kmers<'py>(
    py: Python<'py>,
    paths: Vec<String>,
    k: usize,
    adapters: Option<Vec<(String, String)>>,
    progress: Option<PyRef<Progress>>,
    max_reads: u64,
) -> PyResult<(
    Vec<Bound<'py, PyArray1<u64>>>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u64>>,
    u64,
)> {
    let t = trimmer(adapters);
    let pr = arc(&progress);
    let c = py
        .detach(|| kmers::count(&paths, k, t.as_ref(), pr.as_ref(), max_reads))
        .map_err(to_py_err)?;
    Ok((
        c.by_file
            .into_iter()
            .map(|v| PyArray1::from_vec(py, v))
            .collect(),
        PyArray1::from_vec(py, c.tri),
        PyArray1::from_vec(py, c.di),
        PyArray1::from_vec(py, c.windows_at),
        c.n_reads,
    ))
}

/// up to max_contexts reads holding each seed k-mer, oriented so the seed reads forward, as
/// (read, seed position)
#[pyfunction]
#[pyo3(signature = (paths, k, seeds, max_contexts, adapters=None, progress=None, max_reads=0))]
#[allow(clippy::too_many_arguments)]
fn collect_contexts(
    py: Python<'_>,
    paths: Vec<String>,
    k: usize,
    seeds: Vec<String>,
    max_contexts: usize,
    adapters: Option<Vec<(String, String)>>,
    progress: Option<PyRef<Progress>>,
    max_reads: u64,
) -> PyResult<Vec<Vec<(String, usize)>>> {
    let t = trimmer(adapters);
    let seeds: Vec<Vec<u8>> = seeds.into_iter().map(|s| s.into_bytes()).collect();
    let pr = arc(&progress);
    let got = py
        .detach(|| {
            kmers::contexts(
                &paths,
                k,
                t.as_ref(),
                &seeds,
                max_contexts,
                pr.as_ref(),
                max_reads,
            )
        })
        .map_err(to_py_err)?;
    Ok(got
        .into_iter()
        .map(|v| {
            v.into_iter()
                .map(|(r, p)| (String::from_utf8_lossy(&r).into_owned(), p))
                .collect()
        })
        .collect())
}

/// Places found sequences in reads, allowing for edits (see Matcher in find_overrep.py).
#[pyclass(frozen)]
struct Matcher {
    inner: matcher::Matcher,
}

#[pymethods]
impl Matcher {
    #[new]
    #[pyo3(signature = (seqs, max_diff, min_read_len, min_overlap=0, groups=None))]
    fn new(
        seqs: Vec<String>,
        max_diff: f64,
        min_read_len: usize,
        min_overlap: usize,
        groups: Option<Vec<usize>>,
    ) -> PyResult<Self> {
        if let Some(g) = &groups
            && g.len() != seqs.len()
        {
            return Err(PyValueError::new_err("groups must be as long as seqs"));
        }
        let seqs: Vec<Vec<u8>> = seqs.into_iter().map(|s| s.into_bytes()).collect();
        Ok(Matcher {
            inner: matcher::Matcher::new(&seqs, max_diff, min_read_len, min_overlap, groups.as_deref()),
        })
    }

    /// {seq index: (start in read or None if the read is within it or overlaps an end of
    /// it, is reverse complement, edits, the sequence's bases the read covers)};
    /// `prefilter=False` tries every sequence (for checking the prefilter)
    #[pyo3(signature = (read, prefilter=true))]
    fn hits<'py>(
        &self,
        py: Python<'py>,
        read: &str,
        prefilter: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for h in self.inner.hits(read.as_bytes(), !prefilter) {
            d.set_item(h.seq, (h.start, h.is_rc, h.edits, h.bases))?;
        }
        Ok(d)
    }

    /// counts the found sequences in reads (one file, or two for pairs): ([per-file tally
    /// dicts], pair tally dict or None)
    #[pyo3(signature = (paths, adapters, min_insert_len, progress=None))]
    #[allow(clippy::type_complexity)]
    fn count_in_reads<'py>(
        &self,
        py: Python<'py>,
        paths: Vec<String>,
        adapters: Vec<(String, String)>,
        min_insert_len: usize,
        progress: Option<PyRef<Progress>>,
    ) -> PyResult<(Vec<Bound<'py, PyDict>>, Option<Bound<'py, PyDict>>)> {
        let t = trim::Trimmer::new(&adapters);
        let pr = arc(&progress);
        let (tallies, pairs) = py
            .detach(|| {
                matcher::count_in_reads(&paths, &self.inner, &t, min_insert_len, pr.as_ref())
            })
            .map_err(to_py_err)?;
        let mut out = Vec::new();
        for t in tallies {
            let d = PyDict::new(py);
            d.set_item("reads", t.reads)?;
            d.set_item("reads_with", t.reads_with)?;
            d.set_item("full", t.full)?;
            d.set_item("rc", t.rc)?;
            d.set_item(
                "edits",
                t.edits
                    .into_iter()
                    .map(|m| m.into_iter().collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "bases",
                t.bases
                    .into_iter()
                    .map(|m| m.into_iter().collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "starts",
                t.starts
                    .into_iter()
                    .map(|m| m.into_iter().collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            )?;
            d.set_item("adapters", t.adapters.into_iter().collect::<Vec<_>>())?;
            d.set_item("adapter_only", t.adapter_only)?;
            d.set_item("explained", t.explained)?;
            out.push(d);
        }
        let pairs = match pairs {
            Some(p) => {
                let d = PyDict::new(py);
                d.set_item("pairs", p.pairs)?;
                d.set_item("pairs_with", p.pairs_with)?;
                d.set_item("pairs_both", p.pairs_both)?;
                d.set_item("explained", p.explained)?;
                d.set_item("explained_both", p.explained_both)?;
                d.set_item("adapter_only", p.adapter_only)?;
                Some(d)
            }
            None => None,
        };
        Ok((out, pairs))
    }
}

/// The matching core of `bit match-reads` (see ExactMatcher, EditMatcher and process_batch in
/// match_reads.py, which this follows step for step). With max_edits 0 it's the exact matcher.
#[pyclass(frozen)]
struct ReadMatcher {
    inner: match_reads::Matcher,
    seq_to_ref: Vec<u32>,
    pool: Arc<rayon::ThreadPool>,
}

#[pymethods]
impl ReadMatcher {
    #[new]
    #[pyo3(signature = (seqs, seq_to_ref, k, circular, max_edits, n_samples, step, threads=1))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        seqs: Vec<String>,
        seq_to_ref: Vec<u32>,
        k: usize,
        circular: bool,
        max_edits: usize,
        n_samples: usize,
        step: usize,
        threads: usize,
    ) -> PyResult<Self> {
        if seq_to_ref.len() != seqs.len() {
            return Err(PyValueError::new_err("seq_to_ref must be as long as seqs"));
        }
        if k == 0 || step == 0 || seqs.iter().any(|s| s.is_empty() || !s.is_ascii()) {
            return Err(PyValueError::new_err("sequences must be non-empty ASCII (and k, step > 0)"));
        }
        let seqs: Vec<Vec<u8>> = seqs.into_iter().map(|s| s.into_bytes()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let inner = py.detach(|| {
            pool.install(|| {
                if max_edits > 0 {
                    match_reads::Matcher::Edit(match_reads::Edit::new(seqs, k, circular, max_edits))
                } else {
                    match_reads::Matcher::Exact(match_reads::Exact::new(seqs, k, n_samples, circular, step))
                }
            })
        });
        Ok(ReadMatcher {
            inner,
            seq_to_ref,
            pool: Arc::new(pool),
        })
    }

    /// {seq index: edit distance} for one read (ASCII), as the Python matchers' distances()
    #[pyo3(signature = (seq, min_frac_of_seq=0.0))]
    fn distances<'py>(&self, py: Python<'py>, seq: &str, min_frac_of_seq: f64) -> PyResult<Bound<'py, PyDict>> {
        if !seq.is_ascii() {
            return Err(PyValueError::new_err("reads must be ASCII"));
        }
        let d = PyDict::new(py);
        let mut v = self.inner.distances(seq.as_bytes(), min_frac_of_seq);
        v.sort_unstable();
        for (i, dist) in v {
            d.set_item(i, dist)?;
        }
        Ok(d)
    }

    /// Assigns a batch of reads (or of pairs, given as consecutive mates) across threads:
    /// [(position in the batch, hit refs, hit seqs, best distance, next-best or None), ...]
    /// for those assigned to >= 1 ref, in input order. Reads must be ASCII.
    #[allow(clippy::type_complexity)]
    fn assign(
        &self,
        py: Python<'_>,
        seqs: Vec<String>,
        paired: bool,
        min_frac_of_seq: f64,
    ) -> PyResult<Vec<(usize, Vec<u32>, Vec<u32>, usize, Option<usize>)>> {
        if paired && !seqs.len().is_multiple_of(2) {
            return Err(PyValueError::new_err("paired input needs an even number of mates"));
        }
        if seqs.iter().any(|s| !s.is_ascii()) {
            return Err(PyValueError::new_err("reads must be ASCII"));
        }
        let results = py.detach(|| {
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
            self.pool.install(|| {
                match_reads::assign_batch(&self.inner, &self.seq_to_ref, &refs, paired, min_frac_of_seq)
            })
        });
        Ok(results
            .into_iter()
            .enumerate()
            .filter_map(|(j, r)| r.map(|a| (j, a.hits, a.hit_seqs, a.best, a.next_best)))
            .collect())
    }

    /// Assigns every read (or pair) of a read file (or pair of them) and writes the read-hits
    /// table, as match_reads() in match_reads.py does with the Python matchers: records are read
    /// on threads of their own (readio.rs) and matched a batch at a time across the matcher's
    /// threads; rows go to `hits_path` after `header`, sorted by read length (combined for
    /// pairs), longest first, then edit distance, ties in input order (match_reads::SortedRows,
    /// holding up to `buffer_rows` in memory). `ref_names` name the refs in the table, and
    /// `seq_labels`, if given, the seqs (the matching_seqs column). After each batch,
    /// `on_progress(reads in the batch, bytes of the input files read so far)` is called. With
    /// `reads_out` (for each ref, its output path(s) less the extension: one, or two for
    /// pairs), uniquely assigned reads are written there as they were in the input, gzipped
    /// (".fastq.gz" or ".fasta.gz" added); a ref's files are made only once a read goes to them.
    /// Returns (reads, unique, ambiguous, unique by ref, ambiguous by ref, unique by seq,
    /// ambiguous by seq).
    #[pyo3(signature = (read_1, read_2, min_frac_of_seq, hits_path, header, ref_names, seq_labels,
                        buffer_rows, on_progress, reads_out=None))]
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn assign_files(
        &self,
        py: Python<'_>,
        read_1: String,
        read_2: Option<String>,
        min_frac_of_seq: f64,
        hits_path: String,
        header: String,
        ref_names: Vec<String>,
        seq_labels: Option<Vec<String>>,
        buffer_rows: usize,
        on_progress: Bound<'_, PyAny>,
        reads_out: Option<Vec<Vec<String>>>,
    ) -> PyResult<(u64, u64, u64, Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>)> {
        let paired = read_2.is_some();
        let prog = Arc::new(reads::Progress::default());
        let (batches, fastq, spent) =
            readio::spawn_reader(read_1, read_2, MATCH_BATCH, Some(prog.clone())).map_err(to_py_err)?;
        // (a Mutex only so the receiver can be waited on with the GIL let go)
        let batches = std::sync::Mutex::new(batches);
        let writer = reads_out.map(|outs| spawn_ref_writer(outs, fastq, self.pool.clone()));
        let mut table = match_reads::SortedRows::new(&hits_path, header, buffer_rows);
        let mut t = MatchTally {
            unique_ref: vec![0; ref_names.len()],
            ambig_ref: vec![0; ref_names.len()],
            unique_seq: vec![0; self.seq_to_ref.len()],
            ambig_seq: vec![0; self.seq_to_ref.len()],
            ..Default::default()
        };
        let result = (|| -> PyResult<()> {
            loop {
                let batch = match py.detach(|| batches.lock().expect("reader").recv()) {
                    Err(_) => return Ok(()), // all read
                    Ok(Err(e)) => return Err(to_py_err(e)),
                    Ok(Ok(b)) => b,
                };
                let n = batch.len();
                let writes = py.detach(|| -> io::Result<RefWrites> {
                    let (rows, writes) = self.assign_records(&batch, paired, min_frac_of_seq, writer.is_some());
                    for row in rows {
                        let (key, line) = t.take(row, &ref_names, seq_labels.as_deref());
                        table.add(key, line)?;
                    }
                    Ok(writes)
                }).map_err(to_py_err)?;
                t.reads += n as u64;
                let _ = spent.send(batch);
                if let Some((tx, _)) = &writer
                    && !writes.is_empty()
                    && tx.send(writes).is_err()
                {
                    break; // the writer failed: its error is below
                }
                on_progress.call1((n, prog.bytes.load(Ordering::Relaxed)))?;
                if n < MATCH_BATCH {
                    return Ok(());
                }
            }
            Ok(())
        })();
        drop(batches);
        let wrote = match writer {
            Some((tx, handle)) => {
                drop(tx);
                py.detach(|| handle.join().expect("reads writer thread"))
            }
            None => Ok(()),
        };
        result?;
        wrote.map_err(to_py_err)?;
        py.detach(|| table.finish()).map_err(to_py_err)?;
        Ok((t.reads, t.unique, t.ambig, t.unique_ref, t.ambig_ref, t.unique_seq, t.ambig_seq))
    }
}

/// assign_files()'s counts
#[derive(Default)]
struct MatchTally {
    reads: u64,
    unique: u64,
    ambig: u64,
    unique_ref: Vec<u64>,
    ambig_ref: Vec<u64>,
    unique_seq: Vec<u64>,
    ambig_seq: Vec<u64>,
}

impl MatchTally {
    /// one assigned read (or pair) into the counts; its hits-table row as (sort key, line)
    fn take(&mut self, row: MatchRow, ref_names: &[String], seq_labels: Option<&[String]>) -> ((i64, i64), String) {
        use std::fmt::Write;
        let (id, lengths, total_len, hits, hit_seqs, best, next_best) = row;
        let status = if hits.len() == 1 {
            self.unique += 1;
            self.unique_ref[hits[0] as usize] += 1;
            for &si in &hit_seqs {
                self.unique_seq[si as usize] += 1;
            }
            "unique"
        } else {
            self.ambig += 1;
            for &ri in &hits {
                self.ambig_ref[ri as usize] += 1;
            }
            for &si in &hit_seqs {
                self.ambig_seq[si as usize] += 1;
            }
            "ambiguous"
        };
        let mut line = format!("{id}\t{lengths}\t{status}\t{}\t", hits.len());
        for (j, &ri) in hits.iter().enumerate() {
            if j > 0 {
                line.push(',');
            }
            line.push_str(&ref_names[ri as usize]);
        }
        let _ = write!(line, "\t{best}\t");
        match next_best {
            Some(nb) => {
                let _ = write!(line, "{nb}");
            }
            None => line.push_str("NA"),
        }
        if let Some(labels) = seq_labels {
            line.push('\t');
            for (j, &si) in hit_seqs.iter().enumerate() {
                if j > 0 {
                    line.push(',');
                }
                line.push_str(&labels[si as usize]);
            }
        }
        ((-(total_len as i64), best as i64), line)
    }
}

/// A gzipped text file written from Python (readio::GzOut: compressed in chunks, each its own
/// gzip member, on several threads; gzip files written this way can be joined end to end into
/// one). `buffer`: bytes held before compressing (default 32 MB; smaller when many are open at
/// once, e.g., one per worker process).
#[pyclass]
struct GzWriter {
    out: Option<readio::GzOut>,
}

#[pymethods]
impl GzWriter {
    #[new]
    #[pyo3(signature = (path, buffer=32 << 20))]
    fn new(path: String, buffer: usize) -> PyResult<Self> {
        let out = readio::GzOut::with_buffer(&path, buffer).map_err(to_py_err)?;
        Ok(GzWriter { out: Some(out) })
    }

    /// appends text (str or bytes)
    fn write(&mut self, py: Python<'_>, text: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = self.out.as_mut().ok_or_else(|| PyValueError::new_err("write to a closed GzWriter"))?;
        if let Ok(s) = text.extract::<pyo3::pybacked::PyBackedStr>() {
            py.detach(|| out.push(s.as_bytes())).map_err(to_py_err)
        } else {
            let b: pyo3::pybacked::PyBackedBytes = text.extract()?;
            py.detach(|| out.push(&b)).map_err(to_py_err)
        }
    }

    /// compresses what's left and closes the file (again: nothing)
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        match self.out.take() {
            Some(out) => py.detach(|| out.finish()).map_err(to_py_err),
            None => Ok(()),
        }
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        py: Python<'_>,
        _exc_type: &Bound<'_, PyAny>,
        _exc: &Bound<'_, PyAny>,
        _tb: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.close(py)?;
        Ok(false)
    }
}

// reads (or pairs) read and matched at a time by ReadMatcher.assign_files
const MATCH_BATCH: usize = 16_384;

type MatchRow = (String, String, usize, Vec<u32>, Vec<u32>, usize, Option<usize>);
// a batch's uniquely assigned reads by ref: (ref, read 1's records, read 2's)
type RefWrites = Vec<(u32, Vec<u8>, Vec<u8>)>;

impl ReadMatcher {
    /// one batch's assigned reads (or pairs), in input order: (id, lengths, total length, hit
    /// refs, hit seqs, best distance, next-best or None); and with `keep`, its uniquely assigned
    /// records by ref
    fn assign_records(&self, batch: &[readio::Rec], paired: bool, min_frac: f64, keep: bool) -> (Vec<MatchRow>, RefWrites) {
        let mut seqs: Vec<&[u8]> = Vec::with_capacity(batch.len() * if paired { 2 } else { 1 });
        for r in batch {
            seqs.push(&r.seq1);
            if let Some(s2) = &r.seq2 {
                seqs.push(s2);
            }
        }
        let results = self
            .pool
            .install(|| match_reads::assign_batch(&self.inner, &self.seq_to_ref, &seqs, paired, min_frac));
        let mut rows = Vec::new();
        let mut by_ref: std::collections::BTreeMap<u32, (Vec<u8>, Vec<u8>)> = Default::default();
        for (rec, a) in batch.iter().zip(results) {
            let Some(a) = a else { continue };
            let (id, lengths, total_len) = match &rec.seq2 {
                Some(s2) => (
                    String::from_utf8_lossy(readio::base_name(&rec.raw1)).into_owned(),
                    format!("{},{}", rec.seq1.len(), s2.len()),
                    rec.seq1.len() + s2.len(),
                ),
                None => (readio::name_str(&rec.raw1), rec.seq1.len().to_string(), rec.seq1.len()),
            };
            if keep && a.hits.len() == 1 {
                let (one, two) = by_ref.entry(a.hits[0]).or_default();
                one.extend_from_slice(&rec.raw1);
                if let Some(raw2) = &rec.raw2 {
                    two.extend_from_slice(raw2);
                }
            }
            rows.push((id, lengths, total_len, a.hits, a.hit_seqs, a.best, a.next_best));
        }
        (rows, by_ref.into_iter().map(|(r, (a, b))| (r, a, b)).collect())
    }
}

/// writes uniquely assigned reads to their ref's gzipped file(s), on a thread of its own,
/// compressing each batch's on the matcher's threads (in pieces, as most of a batch's reads can
/// go to one ref) while the next batch is matched; a ref's files are made when its first read
/// comes
#[allow(clippy::type_complexity)]
fn spawn_ref_writer(
    outs: Vec<Vec<String>>,
    fastq: bool,
    pool: Arc<rayon::ThreadPool>,
) -> (std::sync::mpsc::SyncSender<RefWrites>, std::thread::JoinHandle<io::Result<()>>) {
    let (tx, rx) = std::sync::mpsc::sync_channel::<RefWrites>(2);
    let ext = if fastq { "fastq.gz" } else { "fasta.gz" };
    let handle = std::thread::spawn(move || -> io::Result<()> {
        let mut files: Vec<Option<Vec<readio::GzOut>>> = (0..outs.len()).map(|_| None).collect();
        for writes in rx {
            let writes: Vec<(u32, Vec<u8>, Vec<u8>)> = pool.install(|| {
                writes
                    .into_iter()
                    .map(|(r, a, b)| (r, readio::gz_members(&a), readio::gz_members(&b)))
                    .collect()
            });
            for (r, one, two) in writes {
                let r = r as usize;
                if files[r].is_none() {
                    let made: io::Result<Vec<_>> = outs[r]
                        .iter()
                        .map(|p| readio::GzOut::with_buffer(&format!("{p}.{ext}"), 1 << 20))
                        .collect();
                    files[r] = Some(made?);
                }
                let f = files[r].as_mut().unwrap();
                f[0].push_member(&one)?;
                if let Some(f2) = f.get_mut(1) {
                    f2.push_member(&two)?;
                }
            }
        }
        for f in files.into_iter().flatten() {
            for o in f {
                o.finish()?;
            }
        }
        Ok(())
    });
    (tx, handle)
}

type AccTuple = (u64, f64, Option<f64>, Option<f64>, Vec<(f64, u64)>);

fn acc_tuple(a: bam::Acc) -> AccTuple {
    let values = a.order.iter().map(|v| (*v, a.values[&v.to_bits()])).collect();
    (a.count, a.total, a.min, a.max, values)
}

type ReadRowTuple = (String, f64, f64, i64, i64, i64, i64, f64, i64, i64);

/// get_mapped_reads_pids() over a BAM (see bam.rs): (per-reference read rows or None,
/// gap-aware stats, gap-compressed stats, clip totals); stats are (count, total, min, max,
/// [(distinct value, count), ...] in first-seen order). Raises ValueError for anything the
/// Python (pysam) path should handle instead.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn bam_mapped_reads_pids(
    py: Python<'_>,
    path: String,
    include_non_primary: bool,
    store_rows: bool,
) -> PyResult<(
    Option<Vec<(String, Vec<ReadRowTuple>)>>,
    AccTuple,
    AccTuple,
    (u64, i64, i64, i64, f64, u64, f64, Option<i64>, Option<i64>, Vec<(i64, u64)>),
)> {
    let (rows, pid, gc, c) = py
        .detach(|| bam::mapped_reads_pids(&path, include_non_primary, store_rows))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let rows = rows.map(|groups| {
        groups
            .into_iter()
            .map(|(name, rows)| {
                (
                    name,
                    rows.into_iter()
                        .map(|r| {
                            (
                                r.read_id,
                                r.pid,
                                r.gc_pid,
                                r.query_aligned,
                                r.full_read_length,
                                r.soft,
                                r.hard,
                                r.clipped_frac,
                                r.ref_start,
                                r.ref_end,
                            )
                        })
                        .collect(),
                )
            })
            .collect()
    });
    let lens = c.len_order.iter().map(|l| (*l, c.len_values[l])).collect();
    Ok((
        rows,
        acc_tuple(pid),
        acc_tuple(gc),
        (
            c.n_reads,
            c.total_aligned,
            c.total_soft,
            c.total_hard,
            c.clip_frac_sum,
            c.len_count,
            c.len_total,
            c.len_min,
            c.len_max,
            lens,
        ),
    ))
}

/// get_per_contig_pid_stats() over a BAM: [(reference name, stats), ...] in first-seen order
/// (stats as in bam_mapped_reads_pids). Raises ValueError for anything the Python (pysam)
/// path should handle instead.
#[pyfunction]
fn bam_per_contig_pid_stats(
    py: Python<'_>,
    path: String,
    include_non_primary: bool,
    gap_compressed: bool,
) -> PyResult<Vec<(String, AccTuple)>> {
    let accs = py
        .detach(|| bam::per_contig_pid_stats(&path, include_non_primary, gap_compressed))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(accs.into_iter().map(|(n, a)| (n, acc_tuple(a))).collect())
}

/// cov-stats' pass over a bam: (per wanted contig, None (not in the bam) or its coverage
/// totals (total coverage count, bases detected at all, at 10x, [(depth, bases), ...]); the per-contig percent-ID stats as
/// bam_per_contig_pid_stats gives them, or None if `pids` is false). `exclude` is the flag
/// bits of alignments that add no coverage. Raises ValueError for anything the Python
/// (pysam) path should handle instead, and for a bam whose alignments aren't grouped by
/// reference (the message then says so). Alignments under `min_pid` percent identity
/// (gap-compressed; 0: none left out) add to neither the coverage nor the percent IDs.
#[pyfunction]
#[pyo3(signature = (path, wanted, exclude, pids, include_non_primary, gap_compressed, min_pid=0.0, progress=None))]
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn bam_coverage_and_pid_stats(
    py: Python<'_>,
    path: String,
    wanted: Vec<String>,
    exclude: u16,
    pids: bool,
    include_non_primary: bool,
    gap_compressed: bool,
    min_pid: f64,
    progress: Option<PyRef<Progress>>,
) -> PyResult<(Vec<Option<(u128, u128, u128, Vec<(i64, u128)>)>>, Option<Vec<(String, AccTuple)>>)> {
    let pr = arc(&progress);
    let pid = pids.then_some((include_non_primary, gap_compressed));
    let (cov, accs) = py
        .detach(|| bam::coverage_and_pids(&path, &wanted, exclude, pid, min_pid, pr.as_ref()))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok((
        cov.into_iter()
            .map(|o| o.map(|(a, b, c, h)| (a, b, c, h.into_iter().collect())))
            .collect(),
        accs.map(|accs| accs.into_iter().map(|(n, a)| (n, acc_tuple(a))).collect()),
    ))
}

fn val_err(e: io::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// [(id, length, G+C count), ...] for each record, as Biopython's SeqIO "fasta" reads them
/// (ValueError: use the Python path)
#[pyfunction]
fn fasta_records_summary(py: Python<'_>, path: String) -> PyResult<Vec<(String, usize, usize)>> {
    py.detach(|| fasta::records_summary(&path)).map_err(val_err)
}

/// [(name, length, G+C count, [G+C count of each window]), ...] (see fasta::gc_windows)
#[pyfunction]
#[allow(clippy::type_complexity)]
fn fasta_gc_windows(
    py: Python<'_>,
    path: String,
    window: usize,
    step: usize,
) -> PyResult<Vec<(String, usize, usize, Vec<u64>)>> {
    if window == 0 || step == 0 {
        return Err(PyValueError::new_err("window and step need to be positive"));
    }
    py.detach(|| fasta::gc_windows(&path, window, step)).map_err(val_err)
}

/// filter-by-length: writes kept records to out_path; (seqs, kept, bases, bases kept)
#[pyfunction]
#[pyo3(signature = (path, out_path, min_length, max_length))]
fn fasta_filter_by_length(
    py: Python<'_>,
    path: String,
    out_path: String,
    min_length: Option<i64>,
    max_length: i64,
) -> PyResult<(u64, u64, u64, u64)> {
    py.detach(|| fasta::filter_by_length(&path, &out_path, min_length, max_length))
        .map_err(val_err)
}

/// remove-wraps of path, written to the open file descriptor `fd` (the input is checked in
/// full first, so nothing is written if it calls for the Python path); True when done
#[pyfunction]
fn fasta_remove_wraps(py: Python<'_>, path: String, fd: i32) -> PyResult<bool> {
    use std::os::fd::BorrowedFd;
    py.detach(|| -> io::Result<()> {
        fasta::check_plain(&path)?;
        // a duplicate of the descriptor, so closing ours leaves the caller's open
        let dup = unsafe { BorrowedFd::borrow_raw(fd) }.try_clone_to_owned()?;
        let mut out = std::fs::File::from(dup);
        fasta::remove_wraps(&path, &mut out)
    })
    .map_err(val_err)?;
    Ok(true)
}

/// extract-by-primers' search (see fasta::primer_hits): [(header, uppercased seq, [[(start,
/// end exclusive), ...] per primer]), ...], with `equalities` as (code, base) pairs
#[pyfunction]
#[allow(clippy::type_complexity)]
fn fasta_primer_sites(
    py: Python<'_>,
    path: String,
    primers: Vec<String>,
    k: usize,
    equalities: Vec<(String, String)>,
    threads: usize,
) -> PyResult<Vec<(String, String, Vec<Vec<(usize, usize)>>)>> {
    if primers.iter().any(|p| !p.is_ascii()) {
        return Err(PyValueError::new_err("primers must be ASCII"));
    }
    let mut eq = Vec::with_capacity(equalities.len());
    for (a, b) in &equalities {
        match (a.as_bytes(), b.as_bytes()) {
            ([a], [b]) => eq.push((*a, *b)),
            _ => return Err(PyValueError::new_err("equalities must be single ASCII characters")),
        }
    }
    let primers: Vec<Vec<u8>> = primers.into_iter().map(|p| p.into_bytes()).collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    py.detach(|| pool.install(|| fasta::primer_hits(&path, &primers, k, &eq)))
        .map_err(val_err)
}

/// summarize-assembly's counts for a clean FASTA: (records, total length, {char: count},
/// [lengths]) (ValueError: use pyfastx)
#[pyfunction]
#[allow(clippy::type_complexity)]
fn fasta_assembly_summary(
    py: Python<'_>,
    path: String,
) -> PyResult<(u64, u64, Vec<(String, u64)>, Vec<u64>)> {
    let (n, size, comp, lengths) = py.detach(|| fasta::assembly_summary(&path)).map_err(val_err)?;
    let comp = comp
        .into_iter()
        .enumerate()
        .filter(|&(_, c)| c > 0)
        .map(|(b, c)| ((b as u8 as char).to_string(), c))
        .collect();
    Ok((n, size, comp, lengths))
}

fn mt_from(state: Vec<u32>, index: usize, gauss_next: Option<f64>) -> PyResult<gen_reads::Mt> {
    let state: [u32; 624] = state
        .try_into()
        .map_err(|_| PyValueError::new_err("the generator state needs 624 words"))?;
    if index > 624 {
        return Err(PyValueError::new_err("bad generator index"));
    }
    Ok(gen_reads::Mt {
        state,
        index,
        gauss_next,
    })
}

/// Draws from a random.Random state (as getstate() gives it: 624 words, index, gauss_next)
/// for checking this generator against Python's (the tests do, so a seed keeps giving the reads it did): each op is (0, _, _) for random(), (1, a, b)
/// for randint(a, b), (2, mu, sigma) for gauss(mu, sigma). Returns the values (randint's as
/// floats) and the state after.
/// each primer's sites in one sequence: [[(start, end exclusive), ...] per primer] (see fasta::seq_primer_sites)
#[pyfunction]
fn primer_sites_in_seq(seq: String, primers: Vec<String>, k: usize, equalities: Vec<(String, String)>)
                       -> PyResult<Vec<Vec<(usize, usize)>>> {
    if !seq.is_ascii() || primers.iter().any(|p| !p.is_ascii()) {
        return Err(PyValueError::new_err("the sequence and primers must be ASCII"));
    }
    let mut eq = Vec::with_capacity(equalities.len());
    for (a, b) in &equalities {
        match (a.as_bytes(), b.as_bytes()) {
            ([a], [b]) => eq.push((*a, *b)),
            _ => return Err(PyValueError::new_err("equalities must be single ASCII characters")),
        }
    }
    let primers: Vec<Vec<u8>> = primers.into_iter().map(|p| p.into_bytes()).collect();
    fasta::seq_primer_sites(seq.as_bytes(), &primers, k, &eq).map_err(val_err)
}

#[pyfunction]
#[allow(clippy::type_complexity)]
fn rng_sample(
    state: Vec<u32>,
    index: usize,
    gauss_next: Option<f64>,
    ops: Vec<(u8, f64, f64)>,
) -> PyResult<(Vec<f64>, Vec<u32>, usize, Option<f64>)> {
    let mut mt = mt_from(state, index, gauss_next)?;
    let mut out = Vec::with_capacity(ops.len());
    for (kind, a, b) in ops {
        out.push(match kind {
            0 => mt.random(),
            1 => mt
                .randint(a as i64, b as i64)
                .ok_or_else(|| PyValueError::new_err("empty range"))? as f64,
            _ => mt.gauss(a, b),
        });
    }
    Ok((out, mt.state.to_vec(), mt.index, mt.gauss_next))
}

// reads (or pairs) generated at a time by gen_reads_contigs before they're written out
const GEN_CHUNK: u64 = 100_000;

/// The contig loop of gen_one_file_paired / gen_one_file_single in gen_reads.py (see
/// gen_reads.rs), from contig `first` on: each contig's `budgets` reads (or pairs) are made
/// from the random.Random state given (624 words, index, gauss_next) and written straight to
/// the gzipped writers (`r2` for pairs; `tsv` for the per-read sources, whose presence also
/// keeps the coordinates out of the read headers). `seqs` are the contigs' sequences, `ids`
/// their names. A contig that isn't a nucleotide sequence (non-ASCII), or settings that can't
/// make reads from it, raise ValueError. Returns (the read count after, reads written, bases in
/// them, each contig's covered spans (merged, sorted), and the generator state after).
#[pyfunction]
#[pyo3(signature = (seqs, ids, budgets, first, read_count, state, index, gauss_next, paired, long,
                    read_length, mu, sigma, min_len, max_len, circular, include_ns, source_fasta,
                    r1, r2=None, tsv=None))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn gen_reads_contigs(
    py: Python<'_>,
    seqs: Vec<pyo3::pybacked::PyBackedStr>,
    ids: Vec<String>,
    budgets: Vec<u64>,
    first: usize,
    read_count: u64,
    state: Vec<u32>,
    index: usize,
    gauss_next: Option<f64>,
    paired: bool,
    long: bool,
    read_length: i64,
    mu: f64,
    sigma: f64,
    min_len: i64,
    max_len: i64,
    circular: bool,
    include_ns: bool,
    source_fasta: String,
    r1: Bound<'_, GzWriter>,
    r2: Option<Bound<'_, GzWriter>>,
    tsv: Option<Bound<'_, GzWriter>>,
) -> PyResult<(
    u64,
    u64,
    u64,
    Vec<(usize, Vec<(i64, i64)>)>,
    Vec<u32>,
    usize,
    Option<f64>,
)> {
    if ids.len() != seqs.len() || budgets.len() != seqs.len() || paired != r2.is_some() {
        return Err(PyValueError::new_err("seqs, ids and budgets must match, and r2 is for pairs only"));
    }
    let mut mt = mt_from(state, index, gauss_next)?;
    // the writers' outputs, taken for the run and put back after
    let take = |w: &Bound<'_, GzWriter>| -> PyResult<readio::GzOut> {
        w.borrow_mut().out.take().ok_or_else(|| PyValueError::new_err("write to a closed GzWriter"))
    };
    let mut out1 = take(&r1)?;
    let mut out2 = match &r2 { Some(w) => Some(take(w)?), None => None };
    let mut out_tsv = match &tsv { Some(w) => Some(take(w)?), None => None };
    let per_read_tsv = out_tsv.is_some();
    let result = py.detach(|| -> io::Result<_> {
        let (mut count, mut written, mut bases) = (read_count, 0u64, 0u64);
        let mut spans_by_contig = Vec::new();
        for ci in first..seqs.len() {
            let budget = budgets[ci];
            if budget == 0 {
                continue;
            }
            let seq: &str = &seqs[ci];
            if !seq.is_ascii() {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("contig '{}' holds characters that aren't a nucleotide sequence's", ids[ci])));
            }
            let p = gen_reads::Params {
                paired,
                long,
                read_length,
                mu,
                sigma,
                min_len,
                max_len,
                circular,
                include_ns,
                per_read_tsv,
                seq_id: &ids[ci],
                source_fasta: &source_fasta,
            };
            let mut spans = Vec::new();
            let mut done = 0;
            while done < budget {
                let n = GEN_CHUNK.min(budget - done);
                let Some(c) = gen_reads::gen_chunk(&mut mt, seq.as_bytes(), n, count, &p) else {
                    return Err(io::Error::new(io::ErrorKind::InvalidData,
                        format!("couldn't make reads from contig '{}' with these settings", ids[ci])));
                };
                out1.push(c.r1.as_bytes())?;
                if let Some(o) = out2.as_mut() {
                    o.push(c.r2.as_bytes())?;
                }
                if let Some(o) = out_tsv.as_mut() {
                    o.push(c.tsv.as_bytes())?;
                }
                spans.extend(c.spans);
                bases += c.realized_bases;
                written += if paired { 2 * n } else { n };
                count += n;
                done += n;
            }
            spans_by_contig.push((ci, gen_reads::merge(spans)));
        }
        Ok((count, written, bases, spans_by_contig, mt))
    });
    r1.borrow_mut().out = Some(out1);
    if let (Some(w), Some(o)) = (&r2, out2) {
        w.borrow_mut().out = Some(o);
    }
    if let (Some(w), Some(o)) = (&tsv, out_tsv) {
        w.borrow_mut().out = Some(o);
    }
    let (count, written, bases, spans, mt) = result.map_err(to_py_err)?;
    Ok((count, written, bases, spans, mt.state.to_vec(), mt.index, mt.gauss_next))
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(set_threads, m)?)?;
    m.add_function(wrap_pyfunction!(trim_read, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(collect_contexts, m)?)?;
    m.add_function(wrap_pyfunction!(foldback::scan_foldbacks, m)?)?;
    m.add_function(wrap_pyfunction!(foldback::trim_long_folds, m)?)?;
    m.add_function(wrap_pyfunction!(foldback::filter_foldback_files, m)?)?;
    m.add_class::<Matcher>()?;
    m.add_class::<ReadMatcher>()?;
    m.add_class::<GzWriter>()?;
    m.add_function(wrap_pyfunction!(bam_mapped_reads_pids, m)?)?;
    m.add_function(wrap_pyfunction!(bam_per_contig_pid_stats, m)?)?;
    m.add_function(wrap_pyfunction!(bam_coverage_and_pid_stats, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_records_summary, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_gc_windows, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_filter_by_length, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_remove_wraps, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_primer_sites, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_assembly_summary, m)?)?;
    m.add_function(wrap_pyfunction!(rng_sample, m)?)?;
    m.add_function(wrap_pyfunction!(primer_sites_in_seq, m)?)?;
    m.add_function(wrap_pyfunction!(gen_reads_contigs, m)?)?;
    m.add_class::<Progress>()?;
    Ok(())
}
