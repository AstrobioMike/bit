//! bit's compiled core: the per-read loops of `bit find-overrep-seqs` (reading, adapter trimming,
//! k-mer counting, gathering reads around seeds, and placing found sequences in reads
//! allowing for edits), of `bit filter-foldbacks` (finding fold-back reads and pairs, short or long),
//! and of `bit normalize-reads` (counting anchors into a sketch, judging reads and pairs, and writing
//! the kept ones).
//! Everything else stays in Python (bit/modules/find_overrep.py,
//! bit/modules/filter_foldbacks.py, bit/modules/normalize_reads.py), which calls into this as `bit._native`.

mod anchors;
mod edlib;
mod foldback;
mod kmers;
mod matcher;
mod normrun;
mod reads;
mod sketch;
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
    #[pyo3(signature = (seqs, max_diff, min_read_len, groups=None))]
    fn new(
        seqs: Vec<String>,
        max_diff: f64,
        min_read_len: usize,
        groups: Option<Vec<usize>>,
    ) -> PyResult<Self> {
        if let Some(g) = &groups
            && g.len() != seqs.len()
        {
            return Err(PyValueError::new_err("groups must be as long as seqs"));
        }
        let seqs: Vec<Vec<u8>> = seqs.into_iter().map(|s| s.into_bytes()).collect();
        Ok(Matcher {
            inner: matcher::Matcher::new(&seqs, max_diff, min_read_len, groups.as_deref()),
        })
    }

    /// {seq index: (start in read or None if the read is within it, is reverse complement,
    /// edits)}; `prefilter=False` tries every sequence (for checking the prefilter)
    #[pyo3(signature = (read, prefilter=true))]
    fn hits<'py>(
        &self,
        py: Python<'py>,
        read: &str,
        prefilter: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for h in self.inner.hits(read.as_bytes(), !prefilter) {
            d.set_item(h.seq, (h.start, h.is_rc, h.edits))?;
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

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(set_threads, m)?)?;
    m.add_function(wrap_pyfunction!(trim_read, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(collect_contexts, m)?)?;
    m.add_function(wrap_pyfunction!(foldback::scan_foldbacks, m)?)?;
    m.add_function(wrap_pyfunction!(foldback::trim_long_folds, m)?)?;
    m.add_class::<Matcher>()?;
    m.add_function(wrap_pyfunction!(normrun::count_anchors, m)?)?;
    m.add_class::<normrun::AnchorCounts>()?;
    m.add_class::<Progress>()?;
    Ok(())
}
