//! Cutting reads where adapter read-through starts (the rules are set out in trim_read's
//! docstring in bit/modules/find_overrep.py).

use memchr::memmem;

pub const TRIM_PROBE_LEN: usize = 10;
pub const TRIM_MIN_END_OVERLAP: usize = 6;
pub const POLY_G_RUN: usize = 8;

pub struct Adapter {
    pub name: String,
    pub seq: Vec<u8>,
    first: memmem::Finder<'static>,
    second: Option<memmem::Finder<'static>>,
}

pub struct Trimmer {
    pub adapters: Vec<Adapter>,
}

/// which adapter cut a read: an index into Trimmer::adapters, possibly more than one where
/// adapters share a start and the read ends too soon to tell them apart, or poly-G
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Cut {
    Adapters(Vec<usize>),
    PolyG,
}

fn prefix_match(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

impl Trimmer {
    pub fn new(adapters: &[(String, String)]) -> Self {
        let adapters = adapters
            .iter()
            .map(|(name, seq)| {
                let s = seq.as_bytes().to_vec();
                let first = memmem::Finder::new(&s[..s.len().min(TRIM_PROBE_LEN)]).into_owned();
                let second = (s.len() >= 2 * TRIM_PROBE_LEN).then(|| {
                    memmem::Finder::new(&s[TRIM_PROBE_LEN..2 * TRIM_PROBE_LEN]).into_owned()
                });
                Adapter {
                    name: name.clone(),
                    seq: s,
                    first,
                    second,
                }
            })
            .collect();
        Trimmer { adapters }
    }

    /// (length kept, what cut it, if anything)
    pub fn trim(&self, seq: &[u8]) -> (usize, Option<Cut>) {
        let n = seq.len();
        let mut cut = n;
        let mut which: Option<Vec<usize>> = None;
        let mut which_seq: &[u8] = b"";
        for (idx, a) in self.adapters.iter().enumerate() {
            let mut i: Option<usize> = a.first.find(seq);
            // the adapter's next stretch, in case of an error in its first
            if i.is_none()
                && let Some(second) = &a.second
                && let Some(j) = second.find(seq)
                && j >= TRIM_PROBE_LEN
            {
                i = Some(j - TRIM_PROBE_LEN);
            }
            if i.is_none() {
                let probe = &a.seq[..a.seq.len().min(TRIM_PROBE_LEN)];
                let top = (TRIM_PROBE_LEN - 1).min(n);
                for l in (TRIM_MIN_END_OVERLAP..=top).rev() {
                    let start = &probe[..l.min(probe.len())];
                    if seq.ends_with(start) {
                        i = n.checked_sub(l);
                        break;
                    }
                }
            }
            let Some(i) = i else { continue };
            if i < cut {
                cut = i;
                which = Some(vec![idx]);
                which_seq = &a.seq;
            } else if i == cut
                && let Some(w) = which.as_mut()
            {
                // adapters can share a start (TruSeq R1 and R2 do): the longer match wins,
                // and an even one is reported as either
                let mine = prefix_match(&seq[i..], &a.seq);
                let theirs = prefix_match(&seq[i..], which_seq);
                if mine > theirs {
                    *w = vec![idx];
                    which_seq = &a.seq;
                } else if mine == theirs {
                    w.push(idx);
                }
            }
        }
        let g = seq.iter().rev().take_while(|&&b| b == b'G').count();
        if g >= POLY_G_RUN && n - g < cut {
            return (n - g, Some(Cut::PolyG));
        }
        (cut, which.map(Cut::Adapters))
    }

    pub fn cut_name(&self, cut: &Cut) -> String {
        match cut {
            Cut::PolyG => "poly-G".to_string(),
            Cut::Adapters(idx) => idx
                .iter()
                .map(|&i| self.adapters[i].name.as_str())
                .collect::<Vec<_>>()
                .join("|"),
        }
    }
}
