//! Cutting reads where adapter read-through starts (the rules are set out in trim_read's
//! docstring in bit/modules/find_overrep.py).

use memchr::memmem;

pub const TRIM_PROBE_LEN: usize = 10;
pub const TRIM_PROBE_SPAN: usize = 40;
pub const TRIM_MIN_END_OVERLAP: usize = 6;
pub const POLY_G_RUN: usize = 8;

pub struct Adapter {
    pub name: String,
    pub seq: Vec<u8>,
    /// stretches of the adapter looked for in reads, with where each starts in the adapter
    probes: Vec<(usize, memmem::Finder<'static>)>,
}

/// where the stretches looked for start in an adapter `len` long: every TRIM_PROBE_LEN bases
/// over its first TRIM_PROBE_SPAN, and the last full stretch there (all of it, if it's shorter
/// than one)
pub fn probe_offsets(len: usize) -> Vec<usize> {
    let span = len.min(TRIM_PROBE_SPAN);
    if span <= TRIM_PROBE_LEN {
        return vec![0];
    }
    let mut out: Vec<usize> = (0..=span - TRIM_PROBE_LEN)
        .step_by(TRIM_PROBE_LEN)
        .collect();
    if *out.last().unwrap() != span - TRIM_PROBE_LEN {
        out.push(span - TRIM_PROBE_LEN);
    }
    out
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
                let probes = probe_offsets(s.len())
                    .into_iter()
                    .map(|o| {
                        let end = (o + TRIM_PROBE_LEN).min(s.len());
                        (o, memmem::Finder::new(&s[o..end]).into_owned())
                    })
                    .collect();
                Adapter {
                    name: name.clone(),
                    seq: s,
                    probes,
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
            // its first stretch, else a later one (in case of an error in the first, or a read
            // that starts partway into the adapter, which is cut from its start then)
            let mut i: Option<usize> = a
                .probes
                .iter()
                .find_map(|(o, f)| f.find(seq).map(|j| j.saturating_sub(*o)));
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offsets() {
        assert_eq!(probe_offsets(33), vec![0, 10, 20, 23]);
        assert_eq!(probe_offsets(19), vec![0, 9]);
        assert_eq!(probe_offsets(20), vec![0, 10]);
        assert_eq!(probe_offsets(8), vec![0]);
        assert_eq!(probe_offsets(300), vec![0, 10, 20, 30]);
    }
}
