//! A count-min sketch of how many reads (or pairs) hold each anchor, for the counting pass
//! of `bit normalize-reads`: fixed memory whatever the input, filled from many threads at once.
//!
//! ROWS rows of 16-bit counters (saturating at 65,535). Each anchor has one counter per row,
//! picked by remixing its hash with a per-row seed, and each read holding it adds one to all
//! of them; its count is the smallest of them, which can only be too high (when every one of
//! its counters is shared with other anchors), never too low. Counters are atomics raised
//! without locks; as every raise is kept (adds don't depend on each other, unlike conservative
//! update's), the counts come out the same whatever the threads or their timing, and so does
//! what bit normalize-reads keeps.

use std::sync::atomic::{AtomicU16, Ordering};

pub const ROWS: usize = 4;
const SEEDS: [u64; ROWS] = [
    0x9E3779B97F4A7C15,
    0xC2B2AE3D27D4EB4F,
    0x165667B19E3779F9,
    0x27D4EB2F165667C5,
];

#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

pub struct Sketch {
    cells: Vec<AtomicU16>,
    width: usize,
}

impl Sketch {
    /// a sketch taking about `bytes` of memory (at least 1,024 counters a row)
    pub fn new(bytes: usize) -> Self {
        let width = (bytes / (2 * ROWS)).max(1024);
        let mut cells = Vec::with_capacity(width * ROWS);
        cells.resize_with(width * ROWS, || AtomicU16::new(0));
        Sketch { cells, width }
    }

    #[inline(always)]
    fn slots(&self, h: u64) -> [usize; ROWS] {
        let mut out = [0; ROWS];
        for (r, o) in out.iter_mut().enumerate() {
            let x = mix(h ^ SEEDS[r]);
            *o = r * self.width + ((x as u128 * self.width as u128) >> 64) as usize;
        }
        out
    }

    /// counts one more read holding `h`
    #[inline]
    pub fn add(&self, h: u64) {
        for s in self.slots(h) {
            let _ = self.cells[s].fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| {
                (v < u16::MAX).then_some(v + 1)
            });
        }
    }

    /// how many reads hold `h` (possibly too high, see above)
    #[inline]
    pub fn get(&self, h: u64) -> u32 {
        self.slots(h)
            .iter()
            .map(|&s| self.cells[s].load(Ordering::Relaxed))
            .min()
            .unwrap() as u32
    }

    pub fn bytes(&self) -> usize {
        self.cells.len() * 2
    }

    /// share of counters not at zero; past about half, counts of rare anchors start to run
    /// high often enough to matter (a larger sketch, --memory, brings it down)
    pub fn fill(&self) -> f64 {
        let used = self
            .cells
            .iter()
            .filter(|c| c.load(Ordering::Relaxed) > 0)
            .count();
        used as f64 / self.cells.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_exactly_when_roomy() {
        let s = Sketch::new(1 << 20);
        for h in 0..1000u64 {
            for _ in 0..(h % 7) {
                s.add(h.wrapping_mul(0x9E3779B97F4A7C15));
            }
        }
        for h in 0..1000u64 {
            assert_eq!(s.get(h.wrapping_mul(0x9E3779B97F4A7C15)), (h % 7) as u32);
        }
    }

    #[test]
    fn never_counts_low_when_crowded() {
        let s = Sketch::new(4096); // 1,024 counters a row for 5,000 anchors
        for h in 0..5000u64 {
            for _ in 0..3 {
                s.add(mix(h));
            }
        }
        assert!((0..5000u64).all(|h| s.get(mix(h)) >= 3));
    }

    #[test]
    fn saturates() {
        let s = Sketch::new(1 << 12);
        for _ in 0..70_000 {
            s.add(42);
        }
        assert_eq!(s.get(42), u16::MAX as u32);
    }
}
