//! A count-min sketch of how many reads (or pairs) hold each anchor, for the counting pass
//! of `bit normalize-reads`: fixed memory whatever the input, filled from many threads at once.
//!
//! Blocks of 32 16-bit counters (saturating at 65,535), each block one 64-byte cache line, in
//! two lanes of 16. Each anchor has two blocks, picked by remixing its hash, and one counter in
//! each lane of them, picked by further bits: ROWS counters in two memory accesses instead of
//! one per row (blocked count-min). One block for all four would take one access, but anchors
//! sharing a block then share their counters together, and on Zymo Log about 3 times as many
//! counts came out too high, enough to pass more error reads as clean. Each read holding an
//! anchor adds one to all its counters; its count is the smallest of them, which can only be
//! too high (when every one of its counters is shared with other anchors), never too low.
//! Counters are atomics raised without locks; as every raise is kept (adds don't depend on each
//! other, unlike conservative update's), the counts come out the same whatever the threads or
//! their timing, and so does what bit normalize-reads keeps.

use std::sync::atomic::{AtomicU16, Ordering};

pub const ROWS: usize = 4;
const LANE: usize = 16;
const BLOCK: usize = 2 * LANE;
/// blocks per anchor
const SPREAD: usize = ROWS * LANE / BLOCK;
/// how many anchors ahead add_all() and get_all() prefetch
const AHEAD: usize = 8;
const BLOCK_SEEDS: [u64; SPREAD] = [0x9E3779B97F4A7C15, 0xD6E8FEB86659FD93];
const LANE_SEED: u64 = 0xC2B2AE3D27D4EB4F;

#[repr(align(64))]
struct Block([AtomicU16; BLOCK]);

#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

pub struct Sketch {
    blocks: Vec<Block>,
}

impl Sketch {
    /// a sketch taking about `bytes` of memory (at least 4,096 counters)
    pub fn new(bytes: usize) -> Self {
        let n = (bytes / (2 * BLOCK)).max(4096 / BLOCK);
        let mut blocks = Vec::with_capacity(n);
        blocks.resize_with(n, || Block(std::array::from_fn(|_| AtomicU16::new(0))));
        Sketch { blocks }
    }

    /// `h`'s blocks
    #[inline(always)]
    fn blocks_of(&self, h: u64) -> [&Block; SPREAD] {
        BLOCK_SEEDS.map(|seed| &self.blocks[((mix(h ^ seed) as u128 * self.blocks.len() as u128) >> 64) as usize])
    }

    /// `h`'s counters: one in each lane of its blocks
    #[inline(always)]
    fn counters(&self, h: u64) -> [&AtomicU16; ROWS] {
        let blocks = self.blocks_of(h);
        let x = mix(h ^ LANE_SEED);
        std::array::from_fn(|r| &blocks[r / 2].0[(r % 2) * LANE + ((x >> (4 * r)) as usize & (LANE - 1))])
    }

    /// asks for `h`'s counters to be brought into cache ahead of add() or get() (a hint, so the
    /// misses of the anchors coming up overlap instead of waiting one at a time; nothing on
    /// other CPUs)
    #[inline(always)]
    pub fn prefetch(&self, h: u64) {
        for block in self.blocks_of(h) {
            let p = block as *const Block as *const i8;
            #[cfg(target_arch = "x86_64")]
            unsafe {
                std::arch::x86_64::_mm_prefetch(p, std::arch::x86_64::_MM_HINT_T0)
            };
            #[cfg(target_arch = "aarch64")]
            unsafe {
                std::arch::asm!("prfm pldl1keep, [{0}]", in(reg) p, options(nostack, preserves_flags, readonly))
            };
            #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
            let _ = p;
        }
    }

    /// counts one more read holding each of `hs`
    pub fn add_all(&self, hs: &[u64]) {
        for (i, &h) in hs.iter().enumerate() {
            if let Some(&next) = hs.get(i + AHEAD) {
                self.prefetch(next);
            }
            self.add(h);
        }
    }

    /// counts of each of `hs`, into `out`
    pub fn get_all(&self, hs: impl Iterator<Item = u64> + Clone, out: &mut Vec<u32>) {
        out.clear();
        let mut ahead = hs.clone().skip(AHEAD);
        for h in hs {
            if let Some(next) = ahead.next() {
                self.prefetch(next);
            }
            out.push(self.get(h));
        }
    }

    /// counts one more read holding `h`
    #[inline]
    pub fn add(&self, h: u64) {
        for c in self.counters(h) {
            let _ = c.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| {
                (v < u16::MAX).then(|| v + 1)
            });
        }
    }

    /// how many reads hold `h` (possibly too high, see above)
    #[inline]
    pub fn get(&self, h: u64) -> u32 {
        self.counters(h)
            .iter()
            .map(|c| c.load(Ordering::Relaxed))
            .min()
            .unwrap() as u32
    }

    pub fn bytes(&self) -> usize {
        self.blocks.len() * BLOCK * 2
    }

    /// share of counters not at zero; past about 85%, counts of rare anchors run high often
    /// enough to matter (a larger sketch, --memory, brings it down)
    pub fn fill(&self) -> f64 {
        let used = self
            .blocks
            .iter()
            .flat_map(|b| b.0.iter())
            .filter(|c| c.load(Ordering::Relaxed) > 0)
            .count();
        used as f64 / (self.blocks.len() * BLOCK) as f64
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
        let s = Sketch::new(4096); // 4,096 counters for 5,000 anchors
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
