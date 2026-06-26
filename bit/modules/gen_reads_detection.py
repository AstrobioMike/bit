"""
Detection tracking for read generation

Detection is computed incrementally as reads are produced, so gen-reads never
has to store read spans or post-process an alignment to report it.

The tracker keeps the *complement* of what's covered: per contig, a sorted list
of not-yet-covered gaps. Each contig starts as a single gap spanning its whole
length; every read span punches holes in the gaps it overlaps. A contig is fully
covered when it has no gaps left, and the genome reaches detection 1.0 when no
contig has any gap remaining. Because the gap lists only ever shrink, memory is
bounded by the number of remaining uncovered segments (which trends to zero as
coverage fills in) rather than by the number of reads, and no read span is ever
stored or sorted.

Two running counters avoid any rescan: `_uncovered[ci]` is the uncovered base
count for one contig, and `_remaining` is their sum across the genome. Detection
is then (total - _remaining) / total, an O(1) read at any time.

Note on reachability: with fixed-length reads whose start positions are uniform,
the final base(s) of a contig are coverable by very few read placements, so exact
detection 1.0 is often not reached even at high coverage (a handful of terminal
bases linger). Callers that present detection should round for display rather than
relying on the value being exactly 1.0; the `done` flag is an exact "every base
covered" signal kept as a cheap early-out, due to the ends though that won't always happen.
"""

import bisect


class DetectionTracker:
    """
    Track the fraction of a genome's bases covered by >= 1 read, updated online.

    Construct with the genome's contig lengths in the same order spans will be
    reported against (spans are keyed by contig index). Feed read spans with
    `add(contig_index, start, end)` using 0-based, end-exclusive coordinates
    within that contig. Read `detection()` for the current covered fraction, or
    `done` for the exact all-covered flag.
    """

    def __init__(self, contig_lengths):
        self._lengths = list(contig_lengths)
        self._total = sum(self._lengths)

        # per-contig parallel arrays of gap starts/ends (sorted, disjoint).
        # a zero-length contig contributes no gaps and no uncovered bases.
        self._gap_starts = []
        self._gap_ends = []
        self._uncovered = []
        self._live = set()
        for ci, length in enumerate(self._lengths):
            if length > 0:
                self._gap_starts.append([0])
                self._gap_ends.append([length])
                self._uncovered.append(length)
                self._live.add(ci)
            else:
                self._gap_starts.append([])
                self._gap_ends.append([])
                self._uncovered.append(0)

        self._remaining = self._total
        self.done = self._remaining == 0

    def add(self, contig_index, start, end):
        """
        Record that read span [start, end) covers part of contig_index.

        Coordinates are 0-based and end-exclusive. Spans wholly inside already
        covered regions are no-ops. Calls after the genome is fully covered, or
        against an already-complete contig, return in O(1).
        """
        if self.done or contig_index not in self._live:
            return

        gap_starts = self._gap_starts[contig_index]
        gap_ends = self._gap_ends[contig_index]

        # leftmost gap that could overlap the span is the first whose end > start;
        # gaps are disjoint and sorted, so everything earlier ends at/before start.
        lo = bisect.bisect_right(gap_ends, start)
        if lo >= len(gap_starts):
            return  # span lies entirely past the last remaining gap

        i = lo
        removed = 0
        new_starts = []
        new_ends = []
        # walk the gaps the span reaches (those starting before `end`), replacing
        # each with whatever remnants survive the punch.
        while i < len(gap_starts) and gap_starts[i] < end:
            g_s = gap_starts[i]
            g_e = gap_ends[i]
            overlap_s = start if start > g_s else g_s
            overlap_e = end if end < g_e else g_e
            if overlap_s < overlap_e:
                removed += overlap_e - overlap_s
                if g_s < overlap_s:          # left remnant survives
                    new_starts.append(g_s)
                    new_ends.append(overlap_s)
                if overlap_e < g_e:          # right remnant survives
                    new_starts.append(overlap_e)
                    new_ends.append(g_e)
            else:
                # no real overlap; keep the gap as-is (defensive, shouldn't occur
                # given the bisect, but preserves correctness if it ever does)
                new_starts.append(g_s)
                new_ends.append(g_e)
            i += 1

        if removed:
            # splice the surviving remnants in place of the consumed gaps [lo:i)
            gap_starts[lo:i] = new_starts
            gap_ends[lo:i] = new_ends
            self._uncovered[contig_index] -= removed
            self._remaining -= removed
            if self._uncovered[contig_index] <= 0:
                self._live.discard(contig_index)
                if not self._live:
                    self.done = True

    def detection(self):
        """ Current covered fraction in [0, 1]. Empty genome is defined as 1.0. """
        if self._total == 0:
            return 1.0
        return (self._total - self._remaining) / self._total

    def covered_bases(self):
        """ Number of bases covered by >= 1 read so far. """
        return self._total - self._remaining

    def genome_size(self):
        """ Total bases across all contigs (the detection denominator). """
        return self._total
