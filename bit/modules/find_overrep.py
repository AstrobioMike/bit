"""
Core logic for `bit find-overrep-seqs`.

Finds and counts overrepresented sequences in reads without knowing what they are ahead of
time (e.g., primers in a short run of an unamplified library), allowing for sequencing errors.
Intended for short reads (<= ~300 bp): the thresholds are fractions of reads, so on long reads
ordinary sequence passes them; not currently suitable for long reads.
Reads (R1 and R2, if paired) are pooled for finding things, since a single-stranded library
prep puts them in reads in either orientation, and counts are reported per read file.

  1. Adapters. Every k-mer (default k=11) is counted, on both strands combined (canonical).
     At k=11 there are only ~4.2M possible k-mers, so counts are a plain array, no hashing.
     Adapter read-through is trimmed from reads from here on (the built-in Illumina adapters
     and any given), since otherwise every insert's end next to an adapter makes for an
     overrepresented k-mer. With paired reads, unknown adapters are found too: an insert is
     read forward in one read of a pair and reverse complemented in the other, so its k-mers
     are about as common forward in R1 as reverse complemented in R2, while read-through is
     in only one read, and only forward (see mate_asymmetry).

  2. Candidates. K-mers are counted again in the trimmed reads, along with the 2-mers and
     3-mers of an order-2 Markov model of the reads' own composition, as the background. A
     k-mer is a candidate if it's in at least `min_perc` percent of reads and `min_fold` times its
     expected count, or in at least 1% of reads whatever its expected count (the model can
     expect a lot of what it predicts well, like a dominant low-complexity repeat). Rare
     one-off error variants of abundant k-mers are dropped.

  3. Extension. From the most abundant candidate down, each is extended a base at a time
     in both directions while most occurrences of the current k-mer carry on to the same
     next one (a consensus, as in greedy assemblers), skipping candidates already inside an
     extended sequence. Small k-mers recur by chance (e.g., in a genome), which can muddy
     the next base, so the ends are then extended further using the reads that hold them,
     dropping reads that disagree as it goes (see extend_by_reads). Sequences that are just
     a variant of a more abundant one (e.g., missing a base or few at an end) are collapsed.

  4. Counting. A read holds a found sequence if it aligns within the (trimmed) read with at
     most `max_diff` edits per base (substitutions + indels), or, if the read holds none in
     full, the read aligns within it (e.g., a primer missing a few bases, or a read from
     something longer than reads). Start positions in reads are tallied too. Sequences too
     rare once counted, pieces of another found from a seed that stopped short, and rare
     near-copies of a far more common one are left out of what's reported. Each is reported
     in the orientation it's most often in, in the read file it's most common in.

The per-read loops (reading, trimming, k-mer counting, gathering reads around seeds, and
placing found sequences in reads) run in bit's compiled core (rust/, as bit._native);
deciding what's found stays here.
"""

import os
import sys
import math
import threading
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np # type: ignore
import edlib # type: ignore

from bit import _native as native
from bit.modules.general import (check_if_output_dir_exists,
                                 linguistic_complexity,
                                 log_command_run)


DEFAULT_K = 11
DEFAULT_MIN_PERC = 0.1      # percent of reads
DEFAULT_MIN_FOLD = 10
DEFAULT_MAX_DIFF = 0.1
DEFAULT_MAX_LEN = 300
DEFAULT_MAX_READS_FOR_DISCOVERY = 10_000_000   # per file; 0 for all

# finding candidates
DOMINANT_FRAC = 0.01         # a k-mer in at least this share of reads is a candidate whatever its fold
ERROR_VARIANT_RATIO = 20     # a k-mer this many times rarer than one a substitution away is taken for an error

# adapters
ADAPTER_SEARCH_SEEDS = 200   # most abundant candidates extended to look for unknown adapters
MAX_INSERT_ASYMMETRY = 0.5   # k-mers at most this asymmetric between the reads of a pair are insert
ADAPTER_PROBE_LEN = 16       # start of a found adapter compared to others, to tell if it's the same
TRIM_PROBE_LEN = 10          # stretches of an adapter looked for in reads to trim read-through...
TRIM_PROBE_SPAN = 40         # ...every TRIM_PROBE_LEN bases over this much of its start (see probe_offsets)
TRIM_MIN_END_OVERLAP = 6     # shortest adapter start trimmed when it's cut off by a read end
POLY_G_RUN = 8               # runs of G this long at a read end are trimmed (no-signal on 2-color instruments)

# extension
EXTEND_AGREE = 0.8           # share of the current k-mer's occurrences that must carry on to the same next one
EXTEND_MIN_SHARE = 0.05      # stop where the next k-mers' total drops below this share of the seed's count...
EXTEND_MIN_COUNT = 5         # ...or below this many
MAX_CONTEXTS = 2000          # reads gathered per sequence end when extending by reads
REFINE_MIN_COVER = 0.5       # extending by reads stops where fewer than this share of reads still in agreement carry on,
REFINE_MIN_SHARE = 0.6       # ...or the next base is less than this share of them,
REFINE_MIN_RATIO = 3         # ...or less than this many times the next most common

# collapsing and reporting
LOW_COMPLEXITY = 0.4         # see low_complexity(); the same cutoff cov-analyzer uses
MIN_INSERT_LEN = 12          # found sequences shorter than this (after trimming) aren't reported
VARIANT_SLACK = 0.25         # see is_variant()
CONTAINED_SHARE = 0.8        # see the end of find_overrep()
ERROR_VARIANT_SHARE = 0.05   # see the end of find_overrep()
MIN_PARTIAL_LEN = 15         # a (trimmed) read shorter than this isn't checked for being within a found sequence

# given sequences to count
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
MAX_VARIANTS = 4096          # most variants a given degenerate sequence may expand to
COMPARE_VARIANTS = 64        # variants of a given sequence used when comparing it to others

BUILTIN_ADAPTERS = {
    "illumina_truseq_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "illumina_truseq_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "illumina_nextera": "CTGTCTCTTATACACATCT",
    "illumina_pcr_free_me_rc": "ATGTGTATAAGAGACA",
    "illumina_small_rna_3p": "TGGAATTCTCGGGTGCCAAGG",
}

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_ENC = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _ENC[ord(_b)] = _i
    _ENC[ord(_b.lower())] = _i


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


### k-mer encoding ###

def encode(seqs):
    """
    The sequences as one array of 2-bit base codes (4 for anything but A/C/G/T), each followed
    by a separator (also 4), with each sequence's start in it.
    """
    lengths = np.fromiter((len(s) for s in seqs), dtype=np.int64, count=len(seqs))
    joined = ("\n".join(seqs) + "\n").encode()
    enc = _ENC[np.frombuffer(joined, dtype=np.uint8)]
    starts = np.zeros(len(seqs), dtype=np.int64)
    if len(seqs) > 1:
        starts[1:] = np.cumsum(lengths + 1)[:-1]
    return enc, starts, lengths


def window_codes(enc, k):
    """(code, valid) of every k-length window of `enc`; valid if it's all A/C/G/T."""
    m = len(enc) - k + 1
    if m <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=bool)
    codes = np.zeros(m, dtype=np.int64)
    for j in range(k):
        codes = (codes << 2) | (enc[j:j + m] & 3).astype(np.int64)
    bad = np.concatenate([[0], np.cumsum(enc > 3)])
    valid = (bad[k:k + m] - bad[:m]) == 0
    return codes, valid


def rc_codes(codes, k):
    out = np.zeros_like(codes)
    x = codes.copy()
    for _ in range(k):
        out = (out << 2) | (3 - (x & 3))
        x >>= 2
    return out


def code_to_seq(code, k):
    return "".join("ACGT"[(code >> (2 * (k - 1 - i))) & 3] for i in range(k))


def seq_to_code(seq):
    code = 0
    for b in seq:
        code = (code << 2) | "ACGT".index(b)
    return code


### progress ###

def quietly(desc, fn, *args, **kwargs):
    """runs a pass over the reads with no progress shown (`desc` unused)"""
    return fn(*args, **kwargs)


def with_progress_bar(read_paths):
    """
    A runner for passes over the reads (like `quietly`) that shows `desc` and a progress bar
    over the read files' bytes, with reads so far, as match-reads does. The pass runs in a
    thread of its own (the compiled core lets go of the GIL) while this one draws the bar.
    """
    total = sum(os.path.getsize(p) for p in read_paths)

    def run(desc, fn, *args, **kwargs):  # pragma: no cover (drawing only)
        prog = native.Progress()
        result, error = [None], [None]

        def work():
            try:
                result[0] = fn(*args, progress=prog, **kwargs)
            except BaseException as e:  # passed on below
                error[0] = e

        from tqdm import tqdm # type: ignore
        print(f"\n    {desc}\n")
        thread = threading.Thread(target=work, daemon=True)
        thread.start()
        with tqdm(total=total, unit="B", unit_scale=True, unit_divisor=1024, ncols=80, file=sys.stdout,
                  bar_format="      {l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}{postfix}]") as bar:
            while thread.is_alive():
                thread.join(0.2)
                bar.update(min(prog.bytes, total) - bar.n)
                bar.set_postfix_str(f"{prog.reads:,} reads", refresh=False)
            if error[0] is None:
                # (a pass can stop early once it has what it needs; and one quicker than a
                # refresh would otherwise end without its read count)
                bar.set_postfix_str(f"{prog.reads:,} reads", refresh=False)
                bar.update(total - bar.n)
                bar.refresh()
        print()
        if error[0] is not None:
            raise error[0]
        return result[0]

    return run


### pass 1: counting ###

class KmerCounts:
    """
    Forward-strand counts of k-mers (per read file and combined), 3-mers and 2-mers, the
    valid k-mer windows at each read position, and how many reads, over all reads given,
    cut where adapter read-through starts if `adapters` are given.
    """

    def __init__(self, read_paths, k, adapters=None, run=quietly, desc="", max_reads=0):
        self.k = k
        by_file, tri, di, windows_at, self.n_reads = run(
            desc, native.count_kmers, [str(p) for p in read_paths], k, _adapter_list(adapters),
            max_reads=max_reads)
        # signed, as the arithmetic on them expects (e.g., negating to sort)
        by_file = [a.astype(np.int64) for a in by_file]
        self.tri, self.di, self.windows_at = (a.astype(np.int64) for a in (tri, di, windows_at))
        self.by_mate = dict(enumerate(by_file))
        self.kmers = sum(by_file)

    def canonical(self):
        """
        (counts on both strands combined, indexed by either strand's code; the reverse
        complement of every code; whether each code is its own reverse complement)
        """
        all_codes = np.arange(4 ** self.k, dtype=np.int64)
        rc = rc_codes(all_codes, self.k)
        pal = rc == all_codes
        both = self.kmers + self.kmers[rc]
        both[pal] = self.kmers[pal]
        return both, rc, pal

    def expected_probs(self, rc, pal):
        """
        Probability of each k-mer (either strand) under an order-2 Markov model of the reads'
        own composition, so ordinary sequence, even skewed sequence, isn't flagged.
        """
        k = self.k
        p3 = (self.tri + 1) / (self.tri.sum() + 64)
        p2 = (self.di + 1) / (self.di.sum() + 16)
        logp3, logp2 = np.log(p3), np.log(p2)
        codes = np.arange(4 ** k, dtype=np.int64)
        logp = logp3[(codes >> (2 * (k - 3))) & 63]
        for i in range(1, k - 2):
            tri = (codes >> (2 * (k - 3 - i))) & 63
            logp += logp3[tri] - logp2[tri >> 2]
        prob = np.exp(logp)
        both = prob + prob[rc]
        both[pal] = prob[pal]
        return both


### finding candidates ###

def find_candidates(counts, min_frac, min_fold, error_ref=None):
    """
    Canonical codes of candidate k-mers (see module notes), most abundant first, and the
    canonical counts. Error variants are judged against `error_ref` (canonical counts, e.g.,
    from untrimmed reads) if given.
    """
    k = counts.k
    both, rc, pal = counts.canonical()
    probs = counts.expected_probs(rc, pal)
    codes = np.arange(4 ** k, dtype=np.int64)
    canonical = codes <= rc
    total_windows = counts.windows_at.sum()
    frac = both / max(counts.n_reads, 1)

    picked = canonical & (frac >= min_frac) & (both >= min_fold * probs * total_windows)
    # and anything in a large share of reads, whatever the background model expects (it can
    # expect a lot of a sequence it predicts well, like a dominant low-complexity repeat)
    picked |= canonical & (frac >= max(min_frac, DOMINANT_FRAC))
    picked = codes[picked]

    if len(picked):
        picked = not_error_variants(picked, both if error_ref is None else error_ref, k)
    order = sorted(picked.tolist(), key=lambda c: (-both[c], c))
    return order, both


def not_error_variants(codes, both, k, ratio=None):
    """
    The codes left after dropping likely sequencing-error variants: those with a k-mer one
    substitution away at least `ratio` times as abundant. An abundant sequence (say a primer
    in a third of reads) otherwise makes every one-off error k-mer of it look overrepresented.
    """
    ratio = ERROR_VARIANT_RATIO if ratio is None else ratio
    drop = np.zeros(len(codes), dtype=bool)
    for i in range(k):
        shift = 2 * (k - 1 - i)
        base = (codes >> shift) & 3
        for d in range(1, 4):
            other = codes ^ ((base ^ ((base + d) & 3)) << shift)
            drop |= both[other] >= ratio * both[codes]
    return codes[~drop]


### extension ###

def extend(seed, both, k, max_len=DEFAULT_MAX_LEN, agree=EXTEND_AGREE, min_share=EXTEND_MIN_SHARE):
    """
    Extends a k-mer (forward code) base by base in both directions while most occurrences of
    the current k-mer carry on to the same next one (see module notes). Returns (sequence, why the left end stopped, why the
    right end stopped), each "branch", "end" (nothing more in the reads), "repeat" (back on a
    k-mer already in it), or "length".
    """
    mask = 4 ** k - 1
    floor = max(EXTEND_MIN_COUNT, min_share * both[seed])
    seen = {seed}
    seq = code_to_seq(seed, k)

    def step(cur, nexts):
        n = both[nexts]
        total = n.sum()
        best = int(np.argmax(n))
        # where most reads holding this k-mer end (or are trimmed) there's nothing to agree
        # on, even if a minority carry on (e.g., a primer, and a dimer of it)
        if total < floor or n[best] < agree * both[cur]:
            return None, "end" if total < (1 - agree) * both[cur] or total < floor else "branch"
        if n[best] < agree * total:
            return None, "branch"
        return int(nexts[best]), None

    right, cur, right_stop = [], seed, "length"
    while len(seq) + len(right) < max_len:
        nxt, why = step(cur, np.array([((cur << 2) | b) & mask for b in range(4)]))
        if nxt is None:
            right_stop = why
            break
        if nxt in seen:
            right_stop = "repeat"
            break
        seen.add(nxt)
        right.append("ACGT"[nxt & 3])
        cur = nxt

    left, cur, left_stop = [], seed, "length"
    while len(seq) + len(right) + len(left) < max_len:
        nxt, why = step(cur, np.array([(b << (2 * (k - 1))) | (cur >> 2) for b in range(4)]))
        if nxt is None:
            left_stop = why
            break
        if nxt in seen:
            left_stop = "repeat"
            break
        seen.add(nxt)
        left.append("ACGT"[nxt >> (2 * (k - 1))])
        cur = nxt

    return "".join(reversed(left)) + seq + "".join(right), left_stop, right_stop


def canonical_kmers_of(seq, k, rc_all):
    enc, _, _ = encode([seq])
    codes, valid = window_codes(enc, k)
    codes = codes[valid]
    return np.minimum(codes, rc_all[codes])


def extend_all(candidates, both, rc_all, k, max_len=DEFAULT_MAX_LEN):
    """
    Extends candidates in order, skipping any already inside an extended sequence. Returns
    [dict(consensus, seed_count, left_stop, right_stop)].
    """
    used = np.zeros(4 ** k, dtype=bool)
    out = []
    for code in candidates:
        if used[code]:
            continue
        seq, left_stop, right_stop = extend(code, both, k, max_len)
        used[canonical_kmers_of(seq, k, rc_all)] = True
        out.append(dict(consensus=seq, seed_count=int(both[code]), left_stop=left_stop,
                        right_stop=right_stop))
    return out


def extend_by_reads(contexts, k, max_len):
    """
    Bases following the seed in its reads (oriented so the seed reads forward), one column at
    a time, while most reads still in agreement carry on and mostly agree on the next base.
    Reads disagreeing with a chosen base are dropped from then on, so reads that only share
    the seed k-mer by chance (common at small k in a big genome) fall away quickly.
    """
    active = list(contexts)
    out = []
    while active and len(out) < max_len:
        j = len(out)
        col = Counter(seq[p + k + j] for seq, p in active if p + k + j < len(seq))
        covering = sum(col.values())
        if covering < max(EXTEND_MIN_COUNT, REFINE_MIN_COVER * len(active)):
            break
        (base, n), *rest = col.most_common(2) + [(None, 0)]
        second = rest[0][1]
        if n < REFINE_MIN_SHARE * covering or n < REFINE_MIN_RATIO * second:
            break
        out.append(base)
        active = [(seq, p) for seq, p in active if p + k + j < len(seq) and seq[p + k + j] == base]
    return "".join(out)


def refine_ends(extended, read_paths, adapters, k, max_len=DEFAULT_MAX_LEN, run=quietly, desc="",
                max_reads=0):
    """
    Extends each sequence's ends further where the reads holding its end k-mers support it
    (see extend_by_reads). K-mer counts stop where a small k-mer's other occurrences (e.g.,
    in a genome) muddy the next base; the reads themselves don't.
    """
    seeds = []
    for f in extended:
        cons = f["consensus"]
        seeds.append(cons[-k:])
        seeds.append(revcomp(cons[:k]))
    contexts = run(desc, native.collect_contexts, [str(p) for p in read_paths], k, seeds,
                   MAX_CONTEXTS, _adapter_list(adapters), max_reads=max_reads)
    for i, f in enumerate(extended):
        cons = f["consensus"]
        room = max(0, max_len - len(cons))
        right = extend_by_reads(contexts[2 * i], k, room)
        left = extend_by_reads(contexts[2 * i + 1], k, max(0, room - len(right)))
        f["consensus"] = revcomp(left) + cons + right
    return extended


### adapters ###

def _locate(query, target, max_edits):
    """[(start, end inclusive)] of the best placements of `query` within `target`, if within max_edits"""
    if not query or len(target) < len(query) - max_edits:
        return []
    r = edlib.align(query, target, mode="HW", task="locations", k=max_edits)
    return r["locations"] if r["editDistance"] >= 0 else []


def _near(a, b):
    """same length and at most 1 mismatch in 10"""
    return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= len(a) // 10


def probe_offsets(length):
    """
    Where the stretches of an adapter `length` long that trim_read looks for start: every
    TRIM_PROBE_LEN bases over its first TRIM_PROBE_SPAN, and the last full stretch there
    """
    span = min(length, TRIM_PROBE_SPAN)
    if span <= TRIM_PROBE_LEN:
        return [0]
    out = list(range(0, span - TRIM_PROBE_LEN + 1, TRIM_PROBE_LEN))
    if out[-1] != span - TRIM_PROBE_LEN:
        out.append(span - TRIM_PROBE_LEN)
    return out


def adapter_cuts(seq, adapters):
    """
    (left, right, [names]): where adapter sequence ends on the left of `seq` and starts on
    its right. Read-through is the adapter's start following an insert, so on the right it's
    the adapter as given, and on the left its reverse complement (the same molecule seen from
    the other read). An adapter's start cut off by the end of `seq` is found too, and so is
    one missing its first bases where `seq` starts partway into it (by a later stretch of it,
    as in trim_read).
    """
    left, right, names = 0, len(seq), []
    for name, a in adapters.items():
        probe = a[:20]
        # 20-base stretches starting where trim_read's probes do (at least 12 bases long)
        windows = [(0, probe)] + [(o, a[o:o + 20]) for o in probe_offsets(len(a))[1:] if len(a) - o >= 12]
        cut = None
        for o, w in windows:
            hits = _locate(w, seq, 2 if len(w) >= 15 else 1)
            if hits:
                cut = max(0, min(s for s, _ in hits) - o)
                break
        else:
            for L in range(min(len(probe) - 1, len(seq)), TRIM_MIN_END_OVERLAP - 1, -1):
                if _near(seq[-L:], a[:L]):
                    cut = len(seq) - L
                    break
        if cut is not None and cut < right:
            right = cut
            names.append(name)
        cut = None
        for o, w in windows:
            hits = _locate(revcomp(w), seq, 2 if len(w) >= 15 else 1)
            if hits:
                cut = min(len(seq), max(end for _, end in hits) + 1 + o)
                break
        else:
            for L in range(min(len(probe) - 1, len(seq)), TRIM_MIN_END_OVERLAP - 1, -1):
                if _near(seq[:L], revcomp(a[:L])):
                    cut = L
                    break
        if cut is not None and cut > left:
            left = cut
            names.append(name)
    # no-signal G's (and their reverse complement, on the left)
    g = len(seq) - len(seq.rstrip("G"))
    if g >= POLY_G_RUN and len(seq) - g < right:
        right = len(seq) - g
        names.append("poly-G")
    c = len(seq) - len(seq.lstrip("C"))
    if c >= POLY_G_RUN and c > left:
        left = c
        names.append("poly-G")
    return left, max(left, right), sorted(set(names))


def mate_asymmetry(seq, by_mate, k, rc_all):
    """
    For each k-mer of `seq` (in order): (asymmetry, total count, share seen forward). An
    insert is read forward in one read of a pair and reverse complemented in the other, so
    its k-mers are about as common forward in R1 as reverse complemented in R2 (and the
    other way around): asymmetry near 0. Adapter read-through is in only one of the reads,
    and only forward: asymmetry near 1.
    """
    c1, c2 = by_mate[0], by_mate[1]
    enc, _, _ = encode([seq])
    codes, valid = window_codes(enc, k)
    out = []
    for x, v in zip(codes.tolist(), valid.tolist()):
        if not v:
            out.append((1.0, 0, 0.0))
            continue
        r = int(rc_all[x])
        a, b, c, d = int(c1[x]), int(c2[r]), int(c1[r]), int(c2[x])
        total = a + b + c + d
        paired = min(a, b) + min(c, d)
        out.append((1 - 2 * paired / total if total else 1.0, total, (a + d) / total if total else 0.0))
    return out


def split_by_mate_symmetry(seq, by_mate, k, rc_all, max_asym=MAX_INSERT_ASYMMETRY):
    """
    (start, end) of the longest stretch of `seq` whose k-mers are all symmetric between
    the reads of a pair (see mate_asymmetry), i.e., insert, or None if there's none.
    """
    asym = mate_asymmetry(seq, by_mate, k, rc_all)
    best, run_start = None, None
    for i, (a, _, _) in enumerate(asym + [(1.0, 0, 0.0)]):
        if a <= max_asym:
            if run_start is None:
                run_start = i
        elif run_start is not None:
            if best is None or i - run_start > best[1] - best[0]:
                best = (run_start, i)
            run_start = None
    if best is None:
        return None
    return best[0], best[1] - 1 + k


def _is_known(seq, known):
    probe = seq[:ADAPTER_PROBE_LEN]
    return any(_locate(probe, a, 1) or _locate(probe, revcomp(a), 1) for a in known.values())


def _read_through_like(seq, by_mate, k, rc_all):
    """whether all of `seq` is in only one read of a pair, and forward, as read-through is"""
    asym = mate_asymmetry(seq, by_mate, k, rc_all)
    total = sum(t for _, t, _ in asym)
    fwd = sum(t * share for _, t, share in asym)
    return bool(total) and fwd >= 0.9 * total and min(a for a, _, _ in asym) >= 1 - MAX_INSERT_ASYMMETRY


def infer_adapters(found, by_mate, k, rc_all, known):
    """
    Adapter read-through not already in `known`, from paired reads: in a found sequence with
    an insert (k-mers symmetric between the reads of a pair, see mate_asymmetry), stretches
    on either side seen in only one read and only forward. Where adapters share a start (as
    Illumina's R1 and R2 adapters do), extension stops there, so found sequences that are
    all read-through and begin with that start are taken as the fuller adapters. Returns
    {name: sequence oriented as read-through}.
    """
    starts = []
    for f in found:
        cons = f["consensus"]
        span = split_by_mate_symmetry(cons, by_mate, k, rc_all)
        if span is None:
            continue
        for flank in (cons[span[1]:], revcomp(cons[:span[0]])):
            if len(flank) >= TRIM_PROBE_LEN and _read_through_like(flank, by_mate, k, rc_all) \
                    and not _is_known(flank, known) and not _is_known(flank, dict(enumerate(starts))) \
                    and not any(_locate(a[:ADAPTER_PROBE_LEN], flank, 1) for a in starts):
                starts.append(flank)

    adapters = []
    for start in starts:
        fuller = []
        for f in found:
            for oriented in (f["consensus"], revcomp(f["consensus"])):
                if len(oriented) > len(start) and oriented.startswith(start[:TRIM_PROBE_LEN]) \
                        and _read_through_like(oriented, by_mate, k, rc_all) and oriented not in fuller:
                    fuller.append(oriented)
        # the longest of each that isn't just a piece of another
        fuller = [a for a in sorted(fuller, key=len, reverse=True)
                  if not any(b != a and b.startswith(a) for b in fuller)]
        adapters.extend(fuller or [start])
    return {f"inferred_adapter_{i}": a for i, a in enumerate(adapters, start=1)}


### collapsing ###

def low_complexity(seq):
    """
    flagged the same way cov-analyzer flags regions: at most LOW_COMPLEXITY of the possible
    distinct 3-mers present (see general.linguistic_complexity)
    """
    return linguistic_complexity(seq) <= LOW_COMPLEXITY


def allowed_edits(length, max_diff):
    return int(math.floor(max_diff * length))


def within(short, long, max_diff):
    """whether `short` aligns within `long` (either strand) with at most max_diff edits per base"""
    e = allowed_edits(len(short), max_diff)
    return bool(_locate(short, long, e) or _locate(revcomp(short), long, e))


def is_variant(s, t, max_diff):
    """
    Whether `s` is a variant of (more abundant) `t`: within it, at any length (a piece of it,
    or one of its error variants), or nearly the same length with most of `s` within it
    (e.g., missing a base or few at an end, carrying a bit of untrimmed adapter, or shifted).
    A sequence containing `t` with more than a little extra (e.g., a primer dimer) isn't.
    """
    if within(s, t, max_diff):
        return True
    slack = max(3, int(VARIANT_SLACK * min(len(s), len(t))))
    if abs(len(s) - len(t)) > slack:
        return False
    core = s[slack:len(s) - slack]
    return len(core) >= MIN_INSERT_LEN and within(core, t, max_diff)


def collapse(inserts, max_diff, adapters=None):
    """
    Drops inserts that are a variant of a more abundant one (see is_variant), or pieces of
    adapters left by errors that kept a read from being trimmed. Annotates the rest with the
    others they contain (e.g., a primer dimer, both primers).
    """
    kept = []
    for f in sorted(inserts, key=lambda x: -x["seed_count"]):
        s = f["insert"]
        if adapters and any(within(s, a, max_diff) for a in adapters.values()):
            continue
        if not any(is_variant(s, g["insert"], max_diff) for g in kept):
            kept.append(f)
    for f in kept:
        f["contains"] = [g for g in kept if g is not f and len(g["insert"]) < len(f["insert"])
                         and within(g["insert"], f["insert"], max_diff)]
    return kept


### pass 3: counting found sequences in reads ###

def _adapter_list(adapters):
    return None if adapters is None else list(adapters.items())


def trim_read(seq, probes):
    """
    The read cut where adapter read-through starts, and what cut it (an adapter's name,
    names joined by "|" for adapters sharing a start where the read ends too soon to tell
    them apart, "poly-G", or None). `probes` is [(name, adapter sequence)], each as it reads
    following an insert. An adapter is found by its first TRIM_PROBE_LEN bases, or, in case
    of an error there or a read starting partway into the adapter, a later stretch of that
    length (see probe_offsets), cutting where the adapter would start, or from the read's
    start if that's before it; or by at least TRIM_MIN_END_OVERLAP of its start where a read
    ends. A run of POLY_G_RUN or more G's at the end is cut too.
    """
    return native.trim_read(seq, list(probes))


def Matcher(seqs, max_diff, min_read_len=MIN_PARTIAL_LEN, groups=None):
    """
    Places found sequences in (trimmed) reads: `.hits(read)` gives {seq index: (start in
    read, or None if the read is within it; is reverse complement; edits)}, the best
    placement per found sequence. A read holds one if it aligns within the read with at
    most `max_diff` edits per base, or, for one longer than the read, if the read (at least
    min_read_len long) aligns within it. A read holding any in full isn't counted as within
    others, and one only within found sequences counts for the shortest (e.g., a primer
    missing a few bases, not the dimer that also holds it). A q-gram prefilter rules most
    out first, with q chosen so a sequence within its allowed edits must share an exact
    q-mer with the read (pigeonhole), so it can't cause a miss (`.hits(read, False)` skips
    it, for checking). `.count_in_reads(...)` counts over read files (see count_in_reads).
    `groups`, if given, is which sequence each of `seqs` counts as (e.g., the variants of a
    degenerate primer as one), and results are indexed by those.
    """
    return native.Matcher(list(seqs), max_diff, min_read_len, None if groups is None else list(groups))


def count_in_reads(read_paths, seqs, adapters, max_diff, min_insert_len=MIN_INSERT_LEN, run=quietly,
                   desc="", groups=None):
    """
    Counts the found sequences in reads, trimmed with `adapters`. Returns per-file tallies
    (keyed by "R1"/"R2", or "reads" for one file) and, for paired reads, pair-level ones,
    indexed by sequence (or by group, if `groups` is given; see Matcher).
    """
    tallies, pairs = run(desc, Matcher(seqs, max_diff, groups=groups).count_in_reads, [str(p) for p in read_paths],
                         list(adapters.items()), min_insert_len)
    labels = ["R1", "R2"] if len(read_paths) == 2 else ["reads"]
    out = {}
    for lab, t in zip(labels, tallies):
        out[lab] = dict(
            reads=t["reads"],
            reads_with=Counter(dict(enumerate(t["reads_with"]))),
            full=Counter(dict(enumerate(t["full"]))),
            rc=Counter(dict(enumerate(t["rc"]))),
            edits=defaultdict(Counter, {i: Counter(dict(e)) for i, e in enumerate(t["edits"])}),
            starts=defaultdict(Counter, {i: Counter(dict(st)) for i, st in enumerate(t["starts"])}),
            adapters=Counter(dict(t["adapters"])),
            adapter_only=t["adapter_only"],
            explained=t["explained"],
        )
    if pairs is not None:
        pairs = dict(pairs, pairs_with=Counter(dict(enumerate(pairs["pairs_with"]))),
                     pairs_both=Counter(dict(enumerate(pairs["pairs_both"]))))
    return out, pairs


### putting it together ###

def read_fasta(path, unnamed="seq"):
    """[(name, sequence)] from a fasta, names being the first word of headers (`unnamed`_N if blank)"""
    records, name, chunks = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0] if len(line) > 1 else f"{unnamed}_{len(records) + 1}", []
            elif line:
                chunks.append(line)
    if name:
        records.append((name, "".join(chunks).upper()))
    return records


def load_adapters(path=None, use_builtin=True):
    adapters = dict(BUILTIN_ADAPTERS) if use_builtin else {}
    if path:
        for name, seq in read_fasta(path, "adapter"):
            adapters[name] = seq
    return adapters


def expand_degenerate(seq, max_variants=MAX_VARIANTS):
    """every sequence a degenerate one (IUPAC codes) stands for, in a set order"""
    variants = [""]
    for b in seq:
        variants = [v + x for v in variants for x in IUPAC[b]]
        if len(variants) > max_variants:
            raise ValueError(f"A sequence to count stands for more than {max_variants:,} sequences with its "
                             f"degenerate bases (IUPAC codes), which is more than can be counted: {seq}")
    return variants


def load_seqs_to_count(path):
    """
    {name: sequence} from a fasta of sequences to count (U read as T, IUPAC codes allowed);
    raises ValueError for a sequence that's empty or has anything else in it, or a repeated name
    """
    seqs = {}
    for name, seq in read_fasta(path):
        seq = seq.replace("U", "T")
        if name in seqs:
            raise ValueError(f"The sequences to count have more than one named '{name}', but names need to be unique.")
        if not seq:
            raise ValueError(f"The sequence to count named '{name}' is empty.")
        bad = sorted(set(seq) - set(IUPAC))
        if bad:
            raise ValueError(f"The sequence to count named '{name}' has characters that aren't bases or "
                             f"IUPAC codes: {''.join(bad)}")
        seqs[name] = seq
    if not seqs:
        raise ValueError(f"No sequences were found in '{path}'.")
    return seqs


def find_overrep(read_paths, k=DEFAULT_K, min_perc=DEFAULT_MIN_PERC, min_fold=DEFAULT_MIN_FOLD,
                 max_diff=DEFAULT_MAX_DIFF, max_len=DEFAULT_MAX_LEN, adapters=None, jobs=None,
                 max_reads_for_discovery=DEFAULT_MAX_READS_FOR_DISCOVERY, seqs_to_count=None, verbose=False):
    """
    Finds and counts overrepresented sequences in reads (see module notes). Returns a dict
    with the found sequences, the adapters used (built-in, given and inferred), and the
    tallies. `min_perc` is the percent of reads a sequence must be in (as a fraction, `min_frac`,
    everywhere below). `jobs` sets the threads used (once per process; all cores if never set). Only the
    first `max_reads_for_discovery` reads of each file (0: all) are used to find the sequences;
    counting them (the last pass) always uses all reads. `seqs_to_count` ({name: sequence},
    IUPAC codes allowed) are counted too, and always reported, before what's found; a found
    sequence that's just a variant of one of them is left out for it. With
    `verbose`, each pass over the reads shows a progress bar, and what was found along the way
    is printed. Problems
    with the read files (not fasta/fastq, a record cut short, pairs out of step) raise
    ValueError, and ones reading them OSError.
    """
    min_frac = min_perc / 100
    if jobs:
        native.set_threads(jobs)
    adapters = dict(BUILTIN_ADAPTERS) if adapters is None else dict(adapters)
    run = with_progress_bar(read_paths) if verbose else quietly
    say = print if verbose else (lambda *a, **kw: None)

    # first, any adapters not already known, from the most abundant k-mers in the raw reads
    # (otherwise every insert's end next to an adapter makes for an overrepresented k-mer)
    n_disc = max_reads_for_discovery
    raw = KmerCounts(read_paths, k, run=run, desc="Counting k-mers and identifying adapters (pass 1 of 4)...",
                     max_reads=n_disc)
    candidates, raw_both = find_candidates(raw, min_frac, min_fold)
    _, rc_all, _ = raw.canonical()
    inferred = {}
    if len(read_paths) == 2:
        seeds = extend_all(candidates[:ADAPTER_SEARCH_SEEDS], raw_both, rc_all, k, max_len)
        inferred = infer_adapters(seeds, raw.by_mate, k, rc_all, adapters)
    all_adapters = dict(adapters, **inferred)
    if inferred:
        say(f"      Inferred {len(inferred)} adapter(s) not among those given (see identified-adapters.tsv)\n")

    # then counting again with adapter read-through trimmed off
    counts = KmerCounts(read_paths, k, all_adapters, run=run,
                        desc="Counting k-mers excluding adapters (pass 2 of 4)...", max_reads=n_disc)
    candidates, both = find_candidates(counts, min_frac, min_fold, error_ref=raw_both)
    del raw, counts
    say(f"      {len(candidates):,} candidate seed k-mers identified\n")

    extended = refine_ends(extend_all(candidates, both, rc_all, k, max_len), read_paths, all_adapters, k,
                           max_len, run=run, desc="Extending candidates using the reads that hold them (pass 3 of 4)...",
                           max_reads=n_disc)

    inserts = []
    for f in extended:
        left, right, _ = adapter_cuts(f["consensus"], all_adapters)
        insert = f["consensus"][left:right]
        if len(insert) < max(MIN_INSERT_LEN, k):
            continue
        f.update(insert=insert)
        inserts.append(f)
    kept = collapse(inserts, max_diff, all_adapters)
    # one orientation for each while sorting things out (whichever strand sorts first); each is
    # turned to how it's most often seen in reads at the end (orient_as_most_seen)
    for f in kept:
        f["insert"] = min(f["insert"], revcomp(f["insert"]))
        f["low_complexity"] = low_complexity(f["insert"])
        f["given"] = False

    # sequences given to count: a found one that's just a variant of one (e.g., the same
    # primer) is left out for it
    given = []
    for name, seq in (seqs_to_count or {}).items():
        variants = expand_degenerate(seq)
        given.append(dict(id=name, insert=seq, variants=variants, given=True,
                          low_complexity=low_complexity(variants[0])))
    if given:
        kept = [f for f in kept if not any(is_variant(f["insert"], v, max_diff)
                                          for g in given for v in g["variants"][:COMPARE_VARIANTS])]
    n_found = len(kept)
    kept += given
    for f in kept:
        f["contains"] = [g for g in kept if g is not f and len(g["insert"]) < len(f["insert"])
                         and _within_any(g, f, max_diff)]

    seqs, groups = [], []
    for i, f in enumerate(kept):
        for v in f.get("variants", [f["insert"]]):
            seqs.append(v)
            groups.append(i)
    what = f"the {n_found:,} determined overrepresented sequence(s)"
    if given:
        what += f" and {len(given):,} given sequence(s)"
    tallies, pairs = count_in_reads(read_paths, seqs, all_adapters, max_diff, groups=groups,
                                    run=run, desc=f"Counting {what} across all reads (pass 4 of 4)...")

    # ordered by how many reads hold them, leaving out any too rare once counted (a seed's
    # count can't be fully trusted to carry over to its extended sequence)
    def keep(i):
        return max(t["reads_with"][i] / t["reads"] if t["reads"] else 0 for t in tallies.values()) >= min_frac

    reads = lambda i: sum(t["reads_with"][i] for t in tallies.values())
    order = sorted((i for i in range(n_found) if keep(i)), key=lambda i: -reads(i))
    # a piece of another that's in nearly all the same reads is just that one, found from
    # another seed that stopped short; and a rare near-copy of a far more common one is that
    # one with errors in it
    def redundant(i, j):
        s, t = kept[i], kept[j]
        if (len(t["insert"]) > len(s["insert"]) and reads(j) >= CONTAINED_SHARE * reads(i)
                and _within_any(s, t, max_diff)):
            return True
        return reads(i) <= ERROR_VARIANT_SHARE * reads(j) and _within_any(s, t, 2 * max_diff)

    # (given ones are always reported, first, in the order given)
    given_order = list(range(n_found, len(kept)))
    order = [i for i in order if not any(j != i and redundant(i, j) for j in order + given_order)]
    for new, i in enumerate(order, start=1):
        kept[i]["id"] = f"overrep_{new}"
    for i in set(range(n_found)) - set(order):
        kept[i]["id"] = None
    for i in order:
        orient_as_most_seen(i, kept[i], tallies)
    order = given_order + order
    return dict(found=kept, order=order, adapters=all_adapters, inferred=inferred,
                tallies=tallies, pairs=pairs, k=k, n_candidates=len(candidates))


def _within_any(s, t, max_diff):
    """whether found or given sequence `s` is within `t` (any of their variants, for given ones)"""
    return any(within(a, b, max_diff) for a in s.get("variants", [s["insert"]])[:COMPARE_VARIANTS]
               for b in t.get("variants", [t["insert"]])[:COMPARE_VARIANTS])


def orient_as_most_seen(i, f, tallies):
    """
    Turns found sequence `i` around (and its tallies' reverse-complement counts with it) if
    it's more often reverse complemented than not in the read file it's most common in, so
    it's reported as it's most often seen in reads.
    """
    top = max(tallies.values(), key=lambda t: t["reads_with"][i])
    if top["rc"][i] * 2 > top["reads_with"][i]:
        f["insert"] = revcomp(f["insert"])
        for t in tallies.values():
            t["rc"][i] = t["reads_with"][i] - t["rc"][i]


### outputs ###

def output_paths(output_dir, output_prefix=""):
    out = Path(output_dir)
    return {
        "seqs": str(out / f"{output_prefix}overrep-seqs.tsv"),
        "fasta": str(out / f"{output_prefix}overrep-seqs.fasta"),
        "starts": str(out / f"{output_prefix}start-positions.tsv"),
        "adapters": str(out / f"{output_prefix}identified-adapters.tsv"),
        "summary": str(out / f"{output_prefix}summary.txt"),
        "log": str(out / f"{output_prefix}command-execution-info.txt"),
    }


def setup_output_dir(output_dir, output_prefix, force_overwrite, full_cmd_executed):
    check_if_output_dir_exists(output_dir, force_overwrite)
    os.makedirs(output_dir)
    log_command_run(full_cmd_executed, output_dir, output_paths(output_dir, output_prefix)["log"])


def _pct(n, d, places=3):
    return f"{100 * n / d:.{places}f}" if d else "NA"


def summary_lines(result):
    """the summary, as printed (indented) and written to summary.txt; blank lines between blocks"""
    n_given = sum(result["found"][i]["given"] for i in result["order"])
    what = "a found or given sequence" if n_given else "a found sequence"
    rows = [("Overrepresented sequences found:", f"{len(result['order']) - n_given:,}")]
    if n_given:
        rows.append(("Given sequences counted:", f"{n_given:,}"))
    for lab, t in result["tallies"].items():
        name = "Reads" if lab == "reads" else f"{lab} reads"
        rows += [None,
                 (name + ":", f"{t['reads']:,}"),
                 (f"  holding {what}:", f"{t['explained']:,} ({_pct(t['explained'], t['reads'], 2)}%)")]
    p = result["pairs"]
    if p:
        rows += [None,
                 ("Read pairs:", f"{p['pairs']:,}"),
                 (f"  either read holding {what}:", f"{p['explained']:,} ({_pct(p['explained'], p['pairs'], 2)}%)"),
                 (f"  both reads holding {what}:",
                  f"{p['explained_both']:,} ({_pct(p['explained_both'], p['pairs'], 2)}%)")]
    width = max(44, max(len(r[0]) + 2 for r in rows if r))
    return ["" if r is None else f"{r[0]:<{width}}{r[1]}" for r in rows]


def _adapter_source(name, result):
    if name in result["inferred"]:
        return "inferred"
    if BUILTIN_ADAPTERS.get(name) == result["adapters"].get(name):
        return "builtin"
    return "given"


def write_outputs(result, output_dir, output_prefix=""):
    paths = output_paths(output_dir, output_prefix)
    found, tallies, pairs = result["found"], result["tallies"], result["pairs"]
    labels = list(tallies)

    header = ["id", "sequence", "length"]
    for lab in labels:
        pre = "" if lab == "reads" else f"{lab}_"
        header += [f"{pre}reads", f"{pre}pct_of_reads", f"{pre}pct_full_length", f"{pre}pct_revcomp",
                   f"{pre}mean_edits", f"{pre}top_start", f"{pre}pct_at_top_start"]
    if pairs:
        header += ["either_read", "either_read_pct", "both_reads", "both_reads_pct"]
    header += ["low_complexity", "contains", "source"]

    with open(paths["seqs"], "w") as out, open(paths["fasta"], "w") as fa:
        out.write("\t".join(header) + "\n")
        for i in result["order"]:
            f = found[i]
            row = [f["id"], f["insert"], str(len(f["insert"]))]
            for lab in labels:
                t = tallies[lab]
                n = t["reads_with"][i]
                edits = t["edits"][i]
                mean_ed = sum(e * c for e, c in edits.items()) / n if n else 0
                starts = t["starts"][i]
                top, top_n = starts.most_common(1)[0] if starts else ("NA", 0)
                row += [str(n), _pct(n, t["reads"]), _pct(t["full"][i], n), _pct(t["rc"][i], n),
                        f"{mean_ed:.2f}", str(top), _pct(top_n, t["full"][i])]
            if pairs:
                row += [str(pairs["pairs_with"][i]), _pct(pairs["pairs_with"][i], pairs["pairs"]),
                        str(pairs["pairs_both"][i]), _pct(pairs["pairs_both"][i], pairs["pairs"])]
            row += ["yes" if f["low_complexity"] else "no",
                    ",".join(g["id"] for g in f["contains"] if g["id"]) or "none",
                    "given" if f["given"] else "found"]
            out.write("\t".join(row) + "\n")
            fa.write(f">{f['id']}\n{f['insert']}\n")

    with open(paths["starts"], "w") as out:
        out.write("id\tread\tstart\tcount\n")
        for i in result["order"]:
            for lab in labels:
                for start, c in sorted(tallies[lab]["starts"][i].items()):
                    out.write(f"{found[i]['id']}\t{lab}\t{start}\t{c}\n")

    with open(paths["adapters"], "w") as out:
        out.write("adapter\tsource\tsequence" + "".join(
            "\treads\tpct_of_reads" if lab == "reads" else f"\t{lab}_reads\t{lab}_pct_of_reads" for lab in labels) + "\n")
        names = list(result["adapters"]) + ["poly-G"]
        for lab in labels:
            names += [n for n in tallies[lab]["adapters"] if n not in names]
        for name in names:
            if name == "poly-G":
                source = "poly-G"
            else:
                # (a "|" name is read-through that could be any of those adapters: they share a
                # start, and the read ends too soon to tell them apart)
                sources = []
                for n in name.split("|"):
                    src = _adapter_source(n, result)
                    if src not in sources:
                        sources.append(src)
                source = "|".join(sources)
            if name == "poly-G":
                seq = "G" * POLY_G_RUN + "..."
            elif "|" in name:
                # the start the adapters share (all a read ending in it can show)
                seqs = [result["adapters"][n] for n in name.split("|") if n in result["adapters"]]
                seq = os.path.commonprefix(seqs) + "..." if seqs else "NA"
            else:
                seq = result["adapters"].get(name, "NA")
            row = [name, source, seq]
            for lab in labels:
                n = tallies[lab]["adapters"][name]
                row += [str(n), _pct(n, tallies[lab]["reads"])]
            out.write("\t".join(row) + "\n")

    with open(paths["summary"], "w") as out:
        out.write("\n".join(summary_lines(result)) + "\n")

    return paths
