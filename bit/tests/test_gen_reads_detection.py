import random
import pytest # type: ignore

from bit.modules.gen_reads_detection import DetectionTracker


def _brute(contig_lengths, spans_by_contig):
    """ reference detection via a per-base bitmap: covered_bases / total. """
    total = sum(contig_lengths)
    if total == 0:
        return 1.0, True
    covered = 0
    for ci, length in enumerate(contig_lengths):
        bits = bytearray(length)
        for s, e in spans_by_contig[ci]:
            for x in range(s, e):
                bits[x] = 1
        covered += sum(bits)
    return covered / total, (covered == total)


# ─── basic semantics ──────────────────────────────────────────────────────────

def test_empty_genome_is_fully_detected():
    tr = DetectionTracker([])
    assert tr.detection() == 1.0
    assert tr.done is True


def test_all_zero_length_contigs():
    tr = DetectionTracker([0, 0])
    assert tr.detection() == 1.0
    assert tr.done is True


def test_single_contig_full_coverage_sets_done():
    tr = DetectionTracker([10])
    tr.add(0, 0, 10)
    assert tr.detection() == 1.0
    assert tr.done is True


def test_partial_coverage_fraction():
    tr = DetectionTracker([100])
    tr.add(0, 0, 25)
    assert tr.detection() == pytest.approx(0.25)
    assert tr.done is False
    assert tr.covered_bases() == 25
    assert tr.genome_size() == 100


def test_overlapping_spans_not_double_counted():
    tr = DetectionTracker([100])
    tr.add(0, 10, 30)
    first = tr.detection()
    tr.add(0, 10, 30)               # exact repeat
    assert tr.detection() == first
    tr.add(0, 20, 40)               # partial overlap, only 30..40 is new
    assert tr.covered_bases() == 30  # 10..40


def test_adjacent_spans_merge_to_full():
    tr = DetectionTracker([100])
    tr.add(0, 0, 50)
    tr.add(0, 50, 100)
    assert tr.done is True
    assert tr.detection() == 1.0


def test_interior_span_splits_gap():
    # covering the middle of a contig leaves two gaps; detection reflects only
    # the covered middle, and the contig is not yet done
    tr = DetectionTracker([100])
    tr.add(0, 40, 60)
    assert tr.covered_bases() == 20
    assert tr.done is False
    # fill the two flanks
    tr.add(0, 0, 40)
    tr.add(0, 60, 100)
    assert tr.done is True


def test_done_short_circuits_further_adds():
    tr = DetectionTracker([10])
    tr.add(0, 0, 10)
    assert tr.done is True
    # adding more after done must neither error nor change anything
    tr.add(0, 0, 5)
    tr.add(0, 3, 8)
    assert tr.detection() == 1.0


def test_multi_contig_independent_tracking():
    tr = DetectionTracker([50, 50])
    tr.add(0, 0, 50)               # first contig fully covered
    assert tr.detection() == pytest.approx(0.5)
    assert tr.done is False
    tr.add(1, 0, 25)               # half the second
    assert tr.detection() == pytest.approx(0.75)
    tr.add(1, 25, 50)
    assert tr.done is True


def test_zero_length_contig_among_real_ones():
    # a zero-length contig contributes nothing and must not block done
    tr = DetectionTracker([10, 0, 10])
    tr.add(0, 0, 10)
    tr.add(2, 0, 10)
    assert tr.done is True
    assert tr.detection() == 1.0


# ─── randomized agreement with a brute-force bitmap ───────────────────────────

def test_randomized_matches_brute_force():
    rng = random.Random(1234)
    for _ in range(400):
        ncontigs = rng.randint(1, 5)
        lengths = [rng.randint(1, 300) for _ in range(ncontigs)]
        if rng.random() < 0.1:
            lengths[rng.randrange(ncontigs)] = 0

        spans_by_contig = [[] for _ in range(ncontigs)]
        flat = []
        for _ in range(rng.randint(0, 120)):
            ci = rng.randrange(ncontigs)
            length = lengths[ci]
            if length == 0:
                continue
            s = rng.randint(0, length - 1)
            e = min(length, s + rng.randint(1, 40))
            spans_by_contig[ci].append((s, e))
            flat.append((ci, s, e))

        ref_det, ref_done = _brute(lengths, spans_by_contig)

        tr = DetectionTracker(lengths)
        rng.shuffle(flat)            # order independence
        for ci, s, e in flat:
            tr.add(ci, s, e)

        assert tr.detection() == pytest.approx(ref_det)
        assert tr.done == ref_done
        assert tr.covered_bases() == round(ref_det * sum(lengths))


# ─── covered_intervals (the complement of the gap list; basis of the GTA) ──────

def _brute_intervals(length, spans):
    """ reference covered intervals via a per-base bitmap. """
    if length == 0:
        return []
    bits = bytearray(length)
    for s, e in spans:
        for x in range(s, e):
            bits[x] = 1
    out = []
    i = 0
    while i < length:
        if bits[i]:
            j = i
            while j < length and bits[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def test_covered_intervals_empty_when_no_reads():
    tr = DetectionTracker([50])
    assert tr.covered_intervals(0) == []


def test_covered_intervals_full_contig():
    tr = DetectionTracker([20])
    tr.add(0, 0, 20)
    assert tr.covered_intervals(0) == [(0, 20)]


def test_covered_intervals_zero_length_contig():
    tr = DetectionTracker([0, 5])
    tr.add(1, 0, 5)
    assert tr.covered_intervals(0) == []
    assert tr.covered_intervals(1) == [(0, 5)]


def test_covered_intervals_merges_adjacent_and_overlapping():
    tr = DetectionTracker([100])
    tr.add(0, 10, 30)
    tr.add(0, 25, 40)   # overlaps previous -> merge to (10, 40)
    tr.add(0, 40, 50)   # abuts previous    -> merge to (10, 50)
    tr.add(0, 70, 80)   # separate island
    assert tr.covered_intervals(0) == [(10, 50), (70, 80)]


def test_covered_intervals_matches_covered_bases():
    tr = DetectionTracker([200])
    for s, e in [(5, 15), (14, 30), (60, 61), (100, 180)]:
        tr.add(0, s, e)
    total = sum(e - s for s, e in tr.covered_intervals(0))
    assert total == tr.covered_bases()


def test_covered_intervals_randomized_matches_brute_force():
    rng = random.Random(99)
    for _ in range(300):
        ncontigs = rng.randint(1, 4)
        lengths = [rng.randint(1, 200) for _ in range(ncontigs)]
        spans_by_contig = [[] for _ in range(ncontigs)]
        flat = []
        for _ in range(rng.randint(0, 80)):
            ci = rng.randrange(ncontigs)
            length = lengths[ci]
            s = rng.randint(0, length - 1)
            e = min(length, s + rng.randint(1, 30))
            spans_by_contig[ci].append((s, e))
            flat.append((ci, s, e))

        tr = DetectionTracker(lengths)
        rng.shuffle(flat)                 # order independence
        for ci, s, e in flat:
            tr.add(ci, s, e)

        for ci in range(ncontigs):
            assert tr.covered_intervals(ci) == _brute_intervals(
                lengths[ci], spans_by_contig[ci])
