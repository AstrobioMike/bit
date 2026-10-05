import gzip
import random
import edlib # type: ignore
import pytest # type: ignore
from bit.modules.assign_reads import (revcomp,
                                      build_kmer_index,
                                      build_position_index,
                                      get_candidates,
                                      ExactMatcher,
                                      EditMatcher,
                                      summarize_distances,
                                      find_identical_refs,
                                      load_refs,
                                      iter_pairs,
                                      FastxReader,
                                      assign_reads)
from bit.tests.utils import run_cli


K = 15


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def snp(seq, pos):
    alt = "A" if seq[pos] != "A" else "C"
    return seq[:pos] + alt + seq[pos + 1:]


@pytest.fixture
def seqs():
    rng = random.Random(42)
    base = rand_seq(600, rng)
    return {
        "p1": base,
        "p2": snp(base, 300),   # one SNP vs p1
        "p3": rand_seq(800, rng),
    }


def refs_list(seqs):
    return [(name, "refs.fa", s) for name, s in seqs.items()]


def hits_for(read, refs, circular, k=K, n_samples=10, min_frac=0.0):
    matcher = ExactMatcher(refs, k, n_samples, circular)
    return [refs[i][0] for i in matcher.hits(read, min_frac)]


def brute_force_hits(read, refs, circular):
    hits = []
    for i, (_, _, s) in enumerate(refs):
        if circular:
            search = s * (len(read) // len(s) + 2)
        else:
            search = s
        if read in search or revcomp(read) in search:
            hits.append(i)
    return hits


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, s in records:
            f.write(f">{name}\n{s}\n")


def write_fastq(path, records, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as f:
        for name, s in records:
            f.write(f"@{name}\n{s}\n+\n{'I' * len(s)}\n")


### matching logic ###

def test_exact_match_both_strands(seqs):
    refs = refs_list(seqs)
    read = seqs["p3"][100:400]
    assert hits_for(read, refs, circular=False) == ["p3"]
    assert hits_for(revcomp(read), refs, circular=False) == ["p3"]


def test_one_mismatch_is_not_a_match(seqs):
    refs = refs_list(seqs)
    read = snp(seqs["p3"][100:400], 150)
    assert hits_for(read, refs, circular=False) == []


def test_origin_spanning_read_needs_circular(seqs):
    refs = refs_list(seqs)
    p3 = seqs["p3"]
    read = p3[-100:] + p3[:100]
    assert hits_for(read, refs, circular=False) == []
    assert hits_for(read, refs, circular=True) == ["p3"]
    assert hits_for(revcomp(read), refs, circular=True) == ["p3"]


def test_full_length_rotated_read_circular(seqs):
    refs = refs_list(seqs)
    p1 = seqs["p1"]
    read = p1[250:] + p1[:250]  # covers the SNP at 300
    assert hits_for(read, refs, circular=True) == ["p1"]


def test_read_longer_than_circular_ref_matches(seqs):
    refs = refs_list(seqs)
    p3 = seqs["p3"]
    read = (p3 * 3)[50:50 + 2000]  # e.g., from a multimer
    assert hits_for(read, refs, circular=True) == ["p3"]
    assert hits_for(read, refs, circular=False) == []


def test_ambiguous_when_read_misses_snp(seqs):
    refs = refs_list(seqs)
    read = seqs["p1"][0:250]
    assert hits_for(read, refs, circular=False) == ["p1", "p2"]


def test_unique_when_read_covers_snp(seqs):
    refs = refs_list(seqs)
    assert hits_for(seqs["p1"][250:350], refs, circular=False) == ["p1"]
    assert hits_for(seqs["p2"][250:350], refs, circular=False) == ["p2"]


def test_min_frac_of_ref(seqs):
    refs = refs_list(seqs)
    read = seqs["p3"][:400]
    assert hits_for(read, refs, circular=False, min_frac=0.5) == ["p3"]
    assert hits_for(read, refs, circular=False, min_frac=0.9) == []


def test_read_shorter_than_k_still_checked(seqs):
    refs = refs_list(seqs)
    read = seqs["p3"][10:20]
    assert len(read) < K
    assert "p3" in hits_for(read, refs, circular=False)


def test_prefilter_never_excludes_true_matches(seqs):
    """Every substring of a ref (any strand, any rotation) keeps that ref as a candidate."""
    refs = refs_list(seqs)
    index = build_kmer_index(refs, K, circular=True)
    rng = random.Random(1)
    for i, (_, _, s) in enumerate(refs):
        circ = s * 2
        for _ in range(200):
            start = rng.randrange(len(s))
            length = rng.randrange(K, len(s) + 1)
            read = circ[start:start + length]
            if rng.random() < 0.5:
                read = revcomp(read)
            assert i in get_candidates(read, index, K, 25)


@pytest.mark.parametrize("circular", [False, True])
def test_matcher_agrees_with_brute_force(seqs, circular):
    """
    Random reads (exact, with an error, any strand/rotation, short and long) get the
    same hits from the prefilter + anchored comparisons as from brute-force searching.
    """
    refs = refs_list(seqs)
    matcher = ExactMatcher(refs, K, 25, circular)
    rng = random.Random(3)
    n_with_hits = 0
    for _ in range(600):
        _, _, s = rng.choice(refs)
        source = s * 3 if circular else s
        max_len = len(source) if circular else len(s)
        length = rng.choice([rng.randrange(5, 40), rng.randrange(40, max_len + 1)])
        start = rng.randrange(0, (len(s) if circular else len(s) - length + 1))
        read = source[start:start + length]
        if rng.random() < 0.3:
            read = snp(read, rng.randrange(len(read)))
        if rng.random() < 0.5:
            read = revcomp(read)
        expected = brute_force_hits(read, refs, circular)
        n_with_hits += bool(expected)
        assert matcher.hits(read) == expected, (circular, start, length)
    assert n_with_hits > 300


def test_wraparound_matching_at_every_rotation(seqs):
    """
    Full-length, slightly-longer-than-ref, multi-pass, and short origin-spanning reads,
    from every start position, match the circular ref and agree with brute force.
    """
    p3 = seqs["p3"]
    L = len(p3)
    refs = [("p3", "f", p3)]
    matcher = ExactMatcher(refs, K, 25, circular=True)
    source = p3 * 4
    for start in range(L):
        for length in (L, L + 1, 2 * L + 5, 20):
            read = source[start:start + length]
            assert matcher._matches_at(p3, read, start)
            assert matcher._contains(p3, read)
            assert matcher.hits(read) == [0] == brute_force_hits(read, refs, True)


def test_wraparound_rejects_mismatch_in_each_segment(seqs):
    p3 = seqs["p3"]
    L = len(p3)
    matcher = ExactMatcher([("p3", "f", p3)], K, 25, circular=True)
    start = L - 100
    read = (p3 * 3)[start:start + 2 * L + 50]
    # a mismatch before the origin, in the full pass, and in the final partial pass
    for pos in (50, 100 + L // 2, 2 * L + 25):
        assert not matcher._matches_at(p3, snp(read, pos), start)


def test_contains_short_read_at_junction_only_when_circular(seqs):
    p3 = seqs["p3"]
    read = p3[-10:] + p3[:10]
    linear = ExactMatcher([("p3", "f", p3)], K, 25, circular=False)
    circular = ExactMatcher([("p3", "f", p3)], K, 25, circular=True)
    assert not linear._contains(p3, read)
    assert circular._contains(p3, read)
    assert circular._contains(p3, p3[5])  # single base, no junction needed


def test_matcher_stores_refs_once(seqs):
    refs = refs_list(seqs)
    matcher = ExactMatcher(refs, K, 25, circular=True)
    assert all(m is s for m, (_, _, s) in zip(matcher.ref_seqs, refs))
    assert not hasattr(matcher, "searches")


### edit-tolerant matching ###

def add_edits(seq, n_edits, rng):
    """Applies n random substitutions/insertions/deletions."""
    for _ in range(n_edits):
        pos = rng.randrange(len(seq))
        kind = rng.choice(["sub", "ins", "del"])
        if kind == "sub":
            seq = snp(seq, pos)
        elif kind == "ins":
            seq = seq[:pos] + rng.choice("ACGT") + seq[pos:]
        else:
            seq = seq[:pos] + seq[pos + 1:]
    return seq


def brute_force_distances(read, refs, circular, max_edits, min_frac=0.0):
    out = {}
    for i, (_, _, s) in enumerate(refs):
        n, L = len(read), len(s)
        if n < min_frac * L or (not circular and n > L + max_edits):
            continue
        target = s * (n // L + 3) if circular else s
        dists = [edlib.align(r, target, mode="HW", task="distance", k=max_edits)["editDistance"]
                 for r in (read, revcomp(read))]
        dists = [d for d in dists if d != -1]
        if dists:
            out[i] = min(dists)
    return out


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("max_edits", [1, 3])
def test_edit_matcher_agrees_with_brute_force(seqs, circular, max_edits):
    """
    Reads with 0 to max_edits + 1 random edits, any strand/rotation, short (full-alignment
    fallback) and long (prefilter + anchored), get the same distances as brute-force alignment.
    """
    refs = refs_list(seqs)
    matcher = EditMatcher(refs, K, circular, max_edits)
    rng = random.Random(11 + max_edits)
    n_within = 0
    for _ in range(300):
        _, _, s = rng.choice(refs)
        source = s * 3 if circular else s
        length = rng.choice([rng.randrange(20, 60), rng.randrange(60, len(s) + 1)])
        start = rng.randrange(0, (len(s) if circular else len(s) - length + 1))
        read = add_edits(source[start:start + length], rng.randrange(0, max_edits + 2), rng)
        if rng.random() < 0.5:
            read = revcomp(read)
        expected = brute_force_distances(read, refs, circular, max_edits)
        n_within += bool(expected)
        assert matcher.distances(read) == expected, (circular, start, length)
    assert n_within > 150


def test_edit_matcher_respects_min_frac(seqs):
    refs = refs_list(seqs)
    matcher = EditMatcher(refs, K, False, 2)
    read = snp(seqs["p3"][:400], 100)
    assert matcher.distances(read, 0.4) == {2: 1}
    assert matcher.distances(read, 0.9) == {}


def test_closest_ref_wins_and_ties_are_ambiguous(seqs):
    refs = refs_list(seqs)
    matcher = EditMatcher(refs, K, False, 2)

    # covers p1/p2's distinguishing SNP, plus one error elsewhere: p1 at 1, p2 at 2
    read = snp(seqs["p1"][250:550], 200)
    assert summarize_distances(matcher.distances(read)) == ([0], 1, 2)

    # doesn't cover the SNP: tied
    read = snp(seqs["p1"][0:250], 100)
    assert summarize_distances(matcher.distances(read)) == ([0, 1], 1, None)

    assert summarize_distances({}) is None


def test_exact_matcher_distances_are_zero(seqs):
    refs = refs_list(seqs)
    matcher = ExactMatcher(refs, K, 25, False)
    assert matcher.distances(seqs["p1"][0:250]) == {0: 0, 1: 0}


def test_assign_reads_with_edits_end_to_end(tmp_path, seqs):
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    p1, p3 = seqs["p1"], seqs["p3"]
    reads = [
        ("exact_p3", p3[100:500]),
        ("one_edit_p1", snp(p1[250:550], 200)),      # p1 at 1 edit, p2 at 2
        ("tied", snp(p1[0:250], 100)),
        ("too_many", snp(snp(snp(p3[100:500], 50), 150), 250)),  # 3 edits > max of 2
    ]
    reads_fq = tmp_path / "reads.fq"
    write_fastq(reads_fq, reads)

    prefix = str(tmp_path / "ed")
    summary = assign_reads([str(refs_fa)], str(reads_fq), output_prefix=prefix, max_edits=2,
                           k=K, jobs=2, batch_size=1, show_progress=False)
    assert (summary["unique"], summary["ambiguous"]) == (2, 1)

    hits = read_hits_tsv(f"{prefix}-read-hits.tsv")
    assert set(hits) == {"exact_p3", "one_edit_p1", "tied"}
    assert hits["exact_p3"][2:] == ["unique", "1", "p3", "0", "NA"]
    assert hits["one_edit_p1"][2:] == ["unique", "1", "p1", "1", "2"]
    assert hits["tied"][2:] == ["ambiguous", "2", "p1,p2", "1", "NA"]


def test_assign_reads_pairs_ranked_by_combined_edits(tmp_path, seqs):
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    p1 = seqs["p1"]
    # mate 1 has an error and doesn't cover the SNP (tied at 1); mate 2 covers it exactly
    pairs = [("pair", snp(p1[0:150], 60), revcomp(p1[250:400]))]
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    write_fastq(r1, [(f"{n}/1", a) for n, a, _ in pairs])
    write_fastq(r2, [(f"{n}/2", b) for n, _, b in pairs])

    prefix = str(tmp_path / "pe-ed")
    assign_reads([str(refs_fa)], str(r1), read_2=str(r2), output_prefix=prefix,
                 max_edits=1, k=K, jobs=1, show_progress=False)
    hits = read_hits_tsv(f"{prefix}-read-hits.tsv")
    # p1: 1 + 0 = 1, p2: 1 + 1 = 2
    assert hits["pair"][2:] == ["unique", "1", "p1", "1", "2"]


@pytest.mark.parametrize("circular", [False, True])
def test_position_index_gaps_never_exceed_step(seqs, circular):
    refs = refs_list(seqs)
    step = 16
    pos_index = build_position_index(refs, K, circular, step)
    for i, (_, _, s) in enumerate(refs):
        positions = sorted(p for entries in pos_index.values() for j, p in entries if j == i)
        if circular:
            gaps = [b - a for a, b in zip(positions, positions[1:] + [positions[0] + len(s)])]
        else:
            assert positions[0] == 0 and positions[-1] == len(s) - K
            gaps = [b - a for a, b in zip(positions, positions[1:])]
        assert max(gaps) <= step


def test_index_shares_identical_sets(seqs):
    refs = refs_list(seqs)
    index = build_kmer_index(refs, K, circular=False)
    shared = [v for v in index.values() if v == frozenset({0, 1})]
    assert len(shared) > 1
    assert all(v is shared[0] for v in shared)


def test_find_identical_refs(seqs):
    p1 = seqs["p1"]
    refs = [("a", "f", p1), ("b", "f", revcomp(p1[200:] + p1[:200])), ("c", "f", seqs["p3"])]
    assert find_identical_refs(refs, circular=True) == [("a", "b")]
    assert find_identical_refs(refs, circular=False) == []


### input parsing ###

def test_fastx_reader_fasta_and_gz_fastq(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">s1 desc\nACGT\nACGT\n>s2\nTTTT\n")
    assert [(n, s) for n, s, _ in FastxReader(fa)] == [("s1", "ACGTACGT"), ("s2", "TTTT")]

    fq = tmp_path / "x.fq.gz"
    write_fastq(fq, [("r1", "ACGT"), ("r2", "GGCC")], gz=True)
    reader = FastxReader(fq)
    assert [(n, s, q) for n, s, q in reader] == [("r1", "ACGT", "IIII"), ("r2", "GGCC", "IIII")]
    assert reader.bytes_read() == reader.total_bytes


def test_load_refs_rejects_duplicate_names(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    write_fasta(a, [("x", "ACGT")])
    write_fasta(b, [("x", "TTTT")])
    with pytest.raises(SystemExit):
        load_refs([str(a), str(b)])


def test_iter_pairs_checks_names_and_counts(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    write_fastq(r1, [("a/1", "ACGT"), ("b/1", "ACGT")])
    write_fastq(r2, [("a/2", "TTTT"), ("b/2", "TTTT")])
    assert len(list(iter_pairs(FastxReader(r1), FastxReader(r2)))) == 2

    write_fastq(r2, [("a/2", "TTTT"), ("c/2", "TTTT")])
    with pytest.raises(SystemExit):
        list(iter_pairs(FastxReader(r1), FastxReader(r2)))

    write_fastq(r2, [("a/2", "TTTT")])
    with pytest.raises(SystemExit):
        list(iter_pairs(FastxReader(r1), FastxReader(r2)))


### end-to-end ###

def make_reads(seqs, rng):
    """single-end reads with known expected statuses"""
    reads, expected = [], {}
    p1, p2, p3 = seqs["p1"], seqs["p2"], seqs["p3"]
    cases = [
        ("full_p1", p1[250:] + p1[:250], "p1"),
        ("full_p2_rc", revcomp(p2[400:] + p2[:400]), "p2"),
        ("partial_shared", p1[0:200], "p1,p2"),
        ("origin_p3", p3[-150:] + p3[:150], "p3"),
        ("error_read", snp(p3[100:500], 50), None),
        ("junk", rand_seq(300, rng), None),
    ]
    for name, s, exp in cases:
        reads.append((name, s))
        expected[name] = exp
    return reads, expected


def read_hits_tsv(path):
    lines = open(path).read().splitlines()
    return {l.split("\t")[0]: l.split("\t") for l in lines[1:]}


@pytest.mark.parametrize("jobs", [1, 3])
def test_assign_reads_single_end(tmp_path, seqs, jobs):
    rng = random.Random(7)
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    reads, expected = make_reads(seqs, rng)
    reads_fq = tmp_path / "reads.fq.gz"
    write_fastq(reads_fq, reads, gz=True)

    prefix = str(tmp_path / "out")
    summary = assign_reads([str(refs_fa)], str(reads_fq), output_prefix=prefix, circular=True,
                           k=K, write_reads=True, jobs=jobs, batch_size=2, show_progress=False)

    assert summary["total"] == len(reads)
    assert summary["unique"] == 3
    assert summary["ambiguous"] == 1

    hits = read_hits_tsv(f"{prefix}-read-hits.tsv")
    for name, exp in expected.items():
        if exp is None:
            assert name not in hits
        else:
            assert hits[name][4] == exp
            assert hits[name][2] == ("unique" if "," not in exp else "ambiguous")

    # output order follows input order regardless of jobs
    assert list(hits) == [n for n, _ in reads if expected[n]]

    summary_lines = open(f"{prefix}-summary.tsv").read().splitlines()
    assert summary_lines[1].split("\t") == ["p1", "refs.fa", "600", "1", "1"]

    written = (tmp_path / "out-reads" / "p1.fastq").read_text().split("\n")
    assert written[0] == "@full_p1"


def test_assign_reads_paired_end(tmp_path, seqs):
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    p1, p2, p3 = seqs["p1"], seqs["p2"], seqs["p3"]

    # mate 1 forward, mate 2 reverse-complemented from downstream, as in real libraries
    pairs = [
        ("pair_unique", p1[280:380], revcomp(p1[450:550])),       # mate 1 covers the SNP
        ("pair_ambig", p1[0:100], revcomp(p1[450:550])),          # neither covers it
        ("pair_one_bad", p3[0:100], snp(revcomp(p3[300:400]), 20)),
        ("pair_split", p3[0:100], revcomp(p1[450:550])),          # mates match different refs
    ]
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    write_fastq(r1, [(f"{n}/1", a) for n, a, _ in pairs])
    write_fastq(r2, [(f"{n}/2", b) for n, _, b in pairs])

    prefix = str(tmp_path / "pe")
    summary = assign_reads([str(refs_fa)], str(r1), read_2=str(r2), output_prefix=prefix,
                           k=K, write_reads=True, jobs=1, show_progress=False)

    assert summary["unit"] == "pairs"
    hits = read_hits_tsv(f"{prefix}-read-hits.tsv")
    assert set(hits) == {"pair_unique", "pair_ambig"}
    assert hits["pair_unique"][1] == "100,100"
    assert hits["pair_unique"][4] == "p1"
    assert hits["pair_ambig"][4] == "p1,p2"

    assert (tmp_path / "pe-reads" / "p1_R1.fastq").read_text().startswith("@pair_unique/1")
    assert (tmp_path / "pe-reads" / "p1_R2.fastq").read_text().startswith("@pair_unique/2")


def test_cli_runs_and_respects_force(tmp_path, seqs):
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    reads_fq = tmp_path / "reads.fq"
    write_fastq(reads_fq, [("r", seqs["p3"][:300])])
    prefix = str(tmp_path / "cli")

    cmd = ["bit", "assign-reads", "-r", str(refs_fa), "-1", str(reads_fq), "-o", prefix,
           "-k", str(K), "-j", "1"]
    result = run_cli(cmd)
    assert "Uniquely assigned" in result.stdout

    import subprocess
    rerun = subprocess.run(cmd, capture_output=True, text=True)
    assert rerun.returncode != 0

    run_cli(cmd + ["-F"])


def test_cli_rejects_min_frac_with_pairs(tmp_path, seqs):
    import subprocess
    refs_fa = tmp_path / "refs.fa"
    write_fasta(refs_fa, seqs.items())
    reads_fq = tmp_path / "reads.fq"
    write_fastq(reads_fq, [("r", seqs["p3"][:300])])
    result = subprocess.run(["bit", "assign-reads", "-r", str(refs_fa), "-1", str(reads_fq),
                             "-2", str(reads_fq), "--min-frac-of-ref", "0.5"],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "single-end" in result.stderr
