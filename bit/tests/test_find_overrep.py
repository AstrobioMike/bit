import random
import subprocess
import numpy as np
import pytest # type: ignore
from bit.modules.find_overrep import (BUILTIN_ADAPTERS,
                                      Matcher,
                                      encode,
                                      window_codes,
                                      rc_codes,
                                      seq_to_code,
                                      code_to_seq,
                                      revcomp,
                                      not_error_variants,
                                      trim_read,
                                      extend_by_reads,
                                      is_variant,
                                      low_complexity,
                                      mate_asymmetry,
                                      split_by_mate_symmetry,
                                      load_adapters,
                                      find_overrep,
                                      write_outputs,
                                      summary_lines)
from bit.tests.utils import run_cli


A1 = BUILTIN_ADAPTERS["illumina_truseq_r1"]
A2 = BUILTIN_ADAPTERS["illumina_truseq_r2"]
P1 = "GTGCCAGCAGCCGCGGTAATAC"
P2 = "GGACTACACGGGTTTCTAAT"
DIMER = P1 + revcomp(P2)[4:]


def rand_seq(n, rng):
    return "".join(rng.choices("ACGT", k=n))


def mutate(seq, rng, rate=0.003):
    return "".join(rng.choice([x for x in "ACGT" if x != b]) if rng.random() < rate else b for b in seq)


def simulate(path_1, path_2, n, seed=7, read_len=50):
    """
    2x50 reads of a mostly-primer library: P1 (30% of molecules), P2 (15%), a P1-P2 dimer
    (2%), and the rest pieces of a 1-Mb random genome. Inserts go in either orientation, R1
    reads the insert then the R1 adapter, R2 its reverse complement then the R2 adapter.
    """
    rng = random.Random(seed)
    genome = rand_seq(1_000_000, rng)
    with open(path_1, "w") as o1, open(path_2, "w") as o2:
        for i in range(n):
            r = rng.random()
            if r < 0.30:
                insert = P1
            elif r < 0.45:
                insert = P2
            elif r < 0.47:
                insert = DIMER
            else:
                length = rng.randint(30, 200)
                start = rng.randint(0, len(genome) - length)
                insert = genome[start:start + length]
            if rng.random() < 0.5:
                insert = revcomp(insert)
            a = mutate((insert + A1)[:read_len], rng)
            b = mutate((revcomp(insert) + A2)[:read_len], rng)
            o1.write(f"@r{i}/1\n{a}\n+\n{'I' * len(a)}\n")
            o2.write(f"@r{i}/2\n{b}\n+\n{'I' * len(b)}\n")


@pytest.fixture(scope="module")
def reads(tmp_path_factory):
    d = tmp_path_factory.mktemp("overrep")
    r1, r2 = d / "r1.fq", d / "r2.fq"
    simulate(r1, r2, 20_000)
    return r1, r2


def found_counts(result):
    """{sequence: (reads holding it per file, pairs)}"""
    out = {}
    for i in result["order"]:
        f = result["found"][i]
        per_file = tuple(t["reads_with"][i] for t in result["tallies"].values())
        out[f["insert"]] = (per_file, result["pairs"]["pairs_with"][i] if result["pairs"] else None)
    return out


def canon(seq):
    return min(seq, revcomp(seq))


### k-mer basics ###

def test_window_codes_match_sequences():
    seqs = ["ACGTACGTAC", "TTNAAACCC", "GG"]
    enc, starts, lengths = encode(seqs)
    codes, valid = window_codes(enc, 3)
    got = [code_to_seq(int(c), 3) for c, v in zip(codes, valid) if v]
    want = [s[i:i + 3] for s in seqs for i in range(len(s) - 2) if "N" not in s[i:i + 3]]
    assert got == want
    assert starts.tolist() == [0, 11, 21]
    assert lengths.tolist() == [10, 9, 2]


def test_rc_codes():
    rng = random.Random(1)
    kmers = [rand_seq(11, rng) for _ in range(50)]
    codes = np.array([seq_to_code(s) for s in kmers], dtype=np.int64)
    assert [code_to_seq(int(c), 11) for c in rc_codes(codes, 11)] == [revcomp(s) for s in kmers]


def test_not_error_variants_drops_rare_one_offs():
    k = 5
    both = np.zeros(4 ** k, dtype=np.int64)
    main, one_off, other = seq_to_code("ACGTA"), seq_to_code("ACTTA"), seq_to_code("GGGCC")
    both[main], both[one_off], both[other] = 1000, 20, 20
    kept = not_error_variants(np.array([main, one_off, other], dtype=np.int64), both, k, ratio=20)
    assert kept.tolist() == [main, other]
    # not a variant if not rare enough next to it
    both[one_off] = 100
    kept = not_error_variants(np.array([main, one_off], dtype=np.int64), both, k, ratio=20)
    assert kept.tolist() == [main, one_off]


### trimming ###

def test_trim_read():
    probes = list(BUILTIN_ADAPTERS.items())
    # adapter read-through
    assert trim_read(P1 + A1[:28], probes) == (P1, "illumina_truseq_r1")
    assert trim_read(revcomp(P1) + A2[:28], probes) == (revcomp(P1), "illumina_truseq_r2")
    # an error in the adapter's first stretch, found by its next one
    assert trim_read(P1 + "AGTTCGGAAG" + A1[10:28], probes) == (P1, "illumina_truseq_r1")
    # cut off by the read end
    nextera = BUILTIN_ADAPTERS["illumina_nextera"]
    assert trim_read(P1 + nextera[:7], probes) == (P1, "illumina_nextera")
    assert trim_read(P1 + A1[:15], probes) == (P1, "illumina_truseq_r1")
    # too little of it to tell which (TruSeq R1 and R2 share their first 13 bases)
    assert trim_read(P1 + A1[:7], probes) == (P1, "illumina_truseq_r1|illumina_truseq_r2")
    assert trim_read(P1 + A1[:13], probes) == (P1, "illumina_truseq_r1|illumina_truseq_r2")
    # poly-G
    assert trim_read(P1 + "G" * 10, probes) == (P1, "poly-G")
    # nothing to trim
    assert trim_read(P1, probes) == (P1, None)


def test_load_adapters(tmp_path):
    fa = tmp_path / "adapters.fa"
    fa.write_text(">mine some description\nACGTACGTAC\nGGTT\n>other\nTTTTCCCCAAAA\n")
    assert load_adapters(str(fa), use_builtin=False) == {"mine": "ACGTACGTACGGTT", "other": "TTTTCCCCAAAA"}
    both = load_adapters(str(fa))
    assert set(BUILTIN_ADAPTERS) < set(both) and "mine" in both


### extension ###

def test_extend_by_reads_drops_reads_sharing_the_seed_by_chance():
    rng = random.Random(3)
    seed, after = "ACGTTGCAGGA", "TTACGGATCA"
    # 60 reads continue with `after` (then end), 25 more share the seed by chance
    contexts = [("CC" + seed + after, 2) for _ in range(60)]
    contexts += [("CC" + seed + rand_seq(10, rng), 2) for _ in range(25)]
    assert extend_by_reads(contexts, 11, 100) == after


def test_extend_by_reads_stops_where_most_reads_end():
    seed = "ACGTTGCAGGA"
    contexts = [(seed + "TTAC", 0) for _ in range(80)] + [(seed + "TTACGGGG", 0) for _ in range(10)]
    assert extend_by_reads(contexts, 11, 100) == "TTAC"


### variants ###

def test_is_variant():
    # missing a few bases at an end, or with a bit of adapter left on
    assert is_variant(P1[3:], P1, 0.1)
    assert is_variant(P1 + A1[:4], P1, 0.1)
    # either strand
    assert is_variant(revcomp(P1[:-2]), P1, 0.1)
    # a dimer holding a primer isn't a variant of it
    assert not is_variant(DIMER, P1, 0.1)
    assert not is_variant(P2, P1, 0.1)


### counting ###

def test_matcher():
    m = Matcher([P1, DIMER], 0.1)
    # the primer in full, with an error, on the other strand, at an offset
    read = "AC" + P1[:5] + ("A" if P1[5] != "A" else "C") + P1[6:]
    assert m.hits(read) == {0: (2, False, 1)}
    hits = m.hits(revcomp(P1))
    assert list(hits) == [0] and hits[0][1] is True
    # a dimer read holds the dimer and the primer in full
    assert set(m.hits(DIMER)) == {0, 1}
    # a primer missing a few bases is only within the shortest it's within (the primer, not
    # the dimer that also holds it)
    assert m.hits(P1[4:]) == {0: (None, False, 0)}
    # too short to be checked for being within anything
    assert m.hits(P1[:12]) == {}


def test_matcher_prefilter_keeps_what_hits_finds():
    rng = random.Random(9)
    seqs = [rand_seq(n, rng) for n in (12, 20, 25, 60)]
    m = Matcher(seqs, 0.1)
    reads = []
    for s in seqs:
        for _ in range(20):
            piece = mutate(s, rng, rate=0.05)
            reads.append(rand_seq(rng.randint(0, 10), rng) + piece[:50] + rand_seq(rng.randint(0, 10), rng))
    reads += [rand_seq(50, rng) for _ in range(50)]
    hit = 0
    for read in reads:
        got = m.hits(read)
        assert got == m.hits(read, False)
        hit += bool(got)
    # and plenty were placed, so the check above checked something (62 of the 80 pieces at a
    # 5% error rate are within 10% edits, as found with the prefilter off too)
    assert hit == 62


### mate asymmetry ###

def test_mate_asymmetry():
    k = 5
    by_mate = {0: np.zeros(4 ** k, dtype=np.int64), 1: np.zeros(4 ** k, dtype=np.int64)}
    insert, adapter = "ACGGTCATTG", "TTTGCAAC"
    enc, _, _ = encode([insert + adapter])
    codes, _ = window_codes(enc, k)
    enc2, _, _ = encode([revcomp(insert)])
    codes2, _ = window_codes(enc2, k)
    # R1 reads the insert then the adapter, R2 only the insert's reverse complement
    for c in codes:
        by_mate[0][c] += 100
    for c in codes2:
        by_mate[1][c] += 100
    rc_all = rc_codes(np.arange(4 ** k, dtype=np.int64), k)
    asym = [a for a, _, _ in mate_asymmetry(insert + adapter, by_mate, k, rc_all)]
    assert all(a == 0 for a in asym[:len(insert) - k + 1])
    assert all(a == 1 for a in asym[len(insert) - k + 1:])
    assert split_by_mate_symmetry(insert + adapter, by_mate, k, rc_all) == (0, len(insert))


### the whole thing ###

# Expected counts below were taken from a run after checking them against the simulation's
# own draws: 5,994 P1 molecules, 3,053 P2 and 424 dimers. So P1 is in 6,418 pairs (P1 + dimers),
# P2 in 3,477 (P2 + dimers, whose 16 bases of P2 are within its allowed edits; in a few single
# reads errors push that past them), and the dimer in 424.

def test_find_overrep_paired(reads):
    result = find_overrep([str(p) for p in reads])
    assert found_counts(result) == {
        canon(P1): ((6418, 6418), 6418),
        canon(P2): ((3457, 3459), 3477),
        canon(DIMER): ((424, 424), 424),
    }
    assert result["inferred"] == {}
    r1 = result["tallies"]["R1"]
    assert r1["reads"] == 20_000
    # P1 and P2 sit at read starts
    p1 = next(i for i in result["order"] if result["found"][i]["insert"] == canon(P1))
    assert r1["starts"][p1].most_common(1)[0][0] == 0


def test_find_overrep_infers_adapters_from_pairs(reads):
    result = find_overrep([str(p) for p in reads], adapters={})
    assert set(found_counts(result)) == {canon(P1), canon(P2), canon(DIMER)}
    inferred = sorted(result["inferred"].values())
    assert len(inferred) == 2
    assert A1.startswith(inferred[0]) or inferred[0].startswith(A1[:20])
    assert any(a.startswith("AGATCGGAAGAGCACACG") for a in inferred)
    assert any(a.startswith("AGATCGGAAGAGCGTCGTG") for a in inferred)


def test_find_overrep_single_end(reads):
    result = find_overrep([str(reads[0])])
    assert found_counts(result) == {
        canon(P1): ((6418,), None),
        canon(P2): ((3457,), None),
        canon(DIMER): ((424,), None),
    }


def test_write_outputs(reads, tmp_path):
    result = find_overrep([str(p) for p in reads])
    paths = write_outputs(result, tmp_path, "s1-")
    rows = [line.split("\t") for line in open(paths["seqs"]).read().splitlines()]
    header, rows = rows[0], rows[1:]
    assert header[:5] == ["id", "sequence", "length", "R1_reads", "R1_pct_of_reads"]
    assert [r[0] for r in rows] == ["overrep_1", "overrep_2", "overrep_3"]
    assert rows[0][1] == canon(P1) and rows[0][3] == "6418"
    # the dimer is noted as holding both primers
    assert rows[2][header.index("contains")] == "overrep_1,overrep_2"
    fasta = open(paths["fasta"]).read().split()
    assert fasta[:2] == [">overrep_1", canon(P1)]
    starts = [line.split("\t") for line in open(paths["starts"]).read().splitlines()[1:]]
    assert sum(int(c) for i, lab, _, c in starts if i == "overrep_1" and lab == "R1") == result["tallies"]["R1"]["full"][result["order"][0]]
    adapters = {line.split("\t")[0]: line.split("\t") for line in open(paths["adapters"]).read().splitlines()[1:]}
    assert int(adapters["illumina_truseq_r1"][3]) > 9000
    assert any("Overrepresented sequences found:" in line and line.endswith("3") for line in summary_lines(result))


### cli ###

def test_cli(reads, tmp_path):
    out = tmp_path / "cli"
    base = ["bit", "find-overrep-seqs", "-i", str(reads[0]), "-I", str(reads[1]), "-o", str(out)]
    result = run_cli(base)
    assert "Overrepresented sequences found:" in result.stdout
    assert sorted(p.name for p in out.iterdir()) == ["adapters.tsv", "command-execution-info.txt",
                                                     "overrep-seqs.fasta", "overrep-seqs.tsv",
                                                     "start-positions.tsv", "summary.txt"]
    # an existing output dir stops the run without -F
    assert subprocess.run(base, capture_output=True, text=True).returncode != 0
    run_cli(base + ["-F", "-O", "s1-"])
    assert (out / "s1-overrep-seqs.tsv").exists()


@pytest.mark.parametrize("bad, message", [
    (["-k", "13"], "--kmer-size"),
    (["--min-perc", "0"], "--min-perc"),
    (["--min-perc", "101"], "--min-perc"),
    (["--max-diff", "0.6"], "--max-diff"),
    (["--min-fold", "0.5"], "--min-fold"),
])
def test_cli_rejects_bad_values(reads, tmp_path, bad, message):
    result = subprocess.run(["bit", "find-overrep-seqs", "-i", str(reads[0]), "-o", str(tmp_path / "x")] + bad,
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert message in result.stderr


def test_find_overrep_by_position(reads):
    # with --min-perc above the dimer's ~2% of reads, it's still found, since it always
    # starts reads (down to a tenth of min_perc)
    result = find_overrep([str(p) for p in reads], min_perc=5)
    how = {result["found"][i]["insert"]: result["found"][i]["how"] for i in result["order"]}
    assert how == {canon(P1): "overall", canon(P2): "overall", canon(DIMER): "position"}
    # but not below a tenth of it: at 0.3, P2 (~17% of reads) is found by position, the dimer isn't
    result = find_overrep([str(p) for p in reads], min_perc=30)
    how = {result["found"][i]["insert"]: result["found"][i]["how"] for i in result["order"]}
    assert how == {canon(P1): "overall", canon(P2): "position"}


def test_find_overrep_fasta_input(reads, tmp_path):
    fa = tmp_path / "r1.fa"
    with open(reads[0]) as f, open(fa, "w") as out:
        for i, line in enumerate(f):
            if i % 4 == 0:
                out.write(">" + line[1:])
            elif i % 4 == 1:
                out.write(line)
    result = find_overrep([str(fa)])
    assert found_counts(result) == {
        canon(P1): ((6418,), None),
        canon(P2): ((3457,), None),
        canon(DIMER): ((424,), None),
    }


### read input ###

def test_gzipped_input(reads, tmp_path):
    import gzip
    import shutil
    gz = tmp_path / "r1.fq.gz"
    with open(reads[0], "rb") as f, gzip.open(gz, "wb") as out:
        shutil.copyfileobj(f, out)
    assert found_counts(find_overrep([str(gz)])) == found_counts(find_overrep([str(reads[0])]))


def test_bad_inputs(reads, tmp_path):
    not_reads = tmp_path / "not-reads.txt"
    not_reads.write_text("hello\n")
    with pytest.raises(ValueError, match="doesn't look like fasta or fastq"):
        find_overrep([str(not_reads)])
    cut_short = tmp_path / "short.fq"
    cut_short.write_text("@r1\nACGTACGTAC\n+\nIIII\n")
    with pytest.raises(ValueError, match="ends early"):
        find_overrep([str(cut_short)])
    # pairs out of step
    fewer = tmp_path / "fewer.fq"
    with open(reads[1]) as f:
        fewer.write_text("".join(f.readlines()[:400]))
    with pytest.raises(ValueError, match="same number of reads"):
        find_overrep([str(reads[0]), str(fewer)])


def test_low_complexity():
    assert low_complexity("A" * 20)
    assert low_complexity("CA" * 15)
    assert low_complexity("ACG" * 8)
    assert not low_complexity(P1)
    assert not low_complexity(P2)


def test_dominant_kmers_are_candidates_whatever_their_fold(reads):
    # with a fold no k-mer could reach, only those in at least 1% of reads are picked: the two
    # primers and the dimer (~2% of reads), with the same counts as usual
    result = find_overrep([str(p) for p in reads], min_fold=1e12)
    assert found_counts(result) == found_counts(find_overrep([str(p) for p in reads]))


def test_progress_shown_when_verbose(reads, capsys):
    result = find_overrep([str(p) for p in reads], verbose=True)
    out = capsys.readouterr().out
    for n in range(1, 6):
        assert f"pass {n} of 5" in out
    assert "100%" in out
    # and the same results as without
    assert found_counts(result) == found_counts(find_overrep([str(p) for p in reads]))


def test_pair_counts(reads):
    result = find_overrep([str(p) for p in reads])
    pairs, r1, r2 = result["pairs"], result["tallies"]["R1"], result["tallies"]["R2"]
    for i in result["order"]:
        # a pair counts for either_read if either read holds it, for both_reads if both do
        assert pairs["pairs_with"][i] == r1["reads_with"][i] + r2["reads_with"][i] - pairs["pairs_both"][i]
        assert pairs["pairs_both"][i] <= min(r1["reads_with"][i], r2["reads_with"][i])
    assert pairs["explained_both"] <= pairs["explained"]


def test_max_reads_for_discovery(reads):
    # finding from the first 2,000 reads of each file gives the same sequences, counted
    # across all reads, as finding from all of them
    full = find_overrep([str(p) for p in reads])
    capped = find_overrep([str(p) for p in reads], max_reads_for_discovery=2_000)
    assert found_counts(capped) == found_counts(full)
    # the cap holds the discovery passes to that many reads per file
    import bit.modules.find_overrep as fo
    seen = []
    orig = fo.native.count_kmers
    fo.native.count_kmers = lambda *a, **k: (seen.append(k.get("max_reads")), orig(*a, **k))[1]
    try:
        r = find_overrep([str(p) for p in reads], max_reads_for_discovery=2_000)
    finally:
        fo.native.count_kmers = orig
    assert seen == [2_000, 2_000]
    assert r["tallies"]["R1"]["reads"] == 20_000   # counting still used all reads


def test_long_reads_skip_the_positional_look(tmp_path, capsys):
    rng = random.Random(4)
    genome = rand_seq(100_000, rng)
    fq = tmp_path / "long.fq"
    with open(fq, "w") as out:
        for i in range(500):
            s = rng.randint(0, len(genome) - 1500)
            out.write(f"@r{i}\n{genome[s:s + 1500]}\n+\n{'I' * 1500}\n")
    # (only the decision is of interest; stopping right after it keeps this quick, since with
    # thresholds sized for short reads, long reads make for hundreds of candidates to extend)
    import bit.modules.find_overrep as fo

    class Stop(Exception):
        pass

    def bail(*a, **k):
        raise Stop

    orig = fo.extend_all
    fo.extend_all = bail
    try:
        with pytest.raises(Stop):
            find_overrep([str(fq)], verbose=True)
    finally:
        fo.extend_all = orig
    out = capsys.readouterr().out
    assert "pass 3 of 5): skipped, reads are longer than 1,000 bases" in out


def test_cli_help_menus():
    basic = run_cli(["bit", "find-overrep-seqs", "-h"]).stdout
    detailed = subprocess.run(["bit", "find-overrep-seqs", "-H"], capture_output=True, text=True)
    assert detailed.returncode == 0
    detailed = detailed.stderr
    fine = ["--min-perc", "--min-fold", "--max-diff", "--kmer-size", "--max-len", "--adapters",
            "--no-builtin-adapter-trim", "--max-reads-for-discovery"]
    # fine-tuning parameters, and how it works, only in the detailed menu
    assert not any(f in basic for f in fine) and "How it works" not in basic
    assert all(f in detailed for f in fine) and "How it works" in detailed
    assert "Fine-tuning Parameters" in detailed


def test_cli_hidden_parameters_still_parse(reads, tmp_path):
    from bit.cli.find_overrep import build_parser
    args = build_parser().parse_args(["-i", "x", "-m", "5", "-k", "9", "--no-builtin-adapter-trim", "-M", "0"])
    assert (args.min_perc, args.kmer_size, args.no_builtin_adapter_trim, args.max_reads_for_discovery) == (5, 9, True, 0)
    # and logged with the rest of the command
    out = tmp_path / "cli-hidden"
    run_cli(["bit", "find-overrep-seqs", "-i", str(reads[0]), "-o", str(out), "-m", "5"])
    assert "--min-perc 5.0" in (out / "command-execution-info.txt").read_text()
