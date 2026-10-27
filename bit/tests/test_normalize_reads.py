import gzip
import random
import statistics

import pytest # type: ignore

from bit import _native as native
from bit.modules.normalize_reads import default_sketch_bytes, normalize, output_paths
from bit.tests.utils import run_cli


_COMP = str.maketrans("ACGT", "TGCA")
READ = 150


def rc(s):
    return s.translate(_COMP)[::-1]


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def with_errors(seq, rate, rng):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def sample_reads(genome, depth, name, rng, error_rate=0.0):
    """ single reads at `depth`, random positions and strands, as (name, seq, qual) """
    reads = []
    for i in range(int(depth * len(genome) / READ)):
        s = rng.randrange(len(genome) - READ)
        seq = genome[s:s + READ]
        if rng.random() < 0.5:
            seq = rc(seq)
        reads.append((f"{name}_{i}_{s}", with_errors(seq, error_rate, rng), "I" * READ))
    return reads


def sample_pairs(genome, depth, name, rng, insert=350):
    r1s, r2s = [], []
    for i in range(int(depth * len(genome) / (2 * READ))):
        s = rng.randrange(len(genome) - insert)
        frag = genome[s:s + insert]
        if rng.random() < 0.5:
            frag = rc(frag)
        r1s.append((f"{name}_{i}/1", frag[:READ], "I" * READ))
        r2s.append((f"{name}_{i}/2", rc(frag)[:READ], "I" * READ))
    return r1s, r2s


def shuffled_pairs(rng, *sets):
    r1 = [r for s in sets for r in s[0]]
    r2 = [r for s in sets for r in s[1]]
    order = list(range(len(r1)))
    rng.shuffle(order)
    return [r1[i] for i in order], [r2[i] for i in order]


def variant(genome, every, rng):
    """ (genome with a SNP every `every` bp, the SNP positions) """
    out, snps = list(genome), list(range(30, len(genome) - 30, every))
    for p in snps:
        out[p] = rng.choice([b for b in "ACGT" if b != out[p]])
    return "".join(out), snps


def write_fastq(path, recs):
    with open(path, "w") as f:
        for name, seq, qual in recs:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_names(path):
    with gzip.open(path, "rt") as f:
        return [l[1:].strip() for i, l in enumerate(f) if i % 4 == 0]


def read_seqs(path):
    with gzip.open(path, "rt") as f:
        return [l.strip() for i, l in enumerate(f) if i % 4 == 1]


def run(tmp_path, recs, recs2=None, **kwargs):
    tmp_path.mkdir(parents=True, exist_ok=True)
    write_fastq(tmp_path / "r1.fq", recs)
    if recs2 is not None:
        write_fastq(tmp_path / "r2.fq", recs2)
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    return normalize(str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq") if recs2 is not None else None,
                     str(out), show_progress=False, **kwargs)


def depth_by_name(path, glen, prefix):
    return sum(READ for n in read_names(path) if n.startswith(prefix)) / glen


def covered(path, genome):
    """ fraction of the genome's positions covered by the kept reads (error-free reads) """
    hit = bytearray(len(genome))
    for name in read_names(path):
        s = int(name.split("_")[-1])
        for p in range(s, s + READ):
            hit[p] = 1
    return sum(hit) / len(genome)


def minor_allele_reads(paths, minor, snps):
    """ kept reads (either mate) holding the minor allele, at each SNP """
    seqs = read_seqs(paths["r1"]) + read_seqs(paths["r2"]) if "r1" in paths else read_seqs(paths["reads"])
    probes = [minor[p - 10:p + 11] for p in snps]
    return [sum(1 for s in seqs if pr in s or rc(pr) in s) for pr in probes]


@pytest.fixture
def rng():
    return random.Random(11)


def test_deep_genome_capped_and_fully_covered(tmp_path, rng):
    genome = rand_seq(5000, rng)
    reads = sample_reads(genome, 300, "g", rng)
    rng.shuffle(reads)
    counts, paths = run(tmp_path, reads, target_depth=20)
    depth = depth_by_name(paths["reads"], len(genome), "g_")
    assert 15 <= depth <= 45
    assert covered(paths["reads"], genome) > 0.99
    assert counts.reads_or_pairs_in == len(reads)
    assert counts.reads_or_pairs_out == len(read_names(paths["reads"]))


def test_rare_genome_kept_whole(tmp_path, rng):
    deep, rare = rand_seq(5000, rng), rand_seq(5000, rng)
    reads = sample_reads(deep, 300, "deep", rng) + sample_reads(rare, 3, "rare", rng)
    rng.shuffle(reads)
    _, paths = run(tmp_path, reads, target_depth=20)
    n_rare = sum(1 for r in reads if r[0].startswith("rare_"))
    assert sum(1 for n in read_names(paths["reads"]) if n.startswith("rare_")) == n_rare


def test_rare_part_of_a_read_keeps_it(tmp_path, rng):
    # a rare replicon that is mostly a common backbone: its reads are mostly deep sequence,
    # so their median depth is high, but their own stretch keeps them
    backbone, own = rand_seq(4000, rng), rand_seq(1000, rng)
    common = backbone + rand_seq(1000, rng)
    rare = backbone + own
    reads = sample_reads(common, 400, "common", rng) + sample_reads(rare, 8, "rare", rng)
    rng.shuffle(reads)
    _, paths = run(tmp_path, reads, target_depth=20)
    # every rare read reaching 40+ bp into its own stretch is kept
    wanted = [r for r in reads if r[0].startswith("rare_") and int(r[0].split("_")[-1]) + READ >= 4040]
    kept_names = set(read_names(paths["reads"]))
    assert wanted and all(r[0] in kept_names for r in wanted)


def test_minor_variant_kept(tmp_path, rng):
    # a variant at a tenth of the depth, one SNP every 100 bp: kept up to about the cap
    major = rand_seq(5000, rng)
    minor, _ = variant(major, 100, rng)
    reads = sample_reads(major, 400, "major", rng) + sample_reads(minor, 40, "minor", rng)
    rng.shuffle(reads)
    _, paths = run(tmp_path, reads, target_depth=20)
    assert depth_by_name(paths["reads"], 5000, "minor_") >= 15


def test_minor_branches_keep_sparse_variants(tmp_path, rng):
    # pairs, one SNP every 300 bp: a pair holding one SNP has under 10% of its anchors on the
    # variant, so it's kept for the minor branch
    major = rand_seq(6000, rng)
    minor, snps = variant(major, 300, rng)
    r1, r2 = shuffled_pairs(rng, sample_pairs(major, 400, "major", rng),
                            sample_pairs(minor, 40, "minor", rng))
    _, on = run(tmp_path / "on", r1, r2, target_depth=20)
    _, off = run(tmp_path / "off", r1, r2, target_depth=20, minor_branches=False)
    with_branches = minor_allele_reads(on, minor, snps)
    without = minor_allele_reads(off, minor, snps)
    assert statistics.median(with_branches) > statistics.median(without)


def test_every_sparse_variant_kept(tmp_path, rng):
    # judged in one pass, a lone SNP's anchors look like an error's until enough reads with it
    # have come by, and its minor allele could be lost altogether. Counted first, it isn't
    major = rand_seq(6000, rng)
    minor, snps = variant(major, 600, rng)
    r1, r2 = shuffled_pairs(rng, sample_pairs(major, 400, "major", rng),
                            sample_pairs(minor, 40, "minor", rng))
    _, paths = run(tmp_path, r1, r2, target_depth=20)
    kept = minor_allele_reads(paths, minor, snps)
    seqs_in = [r[1] for r in r1 + r2]
    had = [sum(1 for s in seqs_in if pr in s or rc(pr) in s) for pr in (minor[p - 10:p + 11] for p in snps)]
    # every SNP keeps its minor allele: about half the target (the minor branch's cap) where the
    # input had plenty, and most of it where it had few (the first SNP, near the genome's end)
    # (a one-pass version, Oct 2026, same reads: one SNP of the ten kept none)
    assert all(k >= min(5, h - 1) for k, h in zip(kept, had))


def test_same_reads_kept_whatever_the_order(tmp_path, rng):
    genome = rand_seq(5000, rng)
    reads = sample_reads(genome, 300, "g", rng)
    by_pos = sorted(reads, key=lambda r: int(r[0].split("_")[-1]))
    shuffled = reads[:]
    rng.shuffle(shuffled)
    _, p_sorted = run(tmp_path / "a", by_pos, target_depth=20)
    _, p_shuf = run(tmp_path / "b", shuffled, target_depth=20)
    assert sorted(read_names(p_sorted["reads"])) == sorted(read_names(p_shuf["reads"]))


def test_reads_with_errors_not_favored(tmp_path, rng):
    genome = rand_seq(5000, rng)
    clean = sample_reads(genome, 300, "clean", rng)
    noisy = sample_reads(genome, 100, "noisy", rng, error_rate=0.01)
    reads = clean + noisy
    rng.shuffle(reads)
    _, paths = run(tmp_path, reads, target_depth=20)
    names = read_names(paths["reads"])
    share_in = len(noisy) / len(reads)
    share_out = sum(1 for n in names if n.startswith("noisy_")) / len(names)
    assert share_out <= 1.25 * share_in


def test_deep_isolate_errors_not_kept(tmp_path, rng):
    # 1,000x with Illumina-like errors: what's kept is no more error-rich than the input, and
    # the errors are counted as such
    genome = rand_seq(3000, rng)
    reads = sample_reads(genome, 1000, "g", rng, error_rate=0.003)
    counts, paths = run(tmp_path, reads, target_depth=50)

    def mismatches(name, seq):
        s = int(name.split("_")[-1])
        ref = genome[s:s + READ]
        return min(sum(a != b for a, b in zip(seq, ref)), sum(a != b for a, b in zip(rc(seq), ref)))

    rate_in = sum(mismatches(n, s) for n, s, _ in reads) / (READ * len(reads))
    kept = list(zip(read_names(paths["reads"]), read_seqs(paths["reads"])))
    rate_out = sum(mismatches(n, s) for n, s in kept) / (READ * len(kept))
    assert rate_out <= 1.2 * rate_in
    assert 0.02 < counts.error_share < 0.2
    assert len(kept) < 0.15 * len(reads)


def test_pairs_kept_together_in_order(tmp_path, rng):
    genome = rand_seq(5000, rng)
    r1s, r2s = shuffled_pairs(rng, sample_pairs(genome, 200, "g", rng))
    counts, paths = run(tmp_path, r1s, r2s, target_depth=20)
    n1, n2 = read_names(paths["r1"]), read_names(paths["r2"])
    assert [n[:-2] for n in n1] == [n[:-2] for n in n2]
    assert 0 < len(n1) < len(r1s) and counts.reads_or_pairs_out == len(n1)
    pos = {r[0]: i for i, r in enumerate(r1s)}
    assert [pos[n] for n in n1] == sorted(pos[n] for n in n1)


def test_records_written_as_they_were(tmp_path, rng):
    genome = rand_seq(3000, rng)
    recs = [(f"{n} some comment", s, "".join(rng.choice("#5?I") for _ in s))
            for n, s, _ in sample_reads(genome, 5, "g", rng)]
    _, paths = run(tmp_path, recs)
    with gzip.open(paths["reads"], "rt") as f:
        out = f.read()
    with open(tmp_path / "r1.fq") as f:
        assert out == f.read()


def test_reads_without_anchors_dropped(tmp_path, rng):
    genome = rand_seq(3000, rng)
    reads = sample_reads(genome, 5, "g", rng) + [("short", "ACGTACGT", "IIIIIIII"),
                                                 ("ns", "N" * 150, "I" * 150)]
    counts, paths = run(tmp_path, reads, target_depth=20)
    names = read_names(paths["reads"])
    assert "short" not in names and "ns" not in names
    assert counts.no_anchors == 2


def test_fasta_in_fasta_out(tmp_path, rng):
    genome = rand_seq(3000, rng)
    with open(tmp_path / "r.fa", "w") as f:
        for name, seq, _ in sample_reads(genome, 5, "g", rng):
            f.write(f">{name}\n{seq[:70]}\n{seq[70:]}\n")
    (tmp_path / "out").mkdir()
    _, paths = normalize(str(tmp_path / "r.fa"), None, str(tmp_path / "out"), show_progress=False)
    assert paths["reads"].endswith("normalized.fasta.gz")
    with gzip.open(paths["reads"], "rt") as f, open(tmp_path / "r.fa") as g:
        assert f.read() == g.read()


def test_depth_table(tmp_path, rng):
    genome = rand_seq(3000, rng)
    reads = sample_reads(genome, 60, "g", rng)
    counts, paths = run(tmp_path, reads, target_depth=20, depth_table=True)
    with gzip.open(paths["depth_table"], "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    assert rows[0] == ["read", "depth"]
    assert [r[0] for r in rows[1:]] == read_names(paths["reads"])
    # the whole input's depth (60x of 150 bp reads, ~0.84 of positions starting a 25-mer)
    assert 35 <= statistics.median(int(r[1]) for r in rows[1:]) <= 65


def test_no_depth_table_unless_asked(tmp_path, rng):
    _, paths = run(tmp_path, sample_reads(rand_seq(3000, rng), 5, "g", rng))
    assert "depth_table" not in paths
    assert not (tmp_path / "out" / "depths.tsv.gz").exists()


def test_mismatched_mates_stop_the_run(tmp_path, rng):
    genome = rand_seq(3000, rng)
    r1s, r2s = sample_pairs(genome, 5, "g", rng)
    r2s[3] = ("other/2", r2s[3][1], r2s[3][2])
    with pytest.raises(SystemExit):
        run(tmp_path, r1s, r2s)


def test_read_2_short_stops_the_run(tmp_path, rng):
    genome = rand_seq(3000, rng)
    r1s, r2s = sample_pairs(genome, 5, "g", rng)
    with pytest.raises(SystemExit):
        run(tmp_path, r1s, r2s[:-1])


def test_empty_input(tmp_path):
    counts, paths = run(tmp_path, [])
    assert counts.reads_or_pairs_in == 0 and read_names(paths["reads"]) == []


def test_summary_rows(tmp_path, rng):
    counts, paths = run(tmp_path, sample_reads(rand_seq(3000, rng), 20, "g", rng))
    rows = dict(l.rstrip("\n").split("\t") for l in open(paths["summary"]))
    assert rows["measure"] == "value" and int(rows["reads_or_pairs_in"]) == counts.reads_or_pairs_in
    assert "sketch_fill" in rows and "median_depth_in" in rows


def test_counting_pass_counts(tmp_path, rng):
    genome = rand_seq(3000, rng)
    reads = sample_reads(genome, 50, "g", rng)
    write_fastq(tmp_path / "r.fq", reads)
    counted = native.count_anchors(str(tmp_path / "r.fq"), None, 25, 2, 1 << 20)
    assert counted.reads == len(reads)
    # each anchor of a stretch is in about (reads starting within 150 - 25 bp before it)
    mid = counted.counts_of(genome[1400:1600])
    assert 30 <= statistics.median(mid) <= 55


def test_native_rejects_bad_settings(tmp_path):
    write_fastq(tmp_path / "r.fq", [("a", "ACGT" * 40, "I" * 160)])
    with pytest.raises(ValueError):
        native.count_anchors(str(tmp_path / "r.fq"), None, 40, 2, 1 << 20)
    counted = native.count_anchors(str(tmp_path / "r.fq"), None, 25, 2, 1 << 20)
    with pytest.raises(ValueError):
        counted.normalize(str(tmp_path / "r.fq"), None, str(tmp_path / "o.fq.gz"), None, 0)


def test_default_sketch_bytes(tmp_path):
    p = tmp_path / "r.fq"
    p.write_text("@a\nACGT\n+\nIIII\n")
    assert default_sketch_bytes([str(p)]) == 256 * 2**20


def test_output_paths_prefix():
    p = output_paths("d", "s1-", paired=True, fastq=False)
    assert p["r1"].endswith("d/s1-normalized_R1.fasta.gz") and p["summary"].endswith("d/s1-summary.tsv")


def test_cli(tmp_path, rng):
    genome = rand_seq(3000, rng)
    r1s, r2s = sample_pairs(genome, 100, "g", rng)
    write_fastq(tmp_path / "r1.fq", r1s)
    write_fastq(tmp_path / "r2.fq", r2s)
    out = tmp_path / "cli"
    result = run_cli(["bit", "normalize-reads", "-i", str(tmp_path / "r1.fq"), "-I", str(tmp_path / "r2.fq"),
                      "-o", str(out), "-d", "20", "-j", "1", "--depth-table", "-m", "0.01"])
    assert "pairs kept" in result.stdout and "judged sequencing errors" in result.stdout
    assert sorted(p.name for p in out.iterdir()) == ["command-execution-info.txt", "depths.tsv.gz",
                                                     "normalized_R1.fastq.gz", "normalized_R2.fastq.gz",
                                                     "summary.tsv"]


def test_cli_warns_when_counting_runs_out_of_room(tmp_path, rng):
    write_fastq(tmp_path / "r.fq", sample_reads(rand_seq(20000, rng), 30, "g", rng))
    result = run_cli(["bit", "normalize-reads", "-i", str(tmp_path / "r.fq"), "-o", str(tmp_path / "o"),
                      "-j", "1", "-m", "0.00001"])
    assert "--memory" in result.stdout
    summary = result.stdout.index("Summary written to")
    assert result.stdout.index("--memory") > summary


def test_cli_rejects_bad_kmer_size(tmp_path):
    import subprocess
    write_fastq(tmp_path / "r.fq", [("a", "ACGT" * 40, "I" * 160)])
    result = subprocess.run(["bit", "normalize-reads", "-i", str(tmp_path / "r.fq"), "-o",
                             str(tmp_path / "o"), "-k", "40"], capture_output=True, text=True)
    assert result.returncode != 0 and "kmer-size" in result.stderr


def test_cli_detailed_help_shows_fine_tuning():
    result = run_cli(["bit", "normalize-reads", "-H"])
    assert "--novel-fraction" in result.stderr and "HOW IT WORKS" in result.stderr
