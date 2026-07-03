import os
from argparse import Namespace
import shutil
import gzip
import pytest # type: ignore
from bit.modules.gen_reads import parse_proportions_file
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
from bit.modules.gen_reads import (preflight_checks,
                                   gen_reads,
                                   get_proportions,
                                   compute_reads_from_coverage,
                                   extract_subsequence,
                                   parse_coverages_input,
                                   apportion_units,
                                   compute_id_offsets,
                                   distribute_units_across_contigs)


test_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")


def _base_args(**kw):
    d = dict(input_fastas=[], proportions_file=None, coverage=None,
             read_length=10, fragment_size=50, type="paired-end")
    d.update(kw)
    return Namespace(**d)


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _read_fasta_seq(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()


def _count_reads(path):
    with open(path) as f:
        return sum(1 for _ in f) // 4


def _fastq_map(path):
    m = {}
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    for i in range(0, len(lines), 4):
        m[lines[i][1:]] = lines[i + 1]  # strip leading '@'
    return m


def test_gen_reads(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "perfect-reads"),
        "-n", "2",
        "-r", "10",
        "-s", "9"
    ]

    run_cli(cmd)

    R1_path = tmp_path / "perfect-reads_R1.fastq.gz"
    R2_path = tmp_path / "perfect-reads_R2.fastq.gz"
    assert R1_path.exists(), f"R1 file not found at {R1_path}"
    assert R2_path.exists(), f"R2 file not found at {R2_path}"

    with gzip.open(R1_path, 'rt') as f:
        observed = f.read().splitlines()

    expected = [
        "@r1/1 contig=partial-NC_003131.1;start=567;end=577;strand=+",
        "GGAATGCCTG",
        "+",
        "IIIIIIIIII"
    ]

    assert observed == expected, f"R1 content does not match expected:\n{observed}\nExpected:\n{expected}"


def test_parse_proportions_file_logic(tmp_path):

    prop_file = tmp_path / "proportions.tsv"

    # test case if the user provides abundances instead of proportions summing to 1
    prop_file.write_text("file1.fasta\t3\nfile2.fasta\t1\n")

    # total_proportion will be 4.0; file1 should become 0.75, file2 should become 0.25
    observed = parse_proportions_file(str(prop_file), ["file1.fasta", "file2.fasta"])

    assert observed["file1.fasta"] == 0.75
    assert observed["file2.fasta"] == 0.25
    assert len(observed) == 2


def test_parse_proportions_equal_split():

    input_fastas = ["a.fasta", "b.fasta", "c.fasta", "d.fasta"]

    observed = parse_proportions_file(None, input_fastas)

    assert observed["a.fasta"] == 0.25
    assert observed["b.fasta"] == 0.25
    assert observed["c.fasta"] == 0.25
    assert observed["d.fasta"] == 0.25
    assert sum(observed.values()) == pytest.approx(1.0)


def _write_fasta(path, seq_id, seq):
    with open(path, 'w') as f:
        f.write(f">{seq_id}\n")
        f.write(seq + "\n")


def _collect_forward_reads(fastq_path):
    seqs = []
    with open(fastq_path, 'r') as f:
        lines = [l.rstrip('\n') for l in f]
    seqs = lines[1::4]
    return seqs


def _reverse_complement(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def test_wrapped_fragments_generated_with_circularize(tmp_path):

    fasta = tmp_path / "contig.fasta"
    # distinct halves so wrapped fragments are easy to recognize (A.. then C..)
    seq = "A" * 5 + "C" * 5
    _write_fasta(fasta, "contig1", seq)

    out_prefix = tmp_path / "out_circ"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage=None,
        seed=9,
        num_reads=20,
        fragment_size=6,
        fragment_size_range=10,
        read_length=4,
        output_prefix=str(out_prefix),
        circularize=True,
        include_Ns=False,
        per_read_tsv=False,
        type="paired-end",
    )

    gen_reads(args, get_proportions(args))

    forward_fastq = str(out_prefix) + "_R1.fastq"
    assert os.path.exists(forward_fastq)
    seqs = _collect_forward_reads(forward_fastq)

    # compute all possible wrapped forward read sequences for this contig
    seq_len = len(seq)
    frag_len = min(args.fragment_size, seq_len)
    wrapped_fw = set()
    for start in range(seq_len):
        end = start + frag_len
        if end > seq_len:
            fragment = seq[start:] + seq[: end - seq_len]
            wrapped_fw.add(fragment[: args.read_length])

    # fragments are now drawn from either strand ~50/50, so a wrapped forward read
    # may appear as the reverse complement of the forward-orientation fragment
    wrapped_all = wrapped_fw | {_reverse_complement(w) for w in wrapped_fw}

    # ensure at least one wrapped forward read appears
    assert any(s in wrapped_all for s in seqs), "No wrapped fragments found with circularize=True"


def test_gen_reads_single_end(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "se-reads"),
        "-n", "1",
        "-r", "10",
        "-s", "1",
        "--type", "single-end",
    ]

    run_cli(cmd)

    se_path = tmp_path / "se-reads.fastq.gz"
    assert se_path.exists(), f"Single-end file not found at {se_path}"

    # paired-end files should NOT exist
    assert not (tmp_path / "se-reads_R1.fastq.gz").exists()
    assert not (tmp_path / "se-reads_R2.fastq.gz").exists()

    with gzip.open(se_path, 'rt') as f:
        lines = f.read().splitlines()

    # should have 4 lines per read (1 read pair → 1 single-end read)
    assert len(lines) == 4, f"Expected 4 lines, got {len(lines)}"
    # header should NOT have /1 or /2 suffix
    assert not lines[0].endswith("/1") and not lines[0].endswith("/2")
    # read length should match requested
    assert len(lines[1]) == 10


def test_simulate_single_end_reads_unit(tmp_path):

    fasta = tmp_path / "contig.fasta"
    _write_fasta(fasta, "contig1", "ACGTACGTACGTACGT")

    out_prefix = tmp_path / "se_out"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage=None,
        seed=42,
        num_reads=5,
        read_length=8,
        output_prefix=str(out_prefix),
        circularize=False,
        type="single-end",
        include_Ns=False,
        per_read_tsv=False,
    )

    gen_reads(args, get_proportions(args))

    fastq_path = str(out_prefix) + ".fastq"
    assert os.path.exists(fastq_path), "Single-end FASTQ not created"

    with open(fastq_path, 'r') as f:
        lines = [l.rstrip('\n') for l in f]

    # 5 reads × 4 lines each = 20 lines
    assert len(lines) == 20, f"Expected 20 lines, got {len(lines)}"

    for i in range(5):
        header = lines[i * 4]
        seq = lines[i * 4 + 1]
        plus = lines[i * 4 + 2]
        qual = lines[i * 4 + 3]
        assert header.startswith("@")
        assert not header.endswith("/1") and not header.endswith("/2")
        assert len(seq) == 8
        assert plus == "+"
        assert len(qual) == 8


def test_long_reads_variable_length(tmp_path):

    fasta = tmp_path / "contig.fasta"
    # 2000-bp contig so read lengths aren't clamped by seq_length
    _write_fasta(fasta, "contig1", "ACGT" * 500)

    out_prefix = tmp_path / "long_out"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage=None,
        seed=42,
        num_reads=100,
        read_length=1000,
        output_prefix=str(out_prefix),
        circularize=False,
        type="long",
        long_read_length_range=50,
        include_Ns=False,
        per_read_tsv=False,
    )

    gen_reads(args, get_proportions(args))

    fastq_path = str(out_prefix) + ".fastq"
    assert os.path.exists(fastq_path)

    with open(fastq_path, 'r') as f:
        lines = [l.rstrip('\n') for l in f]

    read_lengths = [len(lines[i * 4 + 1]) for i in range(len(lines) // 4)]

    # all reads should be within 500-1500 (50% of 1000)
    assert all(500 <= rl <= 1500 for rl in read_lengths), \
        f"Read lengths outside expected range: min={min(read_lengths)}, max={max(read_lengths)}"

    # quality scores should match read length for each read
    for i in range(len(lines) // 4):
        seq = lines[i * 4 + 1]
        qual = lines[i * 4 + 3]
        assert len(qual) == len(seq), f"Quality length {len(qual)} != read length {len(seq)}"

    # reads should NOT all be the same length (variable lengths)
    assert len(set(read_lengths)) > 1, "All reads have the same length; expected variable lengths"


def test_long_reads_cli(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "long-reads"),
        "-n", "10",
        "-r", "100",
        "--type", "long",
        "-s", "1",
    ]

    run_cli(cmd)

    # should produce single-end output (implied by --long)
    se_path = tmp_path / "long-reads.fastq.gz"
    assert se_path.exists(), f"Single-end file not found at {se_path}"
    assert not (tmp_path / "long-reads_R1.fastq.gz").exists()
    assert not (tmp_path / "long-reads_R2.fastq.gz").exists()

    with gzip.open(se_path, 'rt') as f:
        lines = f.read().splitlines()

    read_lengths = [len(lines[i * 4 + 1]) for i in range(len(lines) // 4)]

    # with default 50% range and read_length=100, expect reads from 50-150
    assert all(50 <= rl <= 150 for rl in read_lengths), \
        f"Read lengths outside expected range: min={min(read_lengths)}, max={max(read_lengths)}"


def test_extract_subsequence_skips_Ns():

    import random
    random.seed(0)

    # sequence with Ns only in the middle
    seq = "ACGTACGT" + "NNNNNNNN" + "ACGTACGT"
    seq_length = len(seq)
    read_len = 8

    # generate many reads; none should contain N when include_Ns=False
    for _ in range(50):
        subseq, _ = extract_subsequence(seq, seq_length, read_len, circularize=False, include_Ns=False)
        assert 'N' not in subseq, f"Found N in read when include_Ns=False: {subseq}"


def test_extract_subsequence_includes_Ns_when_flagged():

    import random
    random.seed(0)

    # sequence that is entirely Ns — reads must contain Ns
    seq = "N" * 20
    seq_length = len(seq)
    read_len = 8

    subseq, _ = extract_subsequence(seq, seq_length, read_len, circularize=False, include_Ns=True)
    assert 'N' in subseq, "Expected Ns in read when include_Ns=True"


def test_gen_paired_reads_skips_Ns(tmp_path):

    fasta = tmp_path / "contig.fasta"
    # 20-char clean regions on each side so all fragment sizes (7-9) fit cleanly
    seq = "ACGTACGTACGTACGTACGT" + "NNNNNNNNNN" + "ACGTACGTACGTACGTACGT"
    _write_fasta(fasta, "contig1", seq)

    out_prefix = tmp_path / "no_ns"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage=None,
        seed=42,
        num_reads=20,
        fragment_size=8,
        fragment_size_range=10,
        read_length=4,
        output_prefix=str(out_prefix),
        circularize=False,
        include_Ns=False,
        per_read_tsv=False,
        type="paired-end",
    )

    gen_reads(args, get_proportions(args))

    for suffix in ("_R1.fastq", "_R2.fastq"):
        with open(str(out_prefix) + suffix, 'r') as f:
            lines = [l.rstrip('\n') for l in f]
        seqs = lines[1::4]
        for s in seqs:
            assert 'N' not in s, f"Found N in {suffix} read when include_Ns=False: {s}"


def test_parse_coverages_input_single_value():

    input_fastas = ["a.fasta", "b.fasta", "c.fasta"]
    result = parse_coverages_input("50", input_fastas)

    assert result == {"a.fasta": 50.0, "b.fasta": 50.0, "c.fasta": 50.0}


def test_parse_coverages_input_float_value():

    input_fastas = ["a.fasta"]
    result = parse_coverages_input("10.5", input_fastas)

    assert result == {"a.fasta": 10.5}


def test_parse_coverages_input_tsv(tmp_path):

    cov_file = tmp_path / "coverages.tsv"
    cov_file.write_text("a.fasta\t100\nb.fasta\t50\n")

    result = parse_coverages_input(str(cov_file), ["a.fasta", "b.fasta"])

    assert result == {"a.fasta": 100.0, "b.fasta": 50.0}


def test_coverage_mode_paired_end(tmp_path):

    fasta = tmp_path / "genome.fasta"
    # 1000-bp genome, 100x coverage, fragment_size=500, read_length=150
    # bases_per_fragment = min(2*150, 500) = 300
    # fragments needed = 100 * 1000 / 300 = 333 (rounded)
    # reads = 2 * 333 = 666
    _write_fasta(fasta, "contig1", "ACGT" * 250)

    out_prefix = tmp_path / "cov_pe"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage="100",
        seed=42,
        num_reads=1000000,  # should be overridden
        fragment_size=500,
        fragment_size_range=10,
        read_length=150,
        output_prefix=str(out_prefix),
        circularize=False,
        include_Ns=False,
        type="paired-end",
        per_read_tsv=False,
    )

    proportions = get_proportions(args)

    assert args.num_reads == 666, f"Expected 666 reads, got {args.num_reads}"
    assert proportions[str(fasta)] == pytest.approx(1.0)


def test_coverage_mode_single_end(tmp_path):

    fasta = tmp_path / "genome.fasta"
    # 1000-bp genome, 50x coverage, read_length=100
    # reads needed = 50 * 1000 / 100 = 500
    _write_fasta(fasta, "contig1", "ACGT" * 250)

    out_prefix = tmp_path / "cov_se"
    args = Namespace(
        input_fastas=[str(fasta)],
        proportions_file=None,
        coverage="50",
        seed=42,
        num_reads=1000000,  # should be overridden
        read_length=100,
        output_prefix=str(out_prefix),
        circularize=False,
        include_Ns=False,
        type="single-end",
        per_read_tsv=False,
    )

    proportions = get_proportions(args)

    assert args.num_reads == 500, f"Expected 500 reads, got {args.num_reads}"
    assert proportions[str(fasta)] == pytest.approx(1.0)


def test_coverage_mode_multiple_genomes(tmp_path):

    fasta_a = tmp_path / "big.fasta"
    fasta_b = tmp_path / "small.fasta"
    _write_fasta(fasta_a, "big", "ACGT" * 500)    # 2000 bp
    _write_fasta(fasta_b, "small", "ACGT" * 125)   # 500 bp

    args = Namespace(
        input_fastas=[str(fasta_a), str(fasta_b)],
        proportions_file=None,
        coverage="100",
        seed=42,
        num_reads=1000000,
        read_length=100,
        output_prefix=str(tmp_path / "cov_multi"),
        circularize=False,
        include_Ns=False,
        type="single-end",
        per_read_tsv=False,
    )

    proportions = get_proportions(args)

    # big genome: 100 * 2000 / 100 = 2000 reads
    # small genome: 100 * 500 / 100 = 500 reads
    # total = 2500
    assert args.num_reads == 2500, f"Expected 2500 reads, got {args.num_reads}"
    assert proportions[str(fasta_a)] == pytest.approx(0.8)
    assert proportions[str(fasta_b)] == pytest.approx(0.2)


def test_coverage_cli_paired_end(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "cov-reads"),
        "-c", "10",
        "-r", "10",
        "-f", "20",
        "-s", "1",
    ]

    run_cli(cmd)

    R1_path = tmp_path / "cov-reads_R1.fastq.gz"
    R2_path = tmp_path / "cov-reads_R2.fastq.gz"
    assert R1_path.exists(), f"R1 file not found at {R1_path}"
    assert R2_path.exists(), f"R2 file not found at {R2_path}"

    # standalone runs also write a per-genome summary
    summary = tmp_path / "cov-reads-per-genome-summary.tsv"
    assert summary.exists(), f"summary not found at {summary}"


def test_per_genome_summary_columns_and_values(tmp_path):
    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "sum-reads"),
        "-c", "30",
        "-r", "100",
        "-s", "1",
    ]
    run_cli(cmd)

    summary = tmp_path / "sum-reads-per-genome-summary.tsv"
    assert summary.exists()
    with open(summary) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]

    assert header == ["input_fasta", "genome_size", "reads_generated",
                      "mean_coverage", "detection"]
    assert len(rows) == 1                       # single input genome -> one row
    row = dict(zip(header, rows[0]))
    assert int(row["genome_size"]) > 0
    assert int(row["reads_generated"]) > 0
    # ~30x requested; realized is close. detection in (0, 1].
    assert float(row["mean_coverage"]) == pytest.approx(30.0, abs=2.0)
    det = float(row["detection"])
    assert 0.0 < det <= 1.0


def test_coverage_cli_single_end(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "cov-se-reads"),
        "-c", "10",
        "-r", "10",
        "--type", "single-end",
        "-s", "1",
    ]

    run_cli(cmd)

    se_path = tmp_path / "cov-se-reads.fastq.gz"
    assert se_path.exists(), f"Single-end file not found at {se_path}"


def test_coverage_tsv_cli(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cov_file = tmp_path / "coverages.tsv"
    cov_file.write_text(f"{str(tmp_path / 'input.fasta')}\t5\n")

    cmd = [
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(tmp_path / "cov-tsv-reads"),
        "-c", str(cov_file),
        "-r", "10",
        "-s", "1",
    ]

    run_cli(cmd)

    R1_path = tmp_path / "cov-tsv-reads_R1.fastq.gz"
    R2_path = tmp_path / "cov-tsv-reads_R2.fastq.gz"
    assert R1_path.exists(), f"R1 file not found at {R1_path}"
    assert R2_path.exists(), f"R2 file not found at {R2_path}"


def test_preflight_missing_input_fasta_exits(tmp_path):
    args = _base_args(input_fastas=[str(tmp_path / "nope.fasta")])
    with pytest.raises(SystemExit):
        preflight_checks(args)


def test_preflight_missing_coverage_file_exits(tmp_path):
    real = tmp_path / "in.fasta"
    real.write_text(">c\nACGT\n")
    args = _base_args(input_fastas=[str(real)], coverage=str(tmp_path / "missing.tsv"))
    with pytest.raises(SystemExit):
        preflight_checks(args)


def test_preflight_missing_proportions_file_exits(tmp_path):
    real = tmp_path / "in.fasta"
    real.write_text(">c\nACGT\n")
    args = _base_args(input_fastas=[str(real)], proportions_file=str(tmp_path / "missing.tsv"))
    with pytest.raises(SystemExit):
        preflight_checks(args)


def test_preflight_proportions_and_coverage_mutually_exclusive(tmp_path):
    real = tmp_path / "in.fasta"
    real.write_text(">c\nACGT\n")
    prop = tmp_path / "p.tsv"
    prop.write_text(f"{real}\t1\n")
    args = _base_args(input_fastas=[str(real)], proportions_file=str(prop), coverage="10")
    with pytest.raises(SystemExit):
        preflight_checks(args)


def test_get_proportions_missing_fasta_in_proportions_exits(tmp_path):
    f1 = tmp_path / "a.fasta"
    f1.write_text(">c\nACGT\n")
    prop = tmp_path / "p.tsv"
    prop.write_text("other.fasta\t1\n")
    args = _base_args(input_fastas=[str(f1)], proportions_file=str(prop), coverage=None)
    with pytest.raises(SystemExit):
        get_proportions(args)


def test_compute_coverage_missing_fasta_in_coverage_exits(tmp_path):
    f1 = tmp_path / "a.fasta"
    f1.write_text(">c\nACGT\n")
    cov = tmp_path / "c.tsv"
    cov.write_text("other.fasta\t10\n")
    args = _base_args(input_fastas=[str(f1)], coverage=str(cov), type="single-end")
    with pytest.raises(SystemExit):
        compute_reads_from_coverage(args)


def test_compute_coverage_zero_total_reads_exits(tmp_path):
    f1 = tmp_path / "a.fasta"
    f1.write_text(">c\nACGT\n")  # 4 bp
    args = _base_args(input_fastas=[str(f1)], coverage="0", type="single-end")
    with pytest.raises(SystemExit):
        compute_reads_from_coverage(args)


def test_per_read_tsv_paired_end_reconstructs_reads(tmp_path):
    fasta = tmp_path / "ref.fasta"
    ref = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"  # 40 bp
    fasta.write_text(f">refcontig\n{ref}\n")
    out = tmp_path / "pe"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=7, num_reads=20, read_length=8, fragment_size=20, fragment_size_range=10,
        output_prefix=str(out), circularize=False, include_Ns=False, per_read_tsv=True,
        type="paired-end")
    gen_reads(args, get_proportions(args))

    tsv = out.parent / (out.name + "-read-sources.tsv.gz")
    assert tsv.exists()

    # headers must be bare ids (no comment field) when per_read_tsv is on
    with open(str(out) + "_R1.fastq") as f:
        first_header = f.readline().rstrip("\n")
    assert " " not in first_header
    assert first_header.endswith("/1")

    r1_seqs = _fastq_map(str(out) + "_R1.fastq")
    r2_seqs = _fastq_map(str(out) + "_R2.fastq")
    refseq = _read_fasta_seq(fasta)

    with gzip.open(tsv, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        assert header == ["read_id", "source_fasta", "contig", "start", "end", "strand", "wrapped"]
        rows = 0
        for line in f:
            rid, _sf, _contig, start, end, strand, wrapped = line.rstrip("\n").split("\t")
            start, end = int(start), int(end)
            emitted = r1_seqs.get(rid) or r2_seqs.get(rid)
            assert emitted is not None, f"no fastq read for {rid}"
            # each non-wrapped read must reconstruct from its reported reference span
            if wrapped == "false":
                expected = refseq[start:end]
                if strand == "-":
                    expected = _revcomp(expected)
                assert emitted == expected, f"{rid}: {emitted} != {expected} (strand {strand})"
            rows += 1
    assert rows == 20  # 10 fragments x 2 mates


def test_per_read_tsv_single_end_reconstructs_reads(tmp_path):
    fasta = tmp_path / "ref.fasta"
    ref = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"
    fasta.write_text(f">refcontig\n{ref}\n")
    out = tmp_path / "se"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=3, num_reads=10, read_length=8, output_prefix=str(out),
        circularize=False, type="single-end", include_Ns=False, per_read_tsv=True)
    gen_reads(args, get_proportions(args))

    tsv = out.parent / (out.name + "-read-sources.tsv.gz")
    assert tsv.exists()
    refseq = _read_fasta_seq(fasta)
    seqs = _fastq_map(str(out) + ".fastq")

    with gzip.open(tsv, "rt") as f:
        f.readline()
        n = 0
        for line in f:
            rid, _sf, _contig, start, end, strand, wrapped = line.rstrip("\n").split("\t")
            if wrapped == "false":
                expected = refseq[int(start):int(end)]
                if strand == "-":
                    expected = _revcomp(expected)
                assert seqs[rid] == expected
            n += 1
    assert n == 10


def test_paired_end_subunity_proportion_yields_full_count(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">c0\n" + "ACGT" * 50 + "\n")  # 200 bp
    out = tmp_path / "pe_tail"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=0, num_reads=20, read_length=10, fragment_size=50, fragment_size_range=10,
        output_prefix=str(out), circularize=False, include_Ns=False, per_read_tsv=False,
        type="paired-end")
    # a single file given a sub-unity proportion is normalized to 1.0 and still
    # receives the full budget (20 reads -> 10 fragments -> 10 R1 + 10 R2)
    gen_reads(args, {str(fasta): 0.5})
    assert _count_reads(str(out) + "_R1.fastq") == 10
    assert _count_reads(str(out) + "_R2.fastq") == 10


def test_single_end_subunity_proportion_yields_full_count(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">c0\n" + "ACGT" * 50 + "\n")
    out = tmp_path / "se_tail"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=0, num_reads=20, read_length=10, output_prefix=str(out),
        circularize=False, type="single-end", include_Ns=False, per_read_tsv=False)
    gen_reads(args, {str(fasta): 0.5})
    assert _count_reads(str(out) + ".fastq") == 20


def test_long_read_subunity_proportion_yields_full_count(tmp_path):
    # exercises long-read length sampling with a sub-unity proportion
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">c0\n" + "ACGT" * 500 + "\n")  # 2000 bp
    out = tmp_path / "se_long_tail"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=0, num_reads=20, read_length=200, output_prefix=str(out),
        circularize=False, type="long", long_read_length_range=50,
        include_Ns=False, per_read_tsv=False)
    gen_reads(args, {str(fasta): 0.5})
    assert _count_reads(str(out) + ".fastq") == 20


def test_subunity_proportion_writes_source_rows_paired_end(tmp_path):
    # sub-unity proportion + per_read_tsv together, so provenance rows are still written
    # for the full normalized budget
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">c0\n" + "ACGT" * 50 + "\n")
    out = tmp_path / "pe_tail_src"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=0, num_reads=20, read_length=10, fragment_size=50, fragment_size_range=10,
        output_prefix=str(out), circularize=False, include_Ns=False, per_read_tsv=True,
        type="paired-end")
    gen_reads(args, {str(fasta): 0.5})
    tsv = out.parent / (out.name + "-read-sources.tsv.gz")
    assert tsv.exists()
    with gzip.open(tsv, "rt") as f:
        rows = [l for l in f.read().splitlines() if l]
    assert len(rows) == 21  # header + 20 reads


def test_subunity_proportion_writes_source_rows_single_end(tmp_path):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">c0\n" + "ACGT" * 50 + "\n")
    out = tmp_path / "se_tail_src"
    args = Namespace(input_fastas=[str(fasta)], proportions_file=None, coverage=None,
        seed=0, num_reads=20, read_length=10, output_prefix=str(out),
        circularize=False, type="single-end", include_Ns=False, per_read_tsv=True)
    gen_reads(args, {str(fasta): 0.5})
    tsv = out.parent / (out.name + "-read-sources.tsv.gz")
    assert tsv.exists()
    with gzip.open(tsv, "rt") as f:
        rows = [l for l in f.read().splitlines() if l]
    assert len(rows) == 21

# ---------------------------------------------------------------------------
# parallel read generation across multiple input genomes
#
# the parallel path generates each input fasta in its own process and
# concatenates the parts in input order. with a fixed --seed (per-file seeds),
# its output must be byte-for-byte identical to the single-process path, so most
# of these tests assert serial == parallel rather than re-deriving sequences.
# ---------------------------------------------------------------------------

def _write_multi_genomes(tmp_path):
    """ three small genomes of differing sizes/contig counts; returns list of paths """
    import random
    rng = random.Random(2024)
    specs = [("g1", [4000, 1500]), ("g2", [9000]), ("g3", [2500, 2500, 2500])]
    paths = []
    for name, lens in specs:
        p = tmp_path / f"{name}.fasta"
        with open(p, "w") as f:
            for i, L in enumerate(lens):
                f.write(f">{name}_c{i}\n" + "".join(rng.choice("ACGT") for _ in range(L)) + "\n")
        paths.append(str(p))
    return paths


def _multi_args(paths, out_prefix, jobs, rtype="paired-end", seed=123, per_read_tsv=True):
    return Namespace(
        input_fastas=list(paths),
        proportions_file=None,
        coverage=None,
        seed=seed,
        num_reads=4000,
        read_length=100,
        fragment_size=400,
        fragment_size_range=10,
        long_read_length_range=50,
        output_prefix=str(out_prefix),
        circularize=False,
        include_Ns=False,
        per_read_tsv=per_read_tsv,
        type=rtype,
        jobs=jobs,
    )


def _read_text(path):
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return f.read()


@pytest.mark.parametrize("rtype,suffixes", [
    ("paired-end", ["_R1.fastq", "_R2.fastq"]),
    ("single-end", [".fastq"]),
    ("long", [".fastq"]),
])
def test_parallel_matches_serial(tmp_path, rtype, suffixes):
    paths = _write_multi_genomes(tmp_path)

    s_args = _multi_args(paths, tmp_path / "serial", jobs=1, rtype=rtype)
    gen_reads(s_args, get_proportions(s_args))

    p_args = _multi_args(paths, tmp_path / "par", jobs=4, rtype=rtype)
    gen_reads(p_args, get_proportions(p_args))

    for suf in suffixes:
        assert _read_text(str(tmp_path / "serial") + suf) == _read_text(str(tmp_path / "par") + suf), \
            f"serial vs parallel mismatch in {suf}"

    # provenance TSV must match too
    assert _read_text(str(tmp_path / "serial") + "-read-sources.tsv.gz") == \
           _read_text(str(tmp_path / "par") + "-read-sources.tsv.gz")


def test_parallel_reproducible_across_job_counts(tmp_path):
    paths = _write_multi_genomes(tmp_path)

    a_args = _multi_args(paths, tmp_path / "j4", jobs=4, seed=77)
    gen_reads(a_args, get_proportions(a_args))
    b_args = _multi_args(paths, tmp_path / "j2", jobs=2, seed=77)
    gen_reads(b_args, get_proportions(b_args))

    assert _read_text(str(tmp_path / "j4") + "_R1.fastq") == _read_text(str(tmp_path / "j2") + "_R1.fastq")
    assert _read_text(str(tmp_path / "j4") + "_R2.fastq") == _read_text(str(tmp_path / "j2") + "_R2.fastq")


def test_single_file_output_independent_of_jobs(tmp_path):
    # a single input fasta always runs serially regardless of --jobs; output identical
    fasta = tmp_path / "solo.fasta"
    import random
    rng = random.Random(5)
    fasta.write_text(">c0\n" + "".join(rng.choice("ACGT") for _ in range(8000)) + "\n")

    a = _multi_args([str(fasta)], tmp_path / "j1", jobs=1, seed=9)
    gen_reads(a, get_proportions(a))
    b = _multi_args([str(fasta)], tmp_path / "j8", jobs=8, seed=9)
    gen_reads(b, get_proportions(b))

    assert _read_text(str(tmp_path / "j1") + "_R1.fastq") == _read_text(str(tmp_path / "j8") + "_R1.fastq")


def test_parallel_total_read_count_and_contiguous_ids(tmp_path):
    paths = _write_multi_genomes(tmp_path)
    p_args = _multi_args(paths, tmp_path / "cnt", jobs=4, rtype="paired-end")
    gen_reads(p_args, get_proportions(p_args))

    r1 = _fastq_map(str(tmp_path / "cnt") + "_R1.fastq")
    r2 = _fastq_map(str(tmp_path / "cnt") + "_R2.fastq")
    # num_reads=4000 -> 2000 fragments -> 2000 R1 + 2000 R2
    assert len(r1) == 2000
    assert len(r2) == 2000

    # base read ids run 1..2000 with no gaps, each present as /1 and /2
    base1 = sorted(int(rid[1:].split("/")[0]) for rid in r1)
    base2 = sorted(int(rid[1:].split("/")[0]) for rid in r2)
    assert base1 == list(range(1, 2001))
    assert base2 == list(range(1, 2001))


def test_parallel_per_read_tsv_maps_to_valid_genomes(tmp_path):
    paths = _write_multi_genomes(tmp_path)
    p_args = _multi_args(paths, tmp_path / "prov", jobs=3, rtype="single-end")
    gen_reads(p_args, get_proportions(p_args))

    valid = set(paths)
    with gzip.open(str(tmp_path / "prov") + "-read-sources.tsv.gz", "rt") as f:
        header = f.readline()
        assert header.startswith("read_id\t")
        rows = 0
        for line in f:
            _rid, source_fasta, _contig, _s, _e, _strand, _w = line.rstrip("\n").split("\t")
            assert source_fasta in valid
            rows += 1
    assert rows == p_args.num_reads


# ---------------------------------------------------------------------------
# apportionment / offset helpers underpinning the parallel split
# ---------------------------------------------------------------------------

def test_apportion_units_sums_exactly_and_is_ordered():
    files = ["a", "b", "c", "d"]
    props = {f: 0.25 for f in files}
    out = apportion_units(files, 999, props)
    assert sum(out.values()) == 999
    assert list(out.keys()) == files  # preserves input order
    # 999 over 4 equal -> three 250s and one 249, deterministic by input order
    assert out == {"a": 250, "b": 250, "c": 250, "d": 249}


def test_apportion_units_normalizes_unnormalized_proportions():
    # weights that don't sum to 1 must still distribute the full budget
    files = ["a", "b"]
    out = apportion_units(files, 100, {"a": 3.0, "b": 1.0})
    assert sum(out.values()) == 100
    assert out == {"a": 75, "b": 25}


def test_apportion_units_zero_and_single_file():
    assert apportion_units(["only"], 4242, {"only": 1.0}) == {"only": 4242}
    assert apportion_units(["a", "b"], 0, {"a": 0.5, "b": 0.5}) == {"a": 0, "b": 0}


def test_apportion_units_nonpositive_weights_fall_back_to_equal():
    out = apportion_units(["a", "b", "c"], 9, {"a": 0.0, "b": 0.0, "c": 0.0})
    assert sum(out.values()) == 9
    # equal fallback -> 3 each
    assert out == {"a": 3, "b": 3, "c": 3}


def test_compute_id_offsets_advances_by_units():
    files = ["a", "b", "c"]
    per_file = {"a": 250, "b": 0, "c": 100}  # b has zero budget
    offs = compute_id_offsets(files, per_file)
    # offsets advance by units; zero-budget file does not advance the counter
    assert offs == {"a": 0, "b": 250, "c": 250}


def test_distribute_units_across_contigs_sums_to_budget():
    counts = distribute_units_across_contigs([4000, 1500, 500], 123)
    assert sum(counts) == 123
    assert len(counts) == 3
    # roughly proportional: the largest contig gets the most
    assert counts[0] >= counts[1] >= counts[2]


def test_distribute_units_handles_zero_budget_and_empty():
    assert distribute_units_across_contigs([1000, 2000], 0) == [0, 0]
    assert distribute_units_across_contigs([0, 0], 10) == [0, 0]


# ─── ground-truth-assembly (GTA) ──────────────────────────────────────────────

from bit.modules.gen_reads_detection import DetectionTracker
from bit.modules.gen_reads import write_gta_records
import io


class _Rec:
    """ minimal stand-in for a parsed fasta record (id + seq). """
    def __init__(self, rid, seq):
        self.id = rid
        self.seq = seq


def _gta_headers(text):
    return [l for l in text.splitlines() if l.startswith(">")]


def _gta_seq_for(text, header_suffix):
    """ pull the sequence under the header whose line ends with header_suffix. """
    lines = text.splitlines()
    idx = next(i for i, l in enumerate(lines)
               if l.startswith(">") and l.endswith(header_suffix))
    seq = ""
    j = idx + 1
    while j < len(lines) and not lines[j].startswith(">"):
        seq += lines[j]
        j += 1
    return seq


def test_gta_linear_splits_on_uncovered_gaps():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 10, 30)
    tr.add(0, 60, 75)
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=False)
    assert n == 2
    assert _gta_headers(buf.getvalue()) == [
        ">g__c1:10-30", ">g__c1:60-75"]
    # sequence is sliced exactly from the reference
    assert _gta_seq_for(buf.getvalue(), "c1:10-30") == seq[10:30]


def test_gta_uncovered_contig_emits_nothing():
    seq = "ACGT" * 25
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])          # no reads added
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=False)
    assert n == 0
    assert buf.getvalue() == ""


def test_gta_provenance_header_uses_source_basename():
    seq = "ACGT" * 25
    rec = _Rec("contigX", seq)
    tr = DetectionTracker([100])
    tr.add(0, 0, 100)
    buf = io.StringIO()
    write_gta_records(buf, "/some/dir/mygenome.fasta", [rec], tr)
    assert _gta_headers(buf.getvalue()) == [">mygenome__contigX:0-100"]


def test_gta_strips_fasta_suffixes_in_header():
    from bit.modules.gen_reads import _strip_fasta_suffix
    # every recognized suffix -> stripped to the bare accession
    for suf in (".fasta.gz", ".fasta", ".fna.gz", ".fna", ".fa.gz", ".fa"):
        assert _strip_fasta_suffix("GCA_1.2" + suf) == "GCA_1.2"
    # case-insensitive
    assert _strip_fasta_suffix("GCA_1.2.FASTA.GZ") == "GCA_1.2"
    # only the trailing suffix goes; internal dots (incl. accession version) stay
    assert _strip_fasta_suffix("weird.name.fasta.gz") == "weird.name"
    assert _strip_fasta_suffix("GCA_1.2") == "GCA_1.2"      # no suffix -> unchanged
    # a bare .gz that is not a fasta suffix is left alone
    assert _strip_fasta_suffix("something.gz") == "something.gz"


def test_gta_header_from_gzipped_source(tmp_path):
    # end-to-end style: a compressed source basename should appear suffix-free
    seq = "ACGT" * 25
    rec = _Rec("CAZNLJ010000001.1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 0, 100)
    buf = io.StringIO()
    write_gta_records(buf, "GCA_964573725.1.fasta.gz", [rec], tr)
    assert _gta_headers(buf.getvalue()) == [
        ">GCA_964573725.1__CAZNLJ010000001.1:0-100"]


def test_gta_circular_stitches_origin_spanning_run():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 80, 100)     # end
    tr.add(0, 0, 30)       # start  -> wraps with the end under circular
    tr.add(0, 50, 60)      # middle island
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=True)
    assert n == 2
    headers = _gta_headers(buf.getvalue())
    # origin-spanning contig has end < start (80-30); middle stays linear
    assert ">g__c1:80-30" in headers
    assert ">g__c1:50-60" in headers
    # stitched sequence is seq[80:100] + seq[0:30]
    assert _gta_seq_for(buf.getvalue(), "c1:80-30") == seq[80:100] + seq[0:30]


def test_gta_circular_same_data_stays_split_when_linear():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 80, 100)
    tr.add(0, 0, 30)
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=False)
    assert n == 2
    assert _gta_headers(buf.getvalue()) == [
        ">g__c1:0-30", ">g__c1:80-100"]


def test_gta_circular_full_contig_stays_single_piece():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 0, 100)
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=True)
    assert n == 1
    assert _gta_headers(buf.getvalue()) == [">g__c1:0-100"]


def test_gta_circular_start_only_does_not_wrap():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    rec = _Rec("c1", seq)
    tr = DetectionTracker([100])
    tr.add(0, 0, 30)       # touches start but not the end
    tr.add(0, 50, 70)
    buf = io.StringIO()
    n = write_gta_records(buf, "g.fasta", [rec], tr, circularize=True)
    assert n == 2
    assert _gta_headers(buf.getvalue()) == [
        ">g__c1:0-30", ">g__c1:50-70"]


def test_gta_multi_genome_keeps_origins_distinct():
    # two genomes with distinct contig ids -> combined GTA tags each by source
    seqA = "".join("ACGT"[i % 4] for i in range(60))
    seqB = "".join("TGCA"[i % 4] for i in range(60))
    trA = DetectionTracker([60]); trA.add(0, 0, 60)
    trB = DetectionTracker([60]); trB.add(0, 0, 60)
    buf = io.StringIO()
    write_gta_records(buf, "A.fasta", [_Rec("cA", seqA)], trA)
    write_gta_records(buf, "B.fasta", [_Rec("cB", seqB)], trB)
    headers = _gta_headers(buf.getvalue())
    assert ">A__cA:0-60" in headers
    assert ">B__cB:0-60" in headers


def test_gta_cli_paired_end_writes_file(tmp_path):
    shutil.copy(test_fasta, tmp_path / "input.fasta")
    out = tmp_path / "gta-reads"
    run_cli([
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(out),
        "-c", "20", "-r", "150",
        "--type", "paired-end",
        "--ground-truth-assembly",
        "-j", "1", "-s", "9",
    ])
    gta = tmp_path / "gta-reads-ground-truth-assembly.fasta"
    assert gta.exists(), f"GTA file not found at {gta}"
    text = gta.read_text()
    headers = _gta_headers(text)
    assert headers, "GTA has no contigs"
    # every GTA piece is a substring of the source contig it names
    src_seq = _read_fasta_seq(str(tmp_path / "input.fasta"))
    for h in headers:
        piece = _gta_seq_for(text, h.split("__", 1)[1])
        assert piece.upper() in src_seq, f"GTA piece {h} not found in source"


def test_gta_cli_off_by_default(tmp_path):
    shutil.copy(test_fasta, tmp_path / "input.fasta")
    out = tmp_path / "noreq"
    run_cli([
        "bit", "gen-reads",
        "-i", str(tmp_path / "input.fasta"),
        "-o", str(out),
        "-c", "5", "-r", "150",
        "--type", "paired-end",
        "-j", "1", "-s", "9",
    ])
    assert not (tmp_path / "noreq-ground-truth-assembly.fasta").exists()
