import os
from argparse import Namespace
import shutil
import gzip
import pytest # type: ignore
from bit.modules.gen_reads import parse_proportions_file
from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
from bit.modules.gen_reads import (gen_paired_reads,
                                   gen_single_reads,
                                   get_proportions,
                                   extract_subsequence,
                                   parse_coverages_input)


test_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")

def test_gen_reads(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit-gen-reads",
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
        expected = f.read().splitlines()

    observed = [
        "@partial-NC_003131.1_1_1529/1",
        "TTTTTGTTTT",
        "+",
        "IIIIIIIIII"
    ]

    assert expected == observed, f"R1 content does not match expected:\n{observed}\nExpected:\n{expected}"


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
    )

    gen_paired_reads(args, get_proportions(args))

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

    # ensure at least one wrapped forward read appears
    assert any(s in wrapped_fw for s in seqs), "No wrapped fragments found with circularize=True"


def test_gen_reads_single_end(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit-gen-reads",
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
    )

    gen_single_reads(args, get_proportions(args))

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
    )

    gen_single_reads(args, get_proportions(args))

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
        "bit-gen-reads",
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
    )

    gen_paired_reads(args, get_proportions(args))

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
        "bit-gen-reads",
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


def test_coverage_cli_single_end(tmp_path):

    shutil.copy(test_fasta, tmp_path / "input.fasta")

    cmd = [
        "bit-gen-reads",
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
        "bit-gen-reads",
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