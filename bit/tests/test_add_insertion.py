import pytest
from argparse import Namespace
from bit.modules.add_insertion import (add_insertion, parse_fasta, parse_insertion_fasta,
                                       parse_position_string, validate_position)
from bit.tests.utils import run_cli


def _write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def _read_fasta(path):
    return parse_fasta(str(path))



class TestParsePositionString:
    def test_valid(self):
        assert parse_position_string("contig_1:1000") == ("contig_1", 1000)

    def test_zero_position(self):
        assert parse_position_string("chr1:0") == ("chr1", 0)

    def test_missing_colon(self):
        with pytest.raises(SystemExit):
            parse_position_string("contig1_1000")

    def test_non_integer_position(self):
        with pytest.raises(SystemExit):
            parse_position_string("contig_1:abc")

    def test_negative_position(self):
        with pytest.raises(SystemExit):
            parse_position_string("contig_1:-5")


class TestParseFasta:
    def test_empty_file(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        with pytest.raises(SystemExit):
            parse_fasta(str(empty))

    def test_single_record(self, tmp_path):
        fa = tmp_path / "one.fasta"
        _write_fasta(fa, [("seq1", "ATCG")])
        result = parse_fasta(str(fa))
        assert list(result.keys()) == ["seq1"]
        assert result["seq1"] == "ATCG"


class TestParseInsertionFasta:
    def test_multiple_records(self, tmp_path):
        fa = tmp_path / "multi.fasta"
        _write_fasta(fa, [("a", "ATCG"), ("b", "GGCC")])
        with pytest.raises(SystemExit):
            parse_insertion_fasta(str(fa))


class TestValidatePosition:
    def test_bad_contig(self, tmp_path):
        fa = tmp_path / "input.fasta"
        _write_fasta(fa, [("contig_1", "ATCGATCG")])
        contigs = parse_fasta(str(fa))
        with pytest.raises(SystemExit):
            validate_position("contig_99:0", contigs)

    def test_position_out_of_range(self, tmp_path):
        fa = tmp_path / "input.fasta"
        _write_fasta(fa, [("contig_1", "ATCG")])
        contigs = parse_fasta(str(fa))
        with pytest.raises(SystemExit):
            validate_position("contig_1:100", contigs)



class TestAddInsertionSingle:
    """single insertion at an explicit position"""

    def test_single_insertion_at_position(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("contig_1", "AAAACCCC")])
        _write_fasta(ins, [("ins", "TTTT")])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="contig_1:4",
            num_insertions=1,
            back_to_back=False,
            seed=42,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert contigs["contig_1"] == "AAAATTTTCCCC"
        assert result["num_insertions"] == 1
        assert result["positions"] == [("contig_1", 4)]
        assert result["insertion_length"] == 4


    def test_position_with_multiple_first_specified_rest_random(self, tmp_path):
        """when -p and -n > 1, first goes at specified position, rest are random"""
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        original_seq = "A" * 200
        insert_seq = "TT"
        n = 3
        _write_fasta(fa, [("contig_1", original_seq)])
        _write_fasta(ins, [("ins", insert_seq)])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="contig_1:100",
            num_insertions=n,
            back_to_back=False,
            seed=42,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert len(contigs["contig_1"]) == len(original_seq) + n * len(insert_seq)
        assert len(result["positions"]) == n
        # first insertion should be at the specified position
        assert result["positions"][0] == ("contig_1", 100)
        # remaining positions should not all be at the same spot (with this seed/length)
        remaining = result["positions"][1:]
        assert any(pos != ("contig_1", 100) for pos in remaining)


class TestAddInsertionRandom:
    """single insertion at a random position"""

    def test_random_insertion_increases_length(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        original_seq = "A" * 100
        insert_seq = "TTTT"
        _write_fasta(fa, [("contig_1", original_seq)])
        _write_fasta(ins, [("ins", insert_seq)])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position=None,
            num_insertions=1,
            back_to_back=False,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert len(contigs["contig_1"]) == len(original_seq) + len(insert_seq)
        assert insert_seq in contigs["contig_1"]


class TestAddInsertionMultipleRandom:
    """multiple insertions spaced randomly"""

    def test_multiple_random_insertions(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        original_seq = "A" * 200
        insert_seq = "GG"
        n = 3

        _write_fasta(fa, [("contig_1", original_seq)])
        _write_fasta(ins, [("ins", insert_seq)])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position=None,
            num_insertions=n,
            back_to_back=False,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert len(contigs["contig_1"]) == len(original_seq) + n * len(insert_seq)
        assert len(result["positions"]) == n


class TestAddInsertionBackToBack:
    """multiple insertions back-to-back at a single position"""

    def test_back_to_back_at_explicit_position(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("contig_1", "AAAACCCC")])
        _write_fasta(ins, [("ins", "TT")])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="contig_1:4",
            num_insertions=3,
            back_to_back=True,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        # 3 copies of "TT" inserted at position 4
        assert contigs["contig_1"] == "AAAA" + "TT" * 3 + "CCCC"
        assert len(result["positions"]) == 3

    def test_back_to_back_random_position(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        original_seq = "A" * 50
        insert_seq = "CC"
        n = 4

        _write_fasta(fa, [("contig_1", original_seq)])
        _write_fasta(ins, [("ins", insert_seq)])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position=None,
            num_insertions=n,
            back_to_back=True,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert len(contigs["contig_1"]) == len(original_seq) + n * len(insert_seq)
        # all positions should be the same since they're back-to-back
        assert len(set(result["positions"])) == 1


class TestAddInsertionMultiContig:
    """input fasta with multiple contigs"""

    def test_insertion_into_multi_contig(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("ctg1", "AAAA"), ("ctg2", "CCCC"), ("ctg3", "GGGG")])
        _write_fasta(ins, [("ins", "TT")])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="ctg2:2",
            num_insertions=1,
            back_to_back=False,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)

        assert contigs["ctg1"] == "AAAA"   # unchanged
        assert contigs["ctg2"] == "CCTTCC" # insertion at position 2
        assert contigs["ctg3"] == "GGGG"   # unchanged


class TestAddInsertionAtBoundaries:
    """insertion at position 0 (start) and at end"""

    def test_insertion_at_start(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("seq1", "ATCG")])
        _write_fasta(ins, [("ins", "XX")])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="seq1:0",
            num_insertions=1,
            back_to_back=False,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)
        assert contigs["seq1"] == "XXATCG"

    def test_insertion_at_end(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("seq1", "ATCG")])
        _write_fasta(ins, [("ins", "XX")])

        args = Namespace(
            input_fasta=str(fa),
            insertion_fasta=str(ins),
            output_fasta=str(out),
            position="seq1:4",
            num_insertions=1,
            back_to_back=False,
            seed=9,
        )

        result = add_insertion(args)
        contigs = _read_fasta(out)
        assert contigs["seq1"] == "ATCGXX"


class TestCLI:
    def test_cli_basic_run(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"

        _write_fasta(fa, [("contig_1", "AAAACCCC")])
        _write_fasta(ins, [("ins", "TTTT")])

        cmd = [
            "bit-add-insertion",
            "-i", str(fa),
            "-I", str(ins),
            "-o", str(out),
            "-p", "contig_1:4",
            "-s", "9",
        ]

        run_cli(cmd)

        assert out.exists()
        contigs = _read_fasta(out)
        assert contigs["contig_1"] == "AAAATTTTCCCC"

    def test_cli_help_exits_zero(self):
        import subprocess
        result = subprocess.run(
            ["bit-add-insertion", "-h"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "insertion" in result.stdout.lower()

    def test_cli_log_file(self, tmp_path):
        fa = tmp_path / "input.fasta"
        ins = tmp_path / "insert.fasta"
        out = tmp_path / "output.fasta"
        log = tmp_path / "insertion.log"

        _write_fasta(fa, [("contig_1", "AAAACCCC")])
        _write_fasta(ins, [("ins", "TTTT")])

        cmd = [
            "bit-add-insertion",
            "-i", str(fa),
            "-I", str(ins),
            "-o", str(out),
            "-p", "contig_1:4",
            "-n", "2",
            "--back-to-back",
            "-s", "9",
            "-l", str(log),
        ]

        run_cli(cmd)

        assert log.exists()
        contents = log.read_text()
        assert "contig_1\t4" in contents

        lines = [l for l in contents.strip().splitlines()]
        assert lines[0] == "insertion_number\tcontig\tposition"
        assert len(lines) == 3
