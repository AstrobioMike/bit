import pytest # type: ignore
from types import SimpleNamespace

from bit.modules.extract_seqs import (
    extract_seqs_by_coords,
    extract_seqs_by_headers,
    extract_seqs_by_primers,
    find_all_primer_hits,
    find_amplicons,
)


# helpers
def write_fasta(path, records):
    """Write a list of (header, seq) tuples to a FASTA file."""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def read_fasta_dict(path):
    """Read a FASTA file and return {header: seq}."""
    seqs = {}
    with open(path) as f:
        header = None
        parts = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
        if header:
            seqs[header] = "".join(parts)
    return seqs


@pytest.fixture
def simple_fasta(tmp_path):

    fasta = tmp_path / "input.fasta"
    write_fasta(fasta, [
        ("seq1", "ATGCATGCATGC"),
        ("seq2", "GGGGCCCCTTTTAAAA"),
        ("seq3", "AAATTTGGGCCC"),
    ])

    return str(fasta)


class TestExtractSeqsByHeaders:

    def test_extract_matching_headers_from_list(self, tmp_path, simple_fasta):
        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=["seq1", "seq3"],
            file_with_headers=None,
            inverse=False,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert set(result.keys()) == {"seq1", "seq3"}
        assert result["seq1"] == "ATGCATGCATGC"
        assert result["seq3"] == "AAATTTGGGCCC"


    def test_extract_matching_headers_from_file(self, tmp_path, simple_fasta):
        header_file = tmp_path / "headers.txt"
        header_file.write_text("seq2\n")
        out = str(tmp_path / "out.fasta")

        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=None,
            file_with_headers=str(header_file),
            inverse=False,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert list(result.keys()) == ["seq2"]
        assert result["seq2"] == "GGGGCCCCTTTTAAAA"


    def test_inverse_flag(self, tmp_path, simple_fasta):
        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=["seq2"],
            file_with_headers=None,
            inverse=True,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert set(result.keys()) == {"seq1", "seq3"}


    def test_partial_match(self, tmp_path, simple_fasta):
        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=["seq1", "missing"],
            file_with_headers=None,
            inverse=False,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert list(result.keys()) == ["seq1"]


class TestExtractSeqsByCoords:

    def test_basic_coordinate_extraction(self, tmp_path):
        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [("chr1", "ATGCATGCATGCATGC")])

        bed = tmp_path / "coords.bed"
        bed.write_text("chr1\t2\t6\n")

        out = tmp_path / "out.fasta"

        args = SimpleNamespace(
            input_fasta=str(fasta),
            bed_file=str(bed),
            output_fasta=str(out),
        )
        extract_seqs_by_coords(args)

        assert out.exists()
        result = read_fasta_dict(str(out))
        assert len(result) == 1
        seq = list(result.values())[0]
        assert seq.upper() == "GCAT"


    def test_multiple_regions(self, tmp_path):
        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [("chr1", "ATGCATGCATGCATGC")])

        bed = tmp_path / "coords.bed"
        bed.write_text("chr1\t0\t4\nchr1\t8\t12\n")

        out = tmp_path / "out.fasta"
        args = SimpleNamespace(
            input_fasta=str(fasta),
            bed_file=str(bed),
            output_fasta=str(out),
        )
        extract_seqs_by_coords(args)

        result = read_fasta_dict(str(out))
        assert len(result) == 2


class TestFindAllPrimerHits:

    def test_exact_match_forward_and_reverse(self):
        #        0         1         2         3
        #        0123456789012345678901234567890123456789
        seq  = "AAAAAATGCATGCCCCCCCCCCCGCATGCATTTTTTT"
        fwd  = "ATGCATGC"
        rev  = "GCATGCAT"

        hits = find_all_primer_hits(seq, fwd, rev, max_mismatches=0)

        assert len(hits) > 0

        labels = {h[0] for h in hits}
        assert "fwd" in labels or "fwd_rc" in labels or "rev" in labels or "rev_rc" in labels


    def test_no_hits(self):
        seq = "AAAAAAAAAA"
        fwd = "CCGG"
        rev = "TTCC"

        hits = find_all_primer_hits(seq, fwd, rev, max_mismatches=0)
        assert hits == []


    def test_mismatches_expand_hits(self):
        seq = "AAAAATGCATGCAAAA"
        fwd = "ATGAATGC"  # 1 mismatch vs ATGCATGC
        rev = "NNNNNNNN"

        hits_strict = find_all_primer_hits(seq, fwd, rev, max_mismatches=0)
        hits_relaxed = find_all_primer_hits(seq, fwd, rev, max_mismatches=1)

        assert len(hits_relaxed) == 1
        assert len(hits_strict) == 0


class TestFindAmplicons:

    def test_simple_amplicon(self):

        fwd = "AAAA"
        rev = "CCCC"  # reverse complement = GGGG
        seq = "AAAATTTTTTTTGGGG"

        amplicons = find_amplicons(seq, fwd, rev, max_mismatches=0)
        assert len(amplicons) >= 1

        # The amplicon should span fwd through revcomp(rev)
        for left_label, right_label, left_start, right_end, amplicon, length in amplicons:
            assert length == right_end - left_start
            assert amplicon == seq[left_start:right_end]


    def test_amplicon_with_mismatch(self):

        fwd = "AAAA"
        rev = "CCCC"
        # one mismatch in each primer site
        seq = "AAGATTTTTTTTGGGC"

        amplicons_strict = find_amplicons(seq, fwd, rev, max_mismatches=0)
        amplicons_relaxed = find_amplicons(seq, fwd, rev, max_mismatches=1)

        assert len(amplicons_relaxed) >= 1
        assert len(amplicons_strict) == 0


class TestExtractSeqsByPrimers:

    def test_basic_primer_extraction(self, tmp_path):
        fwd = "AAAA"
        rev = "CCCC"
        seq = "AAAATTTTTTTTGGGG"

        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [("contig1", seq)])

        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=str(fasta),
            output_fasta=out,
            forward_primer=fwd,
            reverse_primer=rev,
            max_mismatches=0,
        )
        extract_seqs_by_primers(args)

        result = read_fasta_dict(out)
        assert list(result.values())[0] == "AAAATTTTTTTTGGGG"


    def test_output_header_format(self, tmp_path):
        fwd = "AAAA"
        rev = "CCCC"
        seq = "AAAATTTTTTTTGGGG"

        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [("contig1", seq)])

        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=str(fasta),
            output_fasta=out,
            forward_primer=fwd,
            reverse_primer=rev,
            max_mismatches=0,
        )
        extract_seqs_by_primers(args)

        result = read_fasta_dict(out)
        for header in result:
            # headers follow the format:
            # {original_header}|{left_label}-to-{right_label}|{start}-{end}|{length}
            parts = header.split("|")
            assert len(parts) == 4, f"Expected 4 pipe-separated parts in header, got: {header}"
            assert parts[0] == "contig1"
            assert "-to-" in parts[1]


    def test_multi_sequence_input(self, tmp_path):
        fwd = "AAAA"
        rev = "CCCC"

        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [
            ("contig1", "AAAATTTTTTTTGGGG"),
            ("contig2", "TTTTTTTTTTTTTTTT"),  # no hit
            ("contig3", "AAAACCCCCCCCGGGG"),
        ])

        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=str(fasta),
            output_fasta=out,
            forward_primer=fwd,
            reverse_primer=rev,
            max_mismatches=0,
        )
        extract_seqs_by_primers(args)

        result = read_fasta_dict(out)
        returned_headers = {h.split("|")[0] for h in result}
        assert "contig2" not in returned_headers
