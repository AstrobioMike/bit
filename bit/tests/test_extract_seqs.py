import pytest # type: ignore
from types import SimpleNamespace

from bit.modules.extract_seqs import (
    extract_seqs_by_coords,
    extract_seqs_by_headers,
    extract_seqs_by_primers,
    find_all_primer_hits,
    find_amplicons,
    find_invalid_primer_chars,
    count_alt_orientation_amplicons,
    should_report_alt_orientation,
)
from bit.modules.seqs import revcomp


# shared synthetic amplicon pieces for primer tests
PRIMER_F = "ACGTAGCTAG"
PRIMER_R = "TGACTGGATC"                                      # revcomp = GATCCAGTCA
AMPLICON = PRIMER_F + "TTTTTTTTTT" + revcomp(PRIMER_R)       # 30 bp, plus-strand orientation
SPACER = "GGGGGGGGGG"


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
            inverse=False,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert set(result.keys()) == {"seq1", "seq3"}
        assert result["seq1"] == "ATGCATGCATGC"
        assert result["seq3"] == "AAATTTGGGCCC"


    def test_extract_matching_headers_with_leading_gt(self, tmp_path, simple_fasta):
        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=[">seq1", "seq3"],
            inverse=False,
        )
        extract_seqs_by_headers(args)

        result = read_fasta_dict(out)
        assert set(result.keys()) == {"seq1", "seq3"}


    def test_extract_matching_headers_from_file(self, tmp_path, simple_fasta):
        header_file = tmp_path / "headers.txt"
        header_file.write_text("seq2\n")
        out = str(tmp_path / "out.fasta")

        args = SimpleNamespace(
            input_fasta=simple_fasta,
            output_fasta=out,
            headers=[str(header_file)],
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
            bed=[str(bed)],
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
            bed=[str(bed)],
            output_fasta=str(out),
        )
        extract_seqs_by_coords(args)

        result = read_fasta_dict(str(out))
        assert len(result) == 2


    def test_inline_region(self, tmp_path):
        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [("chr1", "ATGCATGCATGCATGC")])

        out = tmp_path / "out.fasta"
        args = SimpleNamespace(
            input_fasta=str(fasta),
            bed=["chr1", "2", "6"],
            output_fasta=str(out),
        )
        extract_seqs_by_coords(args)

        result = read_fasta_dict(str(out))
        assert len(result) == 1
        assert list(result.values())[0].upper() == "GCAT"


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
        rev = "GGGGGGGG"  # no site in seq, either orientation

        hits_strict = find_all_primer_hits(seq, fwd, rev, max_mismatches=0)
        hits_relaxed = find_all_primer_hits(seq, fwd, rev, max_mismatches=1)

        assert len(hits_relaxed) == 1
        assert len(hits_strict) == 0


    def test_fwd_and_rev_at_identical_coords_both_kept(self):
        # fwd and revcomp(rev) land on the same coordinates; both must survive dedup
        seq = "TTTTTACGTAGCTAGGATCCTTTTT"
        fwd = "ACGTAGCTAGGA"
        rev = "TCCTAGCTACGT"  # revcomp of fwd

        hits = find_all_primer_hits(seq, fwd, rev, max_mismatches=0)

        assert hits == [
            ("fwd", 5, 17, "ACGTAGCTAGGA"),
            ("rev_rc", 5, 17, "ACGTAGCTAGGA"),
        ]


    def test_worse_sites_found_alongside_exact_site(self):
        # edlib only reports locations at the best edit distance found, so an exact site used
        # to hide 1-difference sites elsewhere; all sites within max_mismatches must be returned
        seq = "CCCC" + AMPLICON + SPACER + AMPLICON.replace(PRIMER_F, "ACGTAGCTAC") + "CCCC"

        hits_strict = find_all_primer_hits(seq, PRIMER_F, "GGGGGGGGGGGG", max_mismatches=0)
        hits_relaxed = find_all_primer_hits(seq, PRIMER_F, "GGGGGGGGGGGG", max_mismatches=1)

        assert {h[1] for h in hits_strict} == {4}
        assert {h[1] for h in hits_relaxed} == {4, 44}


    @pytest.mark.parametrize("base, should_hit", [("A", True), ("G", True), ("C", False), ("T", False)])
    def test_degenerate_primer_base(self, base, should_hit):
        # R = A or G
        seq = "TTAC" + base + "TAGCTAGTT"

        hits = find_all_primer_hits(seq, "ACRTAGCTAG", "CCCCCCCCCC", max_mismatches=0)

        if should_hit:
            assert hits == [("fwd", 2, 12, "AC" + base + "TAGCTAG")]
        else:
            assert hits == []


    def test_degenerate_primer_hits_via_revcomp(self):
        # the reverse primer's revcomp must turn its ambiguity codes into their complements
        # (Y = C/T -> R = A/G), or the revcomp search would look for the wrong bases
        seq = "CCCC" + AMPLICON + "CCCC"
        rev = "TGACTGGAYC"  # revcomp = GRTCCAGTCA, matches GATCCAGTCA

        hits = find_all_primer_hits(seq, PRIMER_F, rev, max_mismatches=0)

        assert ("rev_rc", 24, 34, "GATCCAGTCA") in hits


    def test_target_ambiguity_codes_do_not_match(self):
        # IUPAC equalities are symmetric in edlib; Ns in the target must not match primer bases
        seq = "CCCC" + "N" * 20 + "CCCC"

        hits = find_all_primer_hits(seq, PRIMER_F, PRIMER_R, max_mismatches=0)

        assert hits == []


    def test_primer_must_be_longer_than_allowed_differences(self):
        with pytest.raises(ValueError):
            find_all_primer_hits("ACGTACGT", "AC", PRIMER_R, max_mismatches=2)


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


    def test_fully_overlapping_primers(self):
        seq = "TTTTTACGTAGCTAGGATCCTTTTT"
        fwd = "ACGTAGCTAGGA"
        rev = "TCCTAGCTACGT"  # revcomp of fwd

        amplicons = find_amplicons(seq, fwd, rev, max_mismatches=0)

        assert amplicons == [("fwd", "rev_rc", 5, 17, "ACGTAGCTAGGA", 12)]


    def test_same_seq_given_as_both_primers(self):
        # a primer found as-given in both slots points the same direction twice, so nothing amplifies
        seq = "TTTTTACGTAGCTAGGATCCTTTTT"
        primer = "ACGTAGCTAGGA"

        amplicons = find_amplicons(seq, primer, primer, max_mismatches=0)

        assert amplicons == []


    def test_partially_overlapping_primers(self):
        # 20 bp target, 12 bp primers overlapping by 4
        seq = "TTTTACGTAGCTAGGATCCAGTCATTTT"
        fwd = "ACGTAGCTAGGA"
        rev = "TGACTGGATCCT"  # revcomp of target[8:]

        amplicons = find_amplicons(seq, fwd, rev, max_mismatches=0)

        assert amplicons == [("fwd", "rev_rc", 4, 24, "ACGTAGCTAGGATCCAGTCA", 20)]


    def test_rev_hit_contained_within_fwd_hit_not_reported(self):
        seq = "TTTTACGTAGCTAGGATCCTTTT"
        fwd = "ACGTAGCTAGGATCC"
        rev = "CTAGCT"  # revcomp lands inside the fwd hit

        amplicons = find_amplicons(seq, fwd, rev, max_mismatches=0)

        assert amplicons == []


    def test_palindromic_primer_reported_once(self):
        # GAATTC hits as both fwd and fwd_rc at the same coords; region reported once
        seq = "TTTTGAATTCCCCCCCCCGGGTTTAAAA"
        fwd = "GAATTC"
        rev = "AAACCC"

        amplicons = find_amplicons(seq, fwd, rev, max_mismatches=0)

        assert amplicons == [("fwd", "rev_rc", 4, 24, "GAATTCCCCCCCCCGGGTTT", 20)]


    def test_multi_copy_target_one_amplicon_per_copy(self):
        # e.g. multiple rRNA operons, one on the minus strand; each copy should give one product,
        # not every combination of forward and reverse sites across copies
        seq = "CCCCC" + AMPLICON + SPACER + AMPLICON + SPACER + revcomp(AMPLICON) + "CCCCC"

        amplicons = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0)

        assert amplicons == [
            ("fwd", "rev_rc", 5, 35, AMPLICON, 30),
            ("fwd", "rev_rc", 45, 75, AMPLICON, 30),
            ("rev", "fwd_rc", 85, 115, revcomp(AMPLICON), 30),
        ]


    def test_minus_strand_amplicon(self):
        seq = "CCCC" + revcomp(AMPLICON) + "CCCC"

        amplicons = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0)

        assert amplicons == [("rev", "fwd_rc", 4, 34, revcomp(AMPLICON), 30)]


    def test_invalid_orientation_not_paired(self):
        # both primers found as given means they point the same direction
        seq = "CCCC" + PRIMER_F + "TTTTTTTTTT" + PRIMER_R + "CCCC"

        amplicons = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0)

        assert amplicons == []


    def test_innermost_left_site_kept(self):
        # two forward sites upstream of one reverse site; only the shorter product is kept
        seq = "CCCC" + PRIMER_F + "GG" + AMPLICON + "CCCC"

        amplicons = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0)

        assert amplicons == [("fwd", "rev_rc", 16, 46, AMPLICON, 30)]


    def test_degenerate_primers_amplicon(self):
        seq = "CCCC" + AMPLICON + "CCCC"

        assert find_amplicons(seq, "ACGTRGCTAG", "TGACTGGATN", max_mismatches=0) == [
            ("fwd", "rev_rc", 4, 34, AMPLICON, 30)
        ]
        # Y = C/T doesn't match the A at that position
        assert find_amplicons(seq, "ACGTYGCTAG", PRIMER_R, max_mismatches=0) == []


    def test_mismatch_copy_found_alongside_exact_copy(self):
        seq = "CCCC" + AMPLICON + SPACER + AMPLICON.replace(PRIMER_F, "ACGTAGCTAC") + "CCCC"

        strict = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0)
        relaxed = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=1)

        assert [(a[2], a[3]) for a in strict] == [(4, 34)]
        assert [(a[2], a[3]) for a in relaxed] == [(4, 34), (44, 74)]


    @pytest.mark.parametrize("max_len, expected_count", [(None, 3), (30, 3), (29, 0)])
    def test_max_amplicon_length(self, max_len, expected_count):
        seq = "CCCCC" + AMPLICON + SPACER + AMPLICON + SPACER + revcomp(AMPLICON) + "CCCCC"

        amplicons = find_amplicons(seq, PRIMER_F, PRIMER_R, max_mismatches=0,
                                   max_amplicon_length=max_len)

        assert len(amplicons) == expected_count


class TestAltOrientation:

    # 3 copies, the last on the minus strand, as rRNA operons can be
    MIXED_STRAND_SEQ = "CCCCC" + AMPLICON + SPACER + AMPLICON + SPACER + revcomp(AMPLICON) + "CCCCC"

    @pytest.mark.parametrize("fwd, rev, expected", [
        (PRIMER_F, PRIMER_R, {"reverse": 2, "forward": 0, "both": 1}),
        (PRIMER_F, revcomp(PRIMER_R), {"reverse": 3, "forward": 1, "both": 0}),
        (revcomp(PRIMER_F), PRIMER_R, {"reverse": 1, "forward": 3, "both": 2}),
        (revcomp(PRIMER_F), revcomp(PRIMER_R), {"reverse": 0, "forward": 2, "both": 3}),
        ("AAAAAAAAAA", "CCCCCAAAAA", {"reverse": 0, "forward": 0, "both": 0}),
    ])
    def test_counts(self, fwd, rev, expected):
        hits = find_all_primer_hits(self.MIXED_STRAND_SEQ, fwd, rev, max_mismatches=0)

        assert count_alt_orientation_amplicons(self.MIXED_STRAND_SEQ, hits) == expected


    def test_correct_primers_always_beat_alternatives(self):
        # with correctly oriented primers, flipped orientations can only produce cross-copy spans
        seq = self.MIXED_STRAND_SEQ
        hits = find_all_primer_hits(seq, PRIMER_F, PRIMER_R, max_mismatches=0)
        n_given = len(find_amplicons(seq, PRIMER_F, PRIMER_R))

        assert not should_report_alt_orientation(count_alt_orientation_amplicons(seq, hits), n_given)


    def test_max_amplicon_length_applies_to_alternatives(self):
        seq = "CCCC" + AMPLICON + "CCCC"
        hits = find_all_primer_hits(seq, PRIMER_F, revcomp(PRIMER_R), max_mismatches=0)

        assert count_alt_orientation_amplicons(seq, hits)["reverse"] == 1
        assert count_alt_orientation_amplicons(seq, hits, max_amplicon_length=29)["reverse"] == 0


    @pytest.mark.parametrize("counts, hit_count, expected", [
        ({"reverse": 3, "forward": 0, "both": 0}, 0, True),
        ({"reverse": 3, "forward": 0, "both": 0}, 2, True),
        ({"reverse": 2, "forward": 0, "both": 1}, 3, False),
        ({"reverse": 3, "forward": 0, "both": 0}, 3, False),
        ({"reverse": 0, "forward": 0, "both": 0}, 0, False),
    ])
    def test_should_report(self, counts, hit_count, expected):
        assert should_report_alt_orientation(counts, hit_count) is expected


class TestFindInvalidPrimerChars:

    def test_iupac_codes_valid(self):
        assert find_invalid_primer_chars("ACGTRYSWKMBDHVN") == []
        assert find_invalid_primer_chars("acgtryswkmbdhvn") == []


    def test_invalid_chars_reported(self):
        assert find_invalid_primer_chars("ACGX-U") == ["-", "U", "X"]


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
            max_amplicon_length=None,
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
            max_amplicon_length=None,
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
            max_amplicon_length=None,
        )
        extract_seqs_by_primers(args)

        result = read_fasta_dict(out)
        returned_headers = {h.split("|")[0] for h in result}
        assert "contig2" not in returned_headers


    def test_max_amplicon_length_filters_output(self, tmp_path):
        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [
            ("short", "CCCC" + AMPLICON + "CCCC"),
            ("long", "CCCC" + PRIMER_F + "T" * 100 + revcomp(PRIMER_R) + "CCCC"),
        ])

        out = str(tmp_path / "out.fasta")
        args = SimpleNamespace(
            input_fasta=str(fasta),
            output_fasta=out,
            forward_primer=PRIMER_F,
            reverse_primer=PRIMER_R,
            max_mismatches=0,
            max_amplicon_length=50,
        )
        extract_seqs_by_primers(args)

        result = read_fasta_dict(out)
        assert {h.split("|")[0] for h in result} == {"short"}


    @pytest.mark.parametrize("fwd, rev, expected_suggestion", [
        (PRIMER_F, PRIMER_R, None),
        (PRIMER_F, revcomp(PRIMER_R), f"-r {PRIMER_R}"),
        (revcomp(PRIMER_F), PRIMER_R, f"-f {PRIMER_F}"),
        (revcomp(PRIMER_F), revcomp(PRIMER_R), f"-f {PRIMER_F} -r {PRIMER_R}"),
        ("AAAAAAAAAA", "CCCCCAAAAA", None),
    ])
    def test_orientation_notice(self, tmp_path, capsys, fwd, rev, expected_suggestion):
        fasta = tmp_path / "input.fasta"
        write_fasta(fasta, [
            ("mixed", "CCCCC" + AMPLICON + SPACER + AMPLICON + SPACER + revcomp(AMPLICON) + "CCCCC"),
        ])

        args = SimpleNamespace(
            input_fasta=str(fasta),
            output_fasta=str(tmp_path / "out.fasta"),
            forward_primer=fwd,
            reverse_primer=rev,
            max_mismatches=0,
            max_amplicon_length=None,
        )
        extract_seqs_by_primers(args)

        out = capsys.readouterr().out
        if expected_suggestion is None:
            assert "NOTICE" not in out
        else:
            assert "NOTICE" in out
            assert expected_suggestion in out


class TestExtractByPrimersCli:

    def _parse(self, *extra):
        from bit.cli.fasta import build_parser
        return build_parser().parse_args(["extract-by-primers", "-i", "x.fa", "-f", PRIMER_F,
                                          "-r", PRIMER_R, *extra])


    def test_default_max_amplicon_length(self):
        assert self._parse().max_amplicon_length == 10000


    @pytest.mark.parametrize("value, expected", [("0", None), ("500", 500)])
    def test_max_amplicon_length_zero_means_no_limit(self, monkeypatch, value, expected):
        import bit.cli.fasta as cli
        import bit.modules.extract_seqs as mod
        seen = {}
        monkeypatch.setattr(mod, "extract_seqs_by_primers", lambda args: seen.update(L=args.max_amplicon_length))
        cli.run_extract_by_primers(self._parse("-L", value))
        assert seen["L"] == expected


    def test_negative_max_amplicon_length_rejected(self):
        import bit.cli.fasta as cli
        with pytest.raises(SystemExit):
            cli.run_extract_by_primers(self._parse("-L", "-1"))
