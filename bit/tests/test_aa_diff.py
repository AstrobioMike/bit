from unittest import mock
from Bio.Seq import Seq  # type: ignore

from bit.modules.aa_diff import (
    run_aa_diff,
    _parse_cs_tag,
    _align,
    _get_gapped_seqs,
    _parse_alignment,
    _collect_mutations,
    _write_tsv,
    _write_mutations,
    _write_summary,
    _write_alignment,
    _report_summary,
)


class TestParseCsTag:

    def test_match_run(self):
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":3", ref, 0)
        assert translated == "MKL"
        assert frameshifts == []

    def test_match_run_with_offset(self):
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":2", ref, 2)
        assert translated == "LR"
        assert frameshifts == []

    def test_substitution(self):
        # *acgR: codon "acg" -> "ACG" -> Thr (T); ref AA was R
        ref = "MR"
        translated, frameshifts = _parse_cs_tag("*acgR", ref, 1)
        assert translated == str(Seq("ACG").translate())
        assert frameshifts == []

    def test_insertion(self):
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":2+KL:1", ref, 0)
        assert translated == "MKS"
        assert frameshifts == []

    def test_deletion_uppercase(self):
        # -AA means 2 AAs deleted from ref (uppercase = protein level)
        ref = "MKAA"
        translated, frameshifts = _parse_cs_tag(":2-AA", ref, 0)
        assert translated == "MK"
        assert frameshifts == []

    def test_frameshift_lowercase_1nt(self):
        # -a: 1 extra nt = +1 frameshift
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":3-a:1", ref, 0)
        assert translated == "MKLR"
        assert len(frameshifts) == 1
        assert frameshifts[0]["type"] == "+1"
        assert frameshifts[0]["ref_pos"] == 4  # 0-based ref_pos=3 -> 1-based 4

    def test_frameshift_lowercase_2nt(self):
        # -ac: 2 extra nts = +2 frameshift
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":3-ac", ref, 0)
        assert len(frameshifts) == 1
        assert frameshifts[0]["type"] == "+2"

    def test_frameshift_f_tag(self):
        # fa: alternate +1 frameshift encoding
        ref = "MKL"
        translated, frameshifts = _parse_cs_tag(":2fa", ref, 0)
        assert len(frameshifts) == 1
        assert frameshifts[0]["type"] == "+1"
        assert frameshifts[0]["ref_pos"] == 3

    def test_frameshift_b_tag(self):
        # bac: alternate +2 frameshift encoding
        ref = "MKL"
        translated, frameshifts = _parse_cs_tag(":1bac", ref, 0)
        assert len(frameshifts) == 1
        assert frameshifts[0]["type"] == "+2"

    def test_combined_operations(self):
        ref = "MKLRST"
        # :2 matches MK, substitution at pos 2 (L -> T via acg), :1 matches R
        translated, frameshifts = _parse_cs_tag(":2*acgL:1", ref, 0)
        thr = str(Seq("ACG").translate())
        assert translated == "MK" + thr + "R"
        assert frameshifts == []

    def test_intron_is_ignored(self):
        ref = "MKLRST"
        translated, frameshifts = _parse_cs_tag(":2~+100gt:2", ref, 0)
        assert translated == "MKLR"
        assert frameshifts == []

    def test_empty_string(self):
        translated, frameshifts = _parse_cs_tag("", "MKLRST", 0)
        assert translated == ""
        assert frameshifts == []


class TestAlignAndGetGappedSeqs:

    def test_identical_sequences(self):
        seq = "MKLRST"
        alignment = _align(seq, seq)
        ref_gapped, qry_gapped = _get_gapped_seqs(alignment, seq, seq)
        assert ref_gapped == seq
        assert qry_gapped == seq

    def test_deletion_in_query(self):
        ref = "MKLRST"
        qry = "MKRST"   # L is deleted
        alignment = _align(ref, qry)
        ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref, qry)
        assert len(ref_gapped) == len(qry_gapped)
        # the gapped query must contain a '-' where the ref has L
        assert "-" in qry_gapped

    def test_insertion_in_query(self):
        ref = "MKRST"
        qry = "MKLRST"   # L is inserted
        alignment = _align(ref, qry)
        ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref, qry)
        assert len(ref_gapped) == len(qry_gapped)
        assert "-" in ref_gapped

    def test_gapped_seqs_same_length(self):
        ref = "MAAAKRST"
        qry = "MKLRST"
        alignment = _align(ref, qry)
        ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref, qry)
        assert len(ref_gapped) == len(qry_gapped)


class TestParseAlignment:

    def test_all_matches(self):
        positions, insertions = _parse_alignment("MKLRS", "MKLRS")
        assert len(positions) == 5
        assert all(p["change_type"] == "match" for p in positions)
        assert insertions == []

    def test_positions_numbered_from_1(self):
        positions, _ = _parse_alignment("MK", "MK")
        assert positions[0]["ref_pos"] == 1
        assert positions[1]["ref_pos"] == 2

    def test_substitution(self):
        positions, insertions = _parse_alignment("MKL", "MRL")
        assert positions[1]["change_type"] == "substitution"
        assert positions[1]["ref_aa"] == "K"
        assert positions[1]["query_aa"] == "R"
        assert insertions == []

    def test_deletion_in_query(self):
        # ref has 3 residues, query has a '-' at position 2
        positions, insertions = _parse_alignment("MKL", "M-L")
        assert positions[1]["change_type"] == "deletion"
        assert positions[1]["query_aa"] == "-"
        assert positions[1]["query_pos"] is None

    def test_insertion_in_query(self):
        # ref has '-' at col 1 (insertion in query)
        positions, insertions = _parse_alignment("M-L", "MKL")
        assert len(insertions) == 1
        assert insertions[0]["after_ref_pos"] == 1
        assert insertions[0]["inserted_seq"] == "K"
        # only 2 ref residues (M and L)
        assert len(positions) == 2

    def test_insertion_at_end(self):
        positions, insertions = _parse_alignment("MK--", "MKXY")
        assert len(positions) == 2
        assert len(insertions) == 1
        assert insertions[0]["after_ref_pos"] == 2
        assert insertions[0]["inserted_seq"] == "XY"

    def test_query_position_tracking(self):
        positions, _ = _parse_alignment("MKL", "MRL")
        assert positions[0]["query_pos"] == 1
        assert positions[1]["query_pos"] == 2
        assert positions[2]["query_pos"] == 3

    def test_query_position_skips_deletion(self):
        # query position should not increment for a deletion
        positions, _ = _parse_alignment("MKL", "M-L")
        assert positions[0]["query_pos"] == 1
        assert positions[1]["query_pos"] is None
        assert positions[2]["query_pos"] == 2


class TestCollectMutations:

    def test_no_mutations(self):
        positions = [
            {"ref_pos": 1, "ref_aa": "M", "query_aa": "M", "change_type": "match"},
            {"ref_pos": 2, "ref_aa": "K", "query_aa": "K", "change_type": "match"},
        ]
        assert _collect_mutations(positions, []) == []

    def test_substitution(self):
        positions = [
            {"ref_pos": 1, "ref_aa": "M", "query_aa": "M", "change_type": "match"},
            {"ref_pos": 2, "ref_aa": "K", "query_aa": "R", "change_type": "substitution"},
        ]
        mutations = _collect_mutations(positions, [])
        assert mutations == ["K2R"]

    def test_deletion(self):
        positions = [
            {"ref_pos": 1, "ref_aa": "M", "query_aa": "-", "change_type": "deletion"},
        ]
        mutations = _collect_mutations(positions, [])
        assert mutations == ["M1del"]

    def test_insertion(self):
        mutations = _collect_mutations([], [{"after_ref_pos": 5, "inserted_seq": "KL"}])
        assert mutations == ["ins5:KL"]

    def test_frameshift(self):
        mutations = _collect_mutations([], [], [{"ref_pos": 10, "type": "+1"}])
        assert mutations == ["fs10+1"]

    def test_frameshift_none(self):
        # frameshifts=None should be treated the same as empty list
        mutations = _collect_mutations([], [], None)
        assert mutations == []

    def test_combined(self):
        positions = [
            {"ref_pos": 1, "ref_aa": "M", "query_aa": "V", "change_type": "substitution"},
            {"ref_pos": 2, "ref_aa": "K", "query_aa": "-", "change_type": "deletion"},
        ]
        insertions = [{"after_ref_pos": 2, "inserted_seq": "XY"}]
        frameshifts = [{"ref_pos": 5, "type": "+2"}]
        mutations = _collect_mutations(positions, insertions, frameshifts)
        assert "M1V" in mutations
        assert "K2del" in mutations
        assert "ins2:XY" in mutations
        assert "fs5+2" in mutations


class TestWriteTsv:

    def test_header_row(self, tmp_path):
        path = tmp_path / "out.tsv"
        _write_tsv([], [], str(path))
        header = path.read_text().splitlines()[0]
        assert header == "ref_pos\tref_aa\tquery_aa\tquery_pos\tchange_type\tinserted_before"

    def test_match_row(self, tmp_path):
        positions = [{"ref_pos": 1, "ref_aa": "M", "query_aa": "M", "query_pos": 1, "change_type": "match"}]
        path = tmp_path / "out.tsv"
        _write_tsv(positions, [], str(path))
        rows = path.read_text().splitlines()
        assert rows[1] == "1\tM\tM\t1\tmatch\t-"

    def test_deletion_row_uses_dash_for_query_pos(self, tmp_path):
        positions = [{"ref_pos": 1, "ref_aa": "M", "query_aa": "-", "query_pos": None, "change_type": "deletion"}]
        path = tmp_path / "out.tsv"
        _write_tsv(positions, [], str(path))
        rows = path.read_text().splitlines()
        assert rows[1].split("\t")[3] == "-"

    def test_insertion_recorded_in_inserted_before_column(self, tmp_path):
        positions = [
            {"ref_pos": 1, "ref_aa": "M", "query_aa": "M", "query_pos": 1, "change_type": "match"},
            {"ref_pos": 2, "ref_aa": "K", "query_aa": "K", "query_pos": 2, "change_type": "match"},
        ]
        insertions = [{"after_ref_pos": 1, "inserted_seq": "XY"}]
        path = tmp_path / "out.tsv"
        _write_tsv(positions, insertions, str(path))
        rows = path.read_text().splitlines()
        # row for ref_pos 2 should show the insertion that came after ref_pos 1
        assert rows[2].split("\t")[5] == "XY"

    def test_row_count_matches_positions(self, tmp_path):
        positions = [
            {"ref_pos": i + 1, "ref_aa": "M", "query_aa": "M", "query_pos": i + 1, "change_type": "match"}
            for i in range(5)
        ]
        path = tmp_path / "out.tsv"
        _write_tsv(positions, [], str(path))
        lines = path.read_text().splitlines()
        assert len(lines) == 6  # 1 header + 5 rows


class TestWriteMutations:

    def test_writes_each_mutation(self, tmp_path):
        path = tmp_path / "mutations.txt"
        _write_mutations(["K2R", "M1del"], str(path))
        content = path.read_text()
        assert "K2R\n" in content
        assert "M1del\n" in content

    def test_empty_mutations_writes_comment(self, tmp_path):
        path = tmp_path / "mutations.txt"
        _write_mutations([], str(path))
        assert "# No mutations detected" in path.read_text()


class TestWriteSummary:

    def test_content_written(self, tmp_path):
        path = tmp_path / "summary.txt"
        _write_summary("some stats text", str(path))
        assert path.read_text() == "some stats text\n"


class TestWriteAlignment:

    def test_creates_file(self, tmp_path):
        path = tmp_path / "aln.txt"
        _write_alignment("MKL", "MKL", "ref_id", "query_id", str(path))
        assert path.exists()

    def test_contains_seq_labels(self, tmp_path):
        path = tmp_path / "aln.txt"
        _write_alignment("MKLRST", "MKLRST", "ref_id", "query_id", str(path))
        content = path.read_text()
        assert "ref_id" in content
        assert "query_id" in content

    def test_contains_sequence_data(self, tmp_path):
        path = tmp_path / "aln.txt"
        _write_alignment("MKLRST", "MKXRST", "ref1", "qry1", str(path))
        content = path.read_text()
        assert "MKLRST" in content
        assert "MKXRST" in content

    def test_match_line_pipes_for_identical(self, tmp_path):
        path = tmp_path / "aln.txt"
        _write_alignment("MKLRST", "MKLRST", "ref1", "qry1", str(path))
        content = path.read_text()
        # all positions match, so the match line should be all pipes
        assert "||||||" in content

    def test_width_parameter_splits_output(self, tmp_path):
        # with width=3 and a 6-residue sequence there should be 2 blocks
        path = tmp_path / "aln.txt"
        _write_alignment("MKLRST", "MKLRST", "ref1", "qry1", str(path), width=3)
        # each block writes 5 lines plus a blank line = 6 lines per block, 2 blocks
        non_empty = [l for l in path.read_text().splitlines() if l.strip()]
        # at least two separate sequence lines each containing "MKL" and "RST"
        seq_lines = [l for l in path.read_text().splitlines() if "MKL" in l or "RST" in l]
        assert len(seq_lines) >= 2


class TestReportSummary:

    def _make_positions(self, n_match=3, n_sub=1, n_del=1):
        positions = []
        ref_pos = 0
        qry_pos = 0
        for _ in range(n_match):
            ref_pos += 1
            qry_pos += 1
            positions.append({"ref_pos": ref_pos, "ref_aa": "M", "query_aa": "M",
                               "query_pos": qry_pos, "change_type": "match"})
        for _ in range(n_sub):
            ref_pos += 1
            qry_pos += 1
            positions.append({"ref_pos": ref_pos, "ref_aa": "K", "query_aa": "R",
                               "query_pos": qry_pos, "change_type": "substitution"})
        for _ in range(n_del):
            ref_pos += 1
            positions.append({"ref_pos": ref_pos, "ref_aa": "L", "query_aa": "-",
                               "query_pos": None, "change_type": "deletion"})
        return positions

    def test_aa_mode_returns_text(self, capsys):
        positions = self._make_positions()
        mutations = ["K4R", "L5del"]
        with mock.patch("bit.modules.general.color_text", side_effect=lambda t, _c: t):
            text = _report_summary(positions, [], mutations, [], "test-prefix")
        assert isinstance(text, str)
        assert len(text) > 0

    def test_aa_mode_contains_ref_length(self, capsys):
        positions = self._make_positions(n_match=3, n_sub=0, n_del=0)
        with mock.patch("bit.modules.general.color_text", side_effect=lambda t, _c: t):
            text = _report_summary(positions, [], [], [], "test-prefix")
        assert "Reference length" in text
        assert "3" in text

    def test_nt_mode_shows_nt_query_len(self, capsys):
        positions = self._make_positions(n_match=5, n_sub=0, n_del=0)
        with mock.patch("bit.modules.general.color_text", side_effect=lambda t, _c: t):
            text = _report_summary(positions, [], [], [], "test-prefix",
                                   translated_path="test-inferred-protein.faa",
                                   nt_query_len=1500)
        assert "Query nt length" in text
        assert "1,500" in text

    def test_frameshifts_reported_in_nt_mode(self, capsys):
        positions = self._make_positions(n_match=5, n_sub=0, n_del=0)
        frameshifts = [{"ref_pos": 3, "type": "+1"}]
        with mock.patch("bit.modules.general.color_text", side_effect=lambda t, _c: t):
            text = _report_summary(positions, [], [], frameshifts, "test-prefix",
                                   translated_path="test-inferred-protein.faa",
                                   nt_query_len=1500)
        assert "Frameshifts" in text

    def test_substitution_count(self, capsys):
        positions = self._make_positions(n_match=3, n_sub=2, n_del=0)
        mutations = ["K4R", "K5R"]
        with mock.patch("bit.modules.general.color_text", side_effect=lambda t, _c: t):
            text = _report_summary(positions, [], mutations, [], "test-prefix")
        assert "Substitutions" in text
        assert "2" in text


# 30-AA test protein and a query with two substitutions (L3V, G9S).
# The NT sequence encodes the same query protein (30 codons + stop = 93 nt).
_REF_PROTEIN = "MKLRSTAEGVDNPQHCWFYIAKLTPFWGRD"
_QUERY_AA    = "MKVRSTAESVDNPQHCWFYIAKLTPFWGRD"   # L3V, G9S vs ref
_QUERY_NT    = (                                    # encodes _QUERY_AA
    "ATGAAAGTTCGTTCTACCGCTGAAAGTGTTGATAATCCTCAACATTGTTGGTTTTATATT"
    "GCTAAACTGACCCCTTTTTGGGGTCGTGATTAA"
)

class TestRunAaDiffAaMode:

    def _write_fastas(self, tmp_path, ref_seq=_REF_PROTEIN, qry_seq=_QUERY_AA,
                     ref_id="ref", qry_id="query"):
        ref = tmp_path / "ref.faa"
        ref.write_text(f">{ref_id}\n{ref_seq}\n")
        qry = tmp_path / "query.faa"
        qry.write_text(f">{qry_id}\n{qry_seq}\n")
        return str(ref), str(qry)

    def test_all_output_files_created(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        assert (tmp_path / "out-all-positions.tsv").exists()
        assert (tmp_path / "out-mutations.txt").exists()
        assert (tmp_path / "out-summary.txt").exists()
        assert (tmp_path / "out-alignment.txt").exists()

    def test_expected_substitutions_detected(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        mutations = (tmp_path / "out-mutations.txt").read_text().splitlines()
        assert "L3V" in mutations
        assert "G9S" in mutations

    def test_no_mutations_for_identical_seqs(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path, qry_seq=_REF_PROTEIN)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        assert "# No mutations detected" in (tmp_path / "out-mutations.txt").read_text()

    def test_tsv_row_count_equals_ref_length(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        lines = (tmp_path / "out-all-positions.tsv").read_text().splitlines()
        assert len(lines) == len(_REF_PROTEIN) + 1  # header + one row per ref residue

    def test_summary_reports_ref_length(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        summary = (tmp_path / "out-summary.txt").read_text()
        assert "Reference length" in summary
        assert str(len(_REF_PROTEIN)) in summary

    def test_only_substitutions_no_indels(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "aa", prefix)
        mutations = (tmp_path / "out-mutations.txt").read_text().splitlines()
        assert not any("del" in m or m.startswith("ins") for m in mutations)


class TestRunAaDiffNtMode:

    def _write_fastas(self, tmp_path):
        ref = tmp_path / "ref.faa"
        ref.write_text(f">ref\n{_REF_PROTEIN}\n")
        qry = tmp_path / "query.fna"
        qry.write_text(f">query_nt\n{_QUERY_NT}\n")
        return str(ref), str(qry)

    def test_all_output_files_created(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "nt", prefix)
        assert (tmp_path / "out-all-positions.tsv").exists()
        assert (tmp_path / "out-mutations.txt").exists()
        assert (tmp_path / "out-summary.txt").exists()
        assert (tmp_path / "out-alignment.txt").exists()
        assert (tmp_path / "out-inferred-protein.faa").exists()
        assert (tmp_path / "out-inferred-cds.fasta").exists()

    def test_expected_substitutions_detected(self, tmp_path):
        ref, qry = self._write_fastas(tmp_path)
        prefix = str(tmp_path / "out")
        run_aa_diff(qry, ref, "nt", prefix)
        mutations = (tmp_path / "out-mutations.txt").read_text().splitlines()
        assert "L3V" in mutations
        assert "G9S" in mutations
