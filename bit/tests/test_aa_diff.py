import os
import shutil
import pytest # type: ignore
from Bio.Seq import Seq # type: ignore
from bit.modules import aa_diff
from bit.modules import general


# one fixed codon per residue, so we can build queries whose translation is known
_CODON = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA",
    "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA",
    "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG",
    "Y": "TAT", "V": "GTT", "*": "TAA",
}


def cds_of(protein):
    """deterministic coding sequence for a protein string"""
    return "".join(_CODON[a] for a in protein)


def make_hit(prot_start, prot_end, strand, nt_start, nt_end, nt_len=None, name="q"):
    """
    build a minimal miniprot PAF row for _extend_over_clips, which only reads
    indices 2,3,4,5,6,7,8. other columns are placeholders.
    """
    if nt_len is None:
        nt_len = nt_end
    row = ["ref", "0", str(prot_start), str(prot_end), strand,
           name, str(nt_len), str(nt_start), str(nt_end)]
    return row


@pytest.fixture
def fail_loud(monkeypatch):
    """
    redirect the module's report_message / notify_premature_exit so the
    'fail-loud' paths raise a catchable exception instead of exiting the process.
    returns (ExceptionType, captured_messages_list).
    """
    class PrematureExit(Exception):
        pass

    messages = []

    def fake_report(msg, **kwargs):
        messages.append(msg)

    def fake_exit(*args, **kwargs):
        raise PrematureExit()

    monkeypatch.setattr(general, "report_message", fake_report)
    monkeypatch.setattr(general, "notify_premature_exit", fake_exit)
    return PrematureExit, messages


class TestParseCsTag:

    REF = "ACDEFG"

    def test_all_matches(self):
        seq, fs, introns = aa_diff._parse_cs_tag(":6", self.REF, 0)
        assert seq == "ACDEFG"
        assert fs == []
        assert introns == []

    def test_substitution_translates_query_codon(self):
        # :2 (AC) then a sub at ref position 3 whose query codon GTT -> V, then :3 (EFG)
        seq, fs, introns = aa_diff._parse_cs_tag(":2*gttD:3", self.REF, 0)
        assert seq == "ACVEFG"          # the D is replaced by the translated query codon V
        assert fs == []
        assert introns == []

    def test_deletion_skips_ref_residues(self):
        # +DE removes ref residues D and E from the query
        seq, fs, introns = aa_diff._parse_cs_tag(":2+DE:2", self.REF, 0)
        assert seq == "ACFG"
        assert fs == []
        assert introns == []

    def test_insertion_whole_codon(self):
        # -gtt inserts one codon (V) that has no counterpart in the reference
        seq, fs, introns = aa_diff._parse_cs_tag(":3-gtt:3", self.REF, 0)
        assert seq == "ACDVEFG"
        assert fs == []
        assert introns == []

    def test_frameshift_plus_one_F_op(self):
        seq, fs, introns = aa_diff._parse_cs_tag(":3-g:3", self.REF, 0)
        assert seq == "ACDEFG"          # extra nt restores frame; no residue emitted
        assert fs == [{"ref_pos": 4, "type": "+1"}]
        assert introns == []

    def test_frameshift_plus_two_F_op(self):
        seq, fs, introns = aa_diff._parse_cs_tag(":3-gc:3", self.REF, 0)
        assert seq == "ACDEFG"
        assert fs == [{"ref_pos": 4, "type": "+2"}]

    def test_frameshift_minus_one_G_op(self):
        # *atE -> 2-nt codon (1-nt deletion); ref AA E kept as placeholder, fs -1
        seq, fs, introns = aa_diff._parse_cs_tag(":3*atE:2", self.REF, 0)
        assert seq == "ACDEFG"
        assert fs == [{"ref_pos": 4, "type": "-1"}]
        assert introns == []

    def test_frameshift_minus_two_G_op(self):
        # *aE -> 1-nt codon (2-nt deletion); fs -2
        seq, fs, introns = aa_diff._parse_cs_tag(":3*aE:2", self.REF, 0)
        assert seq == "ACDEFG"
        assert fs == [{"ref_pos": 4, "type": "-2"}]

    def test_phase0_intron_between_codons(self):
        seq, fs, introns = aa_diff._parse_cs_tag(":3~gt50ag:3", self.REF, 0)
        assert seq == "ACDEFG"
        assert fs == []
        assert introns == [(9, 50)]      # offset = 3 codons * 3 nt; length 50

    def test_phase2_split_codon(self):
        # codon split 2|1 across an intron: head 'gt' + tail 't' -> GTT -> V
        seq, fs, introns = aa_diff._parse_cs_tag(":3*gtE~gt50ag-t:2", self.REF, 0)
        assert seq == "ACDVFG"
        assert fs == []
        assert introns == [(11, 50)]     # offset = 9 (matches) + 2 (head)

    def test_phase1_split_codon(self):
        # codon split 1|2 across an intron: head 'g' + tail 'tt' -> GTT -> V
        seq, fs, introns = aa_diff._parse_cs_tag(":3*gE~gt50ag-tt:2", self.REF, 0)
        assert seq == "ACDVFG"
        assert fs == []
        assert introns == [(10, 50)]     # offset = 9 (matches) + 1 (head)

    def test_ref_start_offset_is_respected(self):
        # a hit that starts partway into the reference (prot_start > 0)
        seq, fs, introns = aa_diff._parse_cs_tag(":3", self.REF, 3)
        assert seq == "EFG"             # ref[3:6]
        assert introns == []

    def test_fail_loud_on_internal_garbage(self, fail_loud):
        PrematureExit, messages = fail_loud
        with pytest.raises(PrematureExit):
            aa_diff._parse_cs_tag(":3ZZZ:3", self.REF, 0)
        assert messages and "unrecognized token" in messages[0]

    def test_fail_loud_on_trailing_garbage(self, fail_loud):
        PrematureExit, _ = fail_loud
        with pytest.raises(PrematureExit):
            aa_diff._parse_cs_tag(":6XYZ", self.REF, 0)


# a real miniprot cs for TP53 genomic vs its protein: phase-0/1/2 introns all
# present, no frameshifts. The whole 393-aa protein must reconstruct exactly.
TP53 = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA"
    "PPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKT"
    "CPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRN"
    "TFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGR"
    "DRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL"
    "KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)
TP53_CS = (
    ":24*ctL~gt117ag-a:7~gt109ag:93~gt757ag:61*gG~gt81ag-gt:37~gt568ag:36"
    "*agS~gt343ag-t:45*gA~gt92ag-ca:24~gt2819ag:35*agS~gt918ag-c:26"
)


def test_parse_cs_tag_real_tp53():
    seq, fs, introns = aa_diff._parse_cs_tag(TP53_CS, TP53, 0)
    assert seq == TP53          # all three intron phases reconstruct the full protein
    assert fs == []
    assert len(introns) == 9


class TestParseAlignment:

    def test_all_match(self):
        positions, insertions = aa_diff._parse_alignment("ACDEF", "ACDEF")
        assert insertions == []
        assert len(positions) == 5
        assert all(p["change_type"] == "match" for p in positions)
        assert positions[0] == {"ref_pos": 1, "ref_aa": "A", "query_pos": 1,
                                "query_aa": "A", "change_type": "match"}

    def test_substitution(self):
        positions, _ = aa_diff._parse_alignment("ACDEF", "AGDEF")
        sub = positions[1]
        assert sub["change_type"] == "substitution"
        assert sub["ref_aa"] == "C" and sub["query_aa"] == "G" and sub["ref_pos"] == 2

    def test_deletion(self):
        positions, _ = aa_diff._parse_alignment("ACDEF", "A-DEF")
        deln = positions[1]
        assert deln["change_type"] == "deletion"
        assert deln["query_aa"] == "-"
        assert deln["query_pos"] is None

    def test_internal_insertion(self):
        positions, insertions = aa_diff._parse_alignment("ACD-EF", "ACDGEF")
        assert insertions == [{"after_ref_pos": 3, "inserted_seq": "G"}]
        assert len(positions) == 5

    def test_leading_insertion(self):
        _, insertions = aa_diff._parse_alignment("-ACDEF", "GACDEF")
        assert insertions == [{"after_ref_pos": 0, "inserted_seq": "G"}]

    def test_trailing_insertion(self):
        _, insertions = aa_diff._parse_alignment("ACDEF-", "ACDEFG")
        assert insertions == [{"after_ref_pos": 5, "inserted_seq": "G"}]


class TestGappedSeqs:

    def _gapped(self, ref, query):
        aln = aa_diff._align(ref, query)
        return aa_diff._get_gapped_seqs(aln, ref, query)

    def test_identical(self):
        r, q = self._gapped("ACDEFGHIK", "ACDEFGHIK")
        assert r == "ACDEFGHIK"
        assert q == "ACDEFGHIK"

    def test_terminal_deletion_keeps_trailing_ref(self):
        # query missing the last two residues -> ref's tail must still appear,
        # opposite gaps in the query (regression for the trailing-portion handling)
        r, q = self._gapped("ACDEFGHIK", "ACDEFGH")
        assert r == "ACDEFGHIK"
        assert q == "ACDEFGH--"

    def test_terminal_insertion_keeps_trailing_query(self):
        r, q = self._gapped("ACDEFG", "ACDEFGHI")
        assert r == "ACDEFG--"
        assert q == "ACDEFGHI"


def _positions(change_types):
    out = []
    for i, ct in enumerate(change_types, 1):
        out.append({"ref_pos": i, "ref_aa": "A",
                    "query_pos": None if ct == "deletion" else i,
                    "query_aa": "-" if ct == "deletion" else "A",
                    "change_type": ct})
    return out


class TestCalcStats:

    def test_all_match(self):
        s = aa_diff._calc_stats_from_positions(_positions(["match"] * 5))
        assert s["perc_id"] == 100.0
        assert s["perc_ref_cov"] == 100.0
        assert s["n_match"] == 5

    def test_one_substitution(self):
        s = aa_diff._calc_stats_from_positions(_positions(["match", "substitution", "match", "match", "match"]))
        assert s["perc_id"] == 80.0
        assert s["perc_ref_cov"] == 100.0

    def test_one_deletion_lowers_coverage(self):
        s = aa_diff._calc_stats_from_positions(_positions(["match", "deletion", "match", "match", "match"]))
        assert s["perc_id"] == 80.0
        assert s["perc_ref_cov"] == 80.0

    def test_empty_is_zero_not_error(self):
        s = aa_diff._calc_stats_from_positions([])
        assert s["perc_id"] == 0.0
        assert s["perc_ref_cov"] == 0.0


class TestThresholds:

    def test_passes_when_above_thresholds(self, fail_loud):
        # should NOT raise
        aa_diff._check_alignment_thresholds(_positions(["match"] * 10), 30, 25)

    def test_fails_on_low_identity(self, fail_loud):
        PrematureExit, _ = fail_loud
        positions = _positions(["substitution"] * 9 + ["match"])   # 10% id
        with pytest.raises(PrematureExit):
            aa_diff._check_alignment_thresholds(positions, 30, 25)

    def test_fails_on_empty(self, fail_loud):
        PrematureExit, _ = fail_loud
        with pytest.raises(PrematureExit):
            aa_diff._check_alignment_thresholds([], 30, 25)


class TestCollectMutations:

    def test_substitution_string(self):
        pos = [{"change_type": "substitution", "ref_aa": "C", "ref_pos": 2, "query_aa": "G"}]
        assert aa_diff._collect_mutations(pos, []) == ["C2G"]

    def test_deletion_string(self):
        pos = [{"change_type": "deletion", "ref_aa": "D", "ref_pos": 4, "query_aa": "-"}]
        assert aa_diff._collect_mutations(pos, []) == ["D4del"]

    def test_insertion_string(self):
        ins = [{"after_ref_pos": 2, "inserted_seq": "KL"}]
        assert aa_diff._collect_mutations([], ins) == ["ins2:KL"]

    def test_frameshift_strings_both_signs(self):
        fs = [{"ref_pos": 10, "type": "+1"}, {"ref_pos": 20, "type": "-2"}]
        assert aa_diff._collect_mutations([], [], fs) == ["fs10+1", "fs20-2"]

    def test_ordering_subs_dels_then_ins_then_fs(self):
        pos = [
            {"change_type": "substitution", "ref_aa": "C", "ref_pos": 2, "query_aa": "G"},
            {"change_type": "deletion", "ref_aa": "D", "ref_pos": 4, "query_aa": "-"},
        ]
        ins = [{"after_ref_pos": 5, "inserted_seq": "K"}]
        fs = [{"ref_pos": 8, "type": "+1"}]
        assert aa_diff._collect_mutations(pos, ins, fs) == ["C2G", "D4del", "ins5:K", "fs8+1"]


class TestExtendOverClips:

    def test_no_clip_appends_stop_codon(self):
        prot = "MKLAV"
        coding = cds_of(prot)
        query_nt = coding + "TAA"
        hit = make_hit(0, 5, "+", 0, 15)
        aa, cds = aa_diff._extend_over_clips(prot, hit, 5, query_nt, [])
        assert aa == "MKLAV"
        assert cds == coding + "TAA"

    def test_intron_is_spliced_out_of_cds(self):
        prot = "MKLAV"
        coding = cds_of(prot)                      # 15 nt
        intron = "GTAAAAAAG"                        # 9 nt, arbitrary content
        query_nt = coding[:6] + intron + coding[6:] + "TAA"
        hit = make_hit(0, 5, "+", 0, 6 + 9 + 9)     # nt_end = 24
        aa, cds = aa_diff._extend_over_clips(prot, hit, 5, query_nt, [(6, 9)])
        assert aa == "MKLAV"
        assert cds == coding + "TAA"                # intron gone, stop appended
        assert str(Seq(cds).translate()).rstrip("*") == prot

    def test_cterminal_clip_is_recovered(self):
        # reference has 6 residues; miniprot only aligned the first 5
        full = "MKLAVF"
        coding = cds_of(full)                       # 18 nt
        query_nt = coding + "TAA"
        hit = make_hit(0, 5, "+", 0, 15)            # prot_end 5 < ref_len 6
        aa, cds = aa_diff._extend_over_clips("MKLAV", hit, 6, query_nt, [])
        assert aa == "MKLAVF"
        assert cds == coding + "TAA"

    def test_reverse_strand_returns_coding_orientation(self):
        coding = cds_of("MKLAV") + "TAA"            # protein-orientation CDS, 18 nt
        query_nt = str(Seq(coding).reverse_complement())
        # aligned region (5 codons) maps to [3, 18] on the forward query
        hit = make_hit(0, 5, "-", 3, 18, nt_len=len(query_nt))
        aa, cds = aa_diff._extend_over_clips("MKLAV", hit, 5, query_nt, [])
        assert aa == "MKLAV"
        assert cds == coding


def _hit_dict(**over):
    """a parsed-hit dict with sensible defaults for exercising the table writer"""
    h = {
        "query_seq": "q", "nt_start": 0, "nt_end": 1137, "strand": "+",
        "score": 700, "aln_aa_len": 379, "perc_id": 100.0, "perc_ref_cov": 100.0,
        "n_total_mut": 0, "n_sub": 0, "n_del_runs": 0, "n_del_aa": 0,
        "n_ins": 0, "n_ins_aa": 0, "n_frameshifts": 0, "n_stops": 0,
        "introns": [], "cds_seq": "ATGAAA", "translated_aa": "MK",
    }
    h.update(over)
    return h


class TestWriteMiniprotHitsTable:

    def _read(self, path):
        rows = [line.split("\t") for line in path.read_text().splitlines()]
        header = rows[0]
        return header, [dict(zip(header, r)) for r in rows[1:]]

    def test_header_and_row_count(self, tmp_path, capsys):
        hits = [_hit_dict(), _hit_dict(strand="-", score=300)]
        out = tmp_path / "alignment-summaries.tsv"
        aa_diff._write_miniprot_hits_table(hits, str(out))

        header, data = self._read(out)
        assert header == ["score_rank", "query_seq", "nt_start", "nt_end", "strand",
                          "score", "aln_aa_len", "perc_id", "perc_ref_cov",
                          "total_mutations", "substitutions", "deletions", "insertions",
                          "frameshifts", "introns", "stops", "inferred_cds", "inferred_protein"]
        assert len(data) == 2
        # the writer enumerates in the order it's given (ranking happens upstream)
        assert [d["score_rank"] for d in data] == ["1", "2"]
        # and prints a notice naming the count
        assert "2 miniprot alignments found" in capsys.readouterr().out

    def test_compound_columns_are_formatted(self, tmp_path):
        hits = [_hit_dict(
            n_del_runs=1, n_del_aa=2, n_ins=1, n_ins_aa=3, n_frameshifts=2,
            introns=[(570, 130), (900, 88)], perc_id=98.7, perc_ref_cov=95.4,
            cds_seq="ATGAAACTT", translated_aa="MKL",
        )]
        out = tmp_path / "tbl.tsv"
        aa_diff._write_miniprot_hits_table(hits, str(out))
        _, data = self._read(out)
        row = data[0]
        assert row["deletions"] == "1 (2 AAs total)"
        assert row["insertions"] == "1 (3 AAs total)"
        assert row["introns"] == "2 (218 NTs total)"      # 130 + 88
        assert row["frameshifts"] == "2"
        assert row["perc_id"] == "98.7"
        assert row["perc_ref_cov"] == "95.4"
        assert row["inferred_cds"] == "ATGAAACTT"
        assert row["inferred_protein"] == "MKL"

    def test_no_introns_renders_zero(self, tmp_path):
        out = tmp_path / "tbl.tsv"
        aa_diff._write_miniprot_hits_table([_hit_dict(introns=[])], str(out))
        _, data = self._read(out)
        assert data[0]["introns"] == "0 (0 NTs total)"


class TestRunAaDiffProtein:

    def _write_fasta(self, path, seq_id, seq):
        with open(path, "w") as f:
            f.write(f">{seq_id}\n{seq}\n")

    def test_writes_outputs_and_finds_substitution(self, tmp_path):
        ref = "ACDEFGHIKLMNPQRSTVWY"
        query = ref[:4] + "K" + ref[5:]      # single substitution at position 5 (G->K)
        ref_fa = tmp_path / "ref.faa"
        qry_fa = tmp_path / "query.faa"
        self._write_fasta(str(ref_fa), "ref", ref)
        self._write_fasta(str(qry_fa), "query", query)

        aa_diff.run_aa_diff(str(qry_fa), str(ref_fa), "prot", str(tmp_path))

        for name in ("all-positions.tsv", "mutations.txt", "summary.txt", "alignment.txt"):
            assert (tmp_path / name).exists()

        muts = (tmp_path / "mutations.txt").read_text().split()
        assert muts == [f"{ref[4]}5{query[4]}"]   # e.g. F5K (exact letters depend on ref)

    def test_output_prefix_is_applied(self, tmp_path):
        ref = "ACDEFGHIKLMNPQRSTVWY"
        ref_fa = tmp_path / "ref.faa"
        qry_fa = tmp_path / "query.faa"
        self._write_fasta(str(ref_fa), "ref", ref)
        self._write_fasta(str(qry_fa), "query", ref)   # identical -> no mutations

        aa_diff.run_aa_diff(str(qry_fa), str(ref_fa), "prot", str(tmp_path), output_prefix="sample1-")

        assert (tmp_path / "sample1-summary.txt").exists()
        assert (tmp_path / "sample1-mutations.txt").read_text().startswith("# No mutations")


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
class TestMiniprotIntegration:

    REF = (
        "MSEPLDLNQLAQKIKQWGLELGFQQVGITDTDLSESEPKLQAWLDKQYHGEMDWMARHGMLRARPHELL"
        "PGTLRVISVRMNYLPANAAFASTLKNPKLGYVSRYALGRDYHKLLRNRLKKLGEMIQQHCVSLNFRPFV"
        "DSAPILERPLAAKAGLGWTGKHSLILNREAGSFFFLGELLVDIPLPVDQPVEEGCGKCIACMTICPTGA"
        "IVEPYTVDARRCISYLTIELEGAIPEELRPLMGNRIYGCDDCQLICPWNRYSQLTTEDDFSPRKPLHAP"
        "ELIELFAWSEEKFLKVTEGSAIRRIGHLRWLRNIAVALGNAPWDETILAALESRKGEHPLLDEHIAWAM"
        "AQQIERRNACIVEVQLPKKQRLVRVIEKGLPRDA"
    )

    def _write(self, path, seq_id, seq):
        with open(path, "w") as f:
            f.write(f">{seq_id}\n{seq}\n")

    def test_clean_nt_query_reconstructs_protein(self, tmp_path):
        ref_fa = tmp_path / "ref.faa"
        qry_fa = tmp_path / "query.fna"
        self._write(str(ref_fa), "ref", self.REF)
        self._write(str(qry_fa), "q", cds_of(self.REF) + "TAA")

        aa_diff.run_aa_diff(str(qry_fa), str(ref_fa), "nt", str(tmp_path))

        inferred = (tmp_path / "inferred-protein.faa").read_text().splitlines()[1]
        assert inferred == self.REF
        cds = (tmp_path / "inferred-cds.fna").read_text().splitlines()[1]
        assert str(Seq(cds).translate()).rstrip("*") == self.REF

    def test_intron_query_splices_and_reports(self, tmp_path):
        coding = cds_of(self.REF)
        intron = "GT" + "TAATAGTGA" * 14 + "AG"     # 128 nt, canonical GT..AG
        cut = (len(self.REF) // 2) * 3              # splice at a codon boundary
        query = coding[:cut] + intron + coding[cut:] + "TAA"
        ref_fa = tmp_path / "ref.faa"
        qry_fa = tmp_path / "query.fna"
        self._write(str(ref_fa), "ref", self.REF)
        self._write(str(qry_fa), "q", query)

        aa_diff.run_aa_diff(str(qry_fa), str(ref_fa), "nt", str(tmp_path))

        inferred = (tmp_path / "inferred-protein.faa").read_text().splitlines()[1]
        assert inferred == self.REF
        cds = (tmp_path / "inferred-cds.fna").read_text().splitlines()[1]
        assert len(cds) == len(self.REF) * 3 + 3   # spliced CDS + stop, intron removed
        assert "Spliced-out introns" in (tmp_path / "summary.txt").read_text()

    def test_multi_hit_writes_summaries_table(self, tmp_path):
        coding = cds_of(self.REF)
        # two copies of the CDS -> two loci -> two hits -> a summaries table
        query = coding + "TAA" + "N" * 60 + coding + "TAA"
        ref_fa = tmp_path / "ref.faa"
        qry_fa = tmp_path / "query.fna"
        self._write(str(ref_fa), "ref", self.REF)
        self._write(str(qry_fa), "q", query)

        aa_diff.run_aa_diff(str(qry_fa), str(ref_fa), "nt", str(tmp_path))

        table = tmp_path / "alignment-summaries.tsv"
        assert table.exists()
        rows = [line.split("\t") for line in table.read_text().splitlines()]
        assert len(rows) == 3                       # header + 2 hits
        assert [r[0] for r in rows[1:]] == ["1", "2"]
        # the primary outputs for the top hit are still produced
        assert (tmp_path / "inferred-protein.faa").exists()
        assert (tmp_path / "summary.txt").exists()


class TestWriteTsv:

    def _rows(self, path):
        lines = path.read_text().splitlines()
        header = lines[0].split("\t")
        return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_insertion_anchor_matches_mutations_file(self, tmp_path):
        # insertion after ref residue 3 -> mutations.txt writes ins3, so TSV row 3
        positions, insertions = aa_diff._parse_alignment("ACD--EFG", "ACDMMEFG")
        out = tmp_path / "all-positions.tsv"
        aa_diff._write_tsv(positions, insertions, str(out))
        header, rows = self._rows(out)
        assert "inserted_after" in header
        by_pos = {r["ref_pos"]: r["inserted_after"] for r in rows}
        assert by_pos["3"] == "MM"      # same anchor number as ins3:MM
        assert by_pos["4"] == "-"

    def test_trailing_insertion_is_visible(self, tmp_path):
        positions, insertions = aa_diff._parse_alignment("ACDEFG--", "ACDEFGMM")
        out = tmp_path / "all-positions.tsv"
        aa_diff._write_tsv(positions, insertions, str(out))
        _, rows = self._rows(out)
        assert rows[-1]["inserted_after"] == "MM"   # after the last residue (ins6)
