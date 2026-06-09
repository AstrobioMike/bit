import pandas as pd  # type: ignore
import pytest # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation  # type: ignore

from bit.modules.genbank import (
    genbank_to_fasta,
    genbank_to_AA_seqs,
    genbank_to_cds_seqs,
    genbank_to_cds_tsv,
    parse_genbank_cds_to_dataframe,
)


def _cds(start, end, strand=1, qualifiers=None, location=None):
    loc = location if location is not None else FeatureLocation(start, end, strand=strand)
    return SeqFeature(loc, type="CDS", qualifiers=qualifiers or {})


def make_record():
    # 60 bp so the simple CDS (0-30) extracts a clean 10-codon ORF
    seq = Seq("ATG" + "AAA" * 8 + "TAA" + "ATGCGT" * 5 + "GGG")
    seq = seq[:60] if len(seq) >= 60 else seq + Seq("A" * (60 - len(seq)))

    features = [
        # 1. clean, complete CDS with all qualifiers -> kept everywhere
        _cds(0, 30, qualifiers={
            "translation": ["MKKKKKKKK"],
            "locus_tag": ["TAG_0001"],
            "gene": ["abcA"],
            "protein_id": ["PROT_001.1"],
            "product": ["alpha subunit"],   # space -> hyphen in FASTA header
        }),
        # 2. note says frameshifted -> excluded from AA/CDS fasta
        _cds(30, 45, qualifiers={
            "translation": ["MKKK"],
            "locus_tag": ["TAG_0002"],
            "note": ["frameshifted internal region"],
        }),
        # 3. transl_except -> excluded from AA/CDS fasta
        _cds(0, 15, qualifiers={
            "translation": ["MKK"],
            "locus_tag": ["TAG_0003"],
            "transl_except": ["(pos:1..3,aa:Sec)"],
        }),
        # 4. pseudo -> excluded from AA/CDS fasta
        _cds(0, 15, qualifiers={
            "translation": ["MKK"],
            "locus_tag": ["TAG_0004"],
            "pseudo": [""],
        }),
        # 5. compound (join) location -> excluded from AA/CDS fasta
        _cds(None, None, qualifiers={
            "translation": ["MKK"],
            "locus_tag": ["TAG_0005"],
        }, location=CompoundLocation([FeatureLocation(0, 9, 1),
                                      FeatureLocation(12, 18, 1)])),
        # 6. CDS without a translation -> ignored by AA/CDS fasta,
        #    but still counted by the dataframe parser
        _cds(0, 15, qualifiers={"locus_tag": ["TAG_0006"]}),
        # 7. non-CDS feature -> ignored everywhere
        SeqFeature(FeatureLocation(0, 60, 1), type="gene",
                   qualifiers={"locus_tag": ["GENE_X"]}),
    ]

    return SeqRecord(seq, id="contig1", name="contig1", description="test",
                     annotations={"molecule_type": "DNA"}, features=features)


@pytest.fixture
def genbank_file(tmp_path):
    from Bio import SeqIO  # type: ignore
    path = tmp_path / "test.gbk"
    SeqIO.write(make_record(), str(path), "genbank")
    return path


# ───────────────────────── genbank_to_fasta ─────────────────────────

def test_genbank_to_fasta_writes_nucleotide(genbank_file, tmp_path):
    out = tmp_path / "out.fasta"
    genbank_to_fasta(str(genbank_file), str(out))
    text = out.read_text()
    # one record, headed by record.name, full contig sequence
    assert text.startswith(">contig1\n")
    assert text.count(">") == 1
    assert "ATG" in text


# ───────────────────────── genbank_to_AA_seqs ─────────────────────────

def test_genbank_to_AA_seqs_keeps_only_clean_cds(genbank_file, tmp_path):
    out = tmp_path / "aa.fasta"
    genbank_to_AA_seqs(str(genbank_file), str(out))
    text = out.read_text()

    # only feature #1 survives all the filters
    assert text.count(">") == 1
    assert ">TAG_0001|PROT_001.1|abcA|alpha-subunit" in text   # spaces -> hyphens
    assert "MKKKKKKKK" in text

    # excluded ones absent
    for tag in ("TAG_0002", "TAG_0003", "TAG_0004", "TAG_0005"):
        assert tag not in text


def test_genbank_to_AA_seqs_qualifier_fallbacks(tmp_path):
    from Bio import SeqIO  # type: ignore
    # a clean CDS missing gene/protein_id/product to hit the .get() defaults
    rec = SeqRecord(
        Seq("ATGAAATAA" + "A" * 21), id="c", name="c",
        annotations={"molecule_type": "DNA"},
        features=[_cds(0, 30, qualifiers={
            "translation": ["MK"],
            "locus_tag": ["ONLY_TAG"],
        })],
    )
    gbk = tmp_path / "fb.gbk"
    SeqIO.write(rec, str(gbk), "genbank")
    out = tmp_path / "fb.faa"
    genbank_to_AA_seqs(str(gbk), str(out))
    assert ">ONLY_TAG|No_protein_id|No_gene|No_product" in out.read_text()


# ───────────────────────── genbank_to_cds_seqs ─────────────────────────

def test_genbank_to_cds_seqs_extracts_nucleotide(genbank_file, tmp_path):
    out = tmp_path / "cds.fasta"
    genbank_to_cds_seqs(str(genbank_file), str(out))
    text = out.read_text()

    assert text.count(">") == 1
    assert ">TAG_0001|PROT_001.1|abcA|alpha-subunit" in text
    # the extracted nucleotide CDS (not the protein) follows the header
    seq_line = text.splitlines()[1]
    assert set(seq_line) <= set("ACGT")
    assert seq_line.startswith("ATG")


# ─────────────────── parse_genbank_cds_to_dataframe ───────────────────

def test_parse_cds_dataframe_counts_all_cds(genbank_file):
    df = parse_genbank_cds_to_dataframe(str(genbank_file))
    assert isinstance(df, pd.DataFrame)
    # the dataframe parser does NOT apply the exclusion filters, so every
    # CDS feature is present (features 1-6); the non-CDS gene is excluded
    assert len(df) == 6
    assert list(df.columns) == ["gene", "protein_id", "locus_tag", "product"]
    assert "TAG_0001" in set(df["locus_tag"])
    assert "GENE_X" not in set(df["locus_tag"])


def test_parse_cds_dataframe_uses_na_defaults(genbank_file):
    df = parse_genbank_cds_to_dataframe(str(genbank_file))
    row = df[df["locus_tag"] == "TAG_0006"].iloc[0]
    assert row["gene"] == "NA"
    assert row["protein_id"] == "NA"
    assert row["product"] == "NA"


# ───────────────────────── genbank_to_cds_tsv ─────────────────────────

def test_genbank_to_cds_tsv_writes_tsv(genbank_file, tmp_path):
    out = tmp_path / "cds.tsv"
    genbank_to_cds_tsv(str(genbank_file), str(out))
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 6
    assert list(df.columns) == ["gene", "protein_id", "locus_tag", "product"]
