from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
from pathlib import Path
import shutil
import pandas as pd # type: ignore
import pytest # type: ignore
from bit.modules import ez_screen as ez


test_targets_fasta = get_package_path("tests/data/ez-screen-targets.fasta")
test_assembly_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")


# =========================================================================== #
# CLI end-to-end
# =========================================================================== #

def test_ez_screen_assembly(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    cmd = [
        "bit", "ez-screen",
        "assembly",
        "-a", str(test_assembly_fasta),
        "-t", str(test_targets_fasta),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-summary.tsv")

    observed = summary_tsv.read_text().splitlines()
    expected = [
        "target\tez-screen-assembly",
        "yopE\t1"
    ]

    assert observed == expected, f"Summary TSV content does not match expected:\n{observed}\nExpected:\n{expected}"


def test_ez_screen_assembly_assemblies_as_rows(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    cmd = [
        "bit", "ez-screen",
        "assembly",
        "-a", str(test_assembly_fasta),
        "-t", str(test_targets_fasta),
        "-o", str(out_prefix),
        "--assemblies-as-rows"
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-summary.tsv")

    observed = summary_tsv.read_text().splitlines()
    expected = [
        "input-assembly\tyopE",
        "ez-screen-assembly\t1"
    ]

    assert observed == expected, f"Summary TSV content does not match expected:\n{observed}\nExpected:\n{expected}"


def test_ez_screen_reads(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    R1 = get_package_path("tests/data/ez-screen-R1.fastq.gz")
    R2 = get_package_path("tests/data/ez-screen-R2.fastq.gz")

    shutil.copy(R1, reads_dir)
    shutil.copy(R2, reads_dir)

    cmd = [
        "bit", "ez-screen",
        "reads",
        "-t", str(test_targets_fasta),
        "-r", str(reads_dir),
        "-o", str(out_prefix),
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-reads-summary.tsv")

    observed = summary_tsv.read_text().splitlines()
    expected = [
        "sample\ttarget\tnum_reads_recruited\tdetection\tmean_perc_id",
        "ez-screen\tyopE\t20\t0.9\t100.0"
    ]

    assert observed == expected, f"Summary TSV content does not match expected:\n{observed}\nExpected:\n{expected}"


# =========================================================================== #
# reads-path units
# =========================================================================== #

class _Read:

    def __init__(self, qlen, nm, ref_id, flags=0):
        self.query_length = qlen
        self._nm = nm
        self.reference_id = ref_id
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False
        self.cigartuples = [(0, qlen)]

    def get_tag(self, _):
        return self._nm


class _Bam:

    def __init__(self, reads):
        self._reads = reads
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc, tb):
        return False
    def fetch(self, until_eof=True):
        for r in self._reads:
            yield r
    def get_reference_name(self, rid):
        return ["yopE","yopK"][rid]


def test_gen_reads_summary_table_outputs(tmp_path, monkeypatch):
    # fake BAM: 3 reads to yopE: nm=0 on qlen=100 => 100% pid
    reads = [_Read(100,0,0) for _ in range(3)]
    monkeypatch.setattr(ez.pysam, "AlignmentFile", lambda *a, **k: _Bam(reads))

    # detection table keeps yopE, drops others/total
    gtab = tmp_path/"global.tab"
    with open(gtab, "w") as f:
        f.write("target\tdepth\tdetection\n")
        f.write("yopE\t1\t0.95\n")
        f.write("total\t1\t0.98\n")

    out = tmp_path/"summary.tsv"
    ez.gen_reads_summary_table(
        input_bam="fake.bam",
        input_global_dist_tab=gtab,
        outpath=out,
        min_perc_id=90.0,
        min_perc_cov=90.0,
    )
    lines = out.read_text().splitlines()
    assert lines[0] == "target\tnum_reads_recruited\tdetection\tmean_perc_id"
    assert lines[1].startswith("yopE\t3\t0.95\t100.00")


def test_gen_reads_summary_table_empty(tmp_path, monkeypatch):
    # No reads retained by coverage table
    monkeypatch.setattr(ez.pysam, "AlignmentFile", lambda *a, **k: _Bam([]))
    gtab = tmp_path/"global.tab"; gtab.write_text("target\tdepth\tdetection\n")
    out = tmp_path/"summary.tsv"
    ez.gen_reads_summary_table("fake.bam", gtab, out, 99.0, 99.0)
    assert "No reads successfully mapped" in out.read_text()


def test_combine_reads_summary_outputs_mixed(tmp_path):
    # two valid per-sample summaries
    good1 = tmp_path / "s1.tsv"
    good2 = tmp_path / "s2.tsv"
    pd.DataFrame([{
        "target": "yopE", "num_reads_recruited": 10, "detection": 0.9, "mean_perc_id": 100.0
    }]).to_csv(good1, sep="\t", index=False)
    pd.DataFrame([{
        "target": "yopE", "num_reads_recruited": 20, "detection": 0.9, "mean_perc_id": 99.5
    }]).to_csv(good2, sep="\t", index=False)

    # skip cases: bad columns, empty table, unparsable (directory)
    bad_cols = tmp_path / "bad_cols.tsv"
    pd.DataFrame([{"foo": 1}]).to_csv(bad_cols, sep="\t", index=False)
    empty = tmp_path / "empty.tsv"
    pd.DataFrame(columns=["target","num_reads_recruited","detection","mean_perc_id"]).to_csv(empty, sep="\t", index=False)
    bad_dir = tmp_path / "bad_dir"
    bad_dir.mkdir()

    samples = {
        "s1": str(good1),
        "s2": str(good2),
        "bad": str(bad_cols),
        "empty": str(empty),
        "oops": str(bad_dir),  # triggers exception path
    }
    out = tmp_path / "combined.tsv"

    ez.combine_reads_summary_outputs(samples, out)

    df = pd.read_csv(out, sep="\t")
    assert df.columns.tolist() == ["sample","target","num_reads_recruited","detection","mean_perc_id"]
    assert df.shape[0] == 2
    assert df["sample"].tolist() == ["s1","s2"]
    assert df.loc[0, "num_reads_recruited"] == 10
    assert df.loc[1, "num_reads_recruited"] == 20


def test_combine_reads_summary_outputs_no_valid(tmp_path):
    # only invalid inputs -> fallback message
    bad_cols = tmp_path / "bad.tsv"
    pd.DataFrame([{"foo": 1}]).to_csv(bad_cols, sep="\t", index=False)
    empty = tmp_path / "empty.tsv"
    pd.DataFrame(columns=["target","num_reads_recruited","detection","mean_perc_id"]).to_csv(empty, sep="\t", index=False)

    out = tmp_path / "combined.tsv"
    ez.combine_reads_summary_outputs({"bad": str(bad_cols), "empty": str(empty)}, out)

    assert "No valid read-mappings were found to any targets." in out.read_text()


# =========================================================================== #
# helpers (loci architecture)
# =========================================================================== #

def make_hsp(qseqid, sseqid, qstart, qend, sstart, send, slen,
             pident=99.0, bitscore=None, qlen=15000):
    """ one HSP row shaped like a (pident-passing) BLAST hit. bitscore defaults
        to the aligned length so longer/contained HSPs win by default. """
    length = abs(qend - qstart) + 1
    if bitscore is None:
        bitscore = length
    return {
        "qseqid": qseqid, "qlen": qlen, "sseqid": sseqid, "slen": slen,
        "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
        "length": length, "qcovs": 0, "qcovhsp": 0, "qcovus": 0,
        "pident": pident, "evalue": 0.0, "bitscore": bitscore,
        "perc-subj-cov": 0.0,
    }


def loci_from(rows, targets_dict, gap_tolerance=200, min_perc_cov=80.0):
    """ run the engine on a list of HSP rows """
    return ez.resolve_assembly_loci(
        pd.DataFrame(rows), targets_dict,
        gap_tolerance=gap_tolerance, min_perc_cov=min_perc_cov)


EMPTY_HSP_COLS = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send",
                  "slen", "pident", "bitscore", "length", "qlen"]


def empty_loci():
    return ez.resolve_assembly_loci(pd.DataFrame(columns=EMPTY_HSP_COLS), {"t": 1000})


LOCI_COLS = ["contig", "target", "q_low", "q_high", "perc_target_cov",
             "pident", "bitscore", "length", "slen"]

REGION_COLS = ["contig", "region_start", "region_end", "region_length",
               "aligned_target", "target_length", "pident", "perc_target_cov",
               "bitscore", "n_overlapping_targets", "other_targets"]


# =========================================================================== #
# resolve_assembly_loci : the engine
# =========================================================================== #

def test_engine_empty():
    out = empty_loci()
    assert out.empty
    assert list(out.columns) == LOCI_COLS


def test_engine_fragmented_single_occurrence_counts_once():
    rows = [
        make_hsp("c1", "t", 100, 579, 1, 480, 1000),
        make_hsp("c1", "t", 600, 1090, 510, 1000, 1000),
    ]
    loci = loci_from(rows, {"t": 1000})
    assert len(loci) == 1
    assert loci.iloc[0]["perc_target_cov"] == 97.1
    assert loci.iloc[0]["q_low"] == 100
    assert loci.iloc[0]["q_high"] == 1090


def test_engine_long_backbone_fragmented_counts_once():
    rows = [
        make_hsp("c1", "bb", 100, 2100, 1, 2000, 12000, bitscore=1000),
        make_hsp("c1", "bb", 2200, 4200, 2100, 4100, 12000, bitscore=1000),
        make_hsp("c1", "bb", 4300, 6300, 4200, 6200, 12000, bitscore=1000),
        make_hsp("c1", "bb", 6400, 8400, 6300, 8300, 12000, bitscore=1000),
        make_hsp("c1", "bb", 8500, 10100, 8400, 10000, 12000, bitscore=1100),
    ]
    loci = loci_from(rows, {"bb": 12000})
    assert len(loci) == 1
    assert loci.iloc[0]["perc_target_cov"] == 80.0
    assert loci.iloc[0]["bitscore"] == 1100


def test_engine_incidental_low_coverage_dropped():
    rows = [make_hsp("c1", "bb", 100, 300, 4000, 4200, 12000)]
    assert loci_from(rows, {"bb": 12000}).empty


def test_engine_two_loci_same_target_separate_places():
    rows = [
        make_hsp("c1", "t", 100, 1000, 1, 900, 1000),
        make_hsp("c1", "t", 5000, 5900, 1, 900, 1000),
    ]
    assert len(loci_from(rows, {"t": 1000})) == 2


def test_engine_same_target_two_contigs():
    rows = [
        make_hsp("c1", "t", 100, 1000, 1, 900, 1000),
        make_hsp("c2", "t", 100, 1000, 1, 900, 1000),
    ]
    loci = loci_from(rows, {"t": 1000})
    assert len(loci) == 2
    assert set(loci["contig"]) == {"c1", "c2"}


def test_engine_gap_tolerance_controls_locus_joining():
    rows = [
        make_hsp("c1", "t", 100, 549, 1, 450, 900),
        make_hsp("c1", "t", 700, 1149, 451, 900, 900),
    ]
    assert len(loci_from(rows, {"t": 900}, gap_tolerance=200)) == 1
    assert len(loci_from(rows, {"t": 900}, gap_tolerance=100)) == 0


def test_engine_subject_overlap_no_double_count():
    rows = [
        make_hsp("c1", "t", 100, 700, 1, 600, 1000),
        make_hsp("c1", "t", 120, 720, 1, 600, 1000),
    ]
    loci = loci_from(rows, {"t": 1000}, min_perc_cov=50)
    assert len(loci) == 1
    assert loci.iloc[0]["perc_target_cov"] == 60.0


def test_engine_unknown_target_skipped():
    rows = [make_hsp("c1", "ghost", 100, 700, 1, 600, 600)]
    assert loci_from(rows, {"t": 1000}).empty


def test_engine_coverage_threshold_boundary():
    rows = [make_hsp("c1", "t", 100, 899, 1, 800, 1000)]
    assert len(loci_from(rows, {"t": 1000}, min_perc_cov=80)) == 1
    rows2 = [make_hsp("c1", "t", 100, 898, 1, 799, 1000)]
    assert len(loci_from(rows2, {"t": 1000}, min_perc_cov=80)) == 0


# =========================================================================== #
# _reciprocal_overlap_frac
# =========================================================================== #

def test_reciprocal_overlap_full_containment():
    assert ez._reciprocal_overlap_frac(100, 700, 200, 400) == 1.0


def test_reciprocal_overlap_identical_spans():
    assert ez._reciprocal_overlap_frac(100, 200, 100, 200) == 1.0


def test_reciprocal_overlap_disjoint():
    assert ez._reciprocal_overlap_frac(100, 200, 300, 400) == 0.0


def test_reciprocal_overlap_one_base_touch():
    # 1-based inclusive: sharing a single boundary coordinate is a 1 bp overlap.
    # shorter span = 201 bp -> 1/201
    assert ez._reciprocal_overlap_frac(100, 300, 300, 500) == pytest.approx(1 / 201)


def test_reciprocal_overlap_exactly_half():
    assert ez._reciprocal_overlap_frac(1, 200, 101, 400) == pytest.approx(0.5)


def test_reciprocal_overlap_uses_shorter_span():
    # denominator must be the SHORTER span: overlap 51 / shorter 200 = 0.255
    assert ez._reciprocal_overlap_frac(1, 200, 150, 1000) == pytest.approx(51 / 200)


# =========================================================================== #
# resolve_contig_regions  (operates on loci for ONE contig)
# =========================================================================== #

def test_regions_empty():
    assert ez.resolve_contig_regions(empty_loci()) == []


def test_regions_single_locus():
    rows = [make_hsp("c", "pUC19_ori", 100, 700, 1, 600, 600, bitscore=1100)]
    loci = loci_from(rows, {"pUC19_ori": 600})
    regions = ez.resolve_contig_regions(loci)
    assert len(regions) == 1
    r = regions[0]
    assert r["aligned_target"] == "pUC19_ori"
    assert r["n_overlapping_targets"] == 1
    assert r["other_targets"] == []
    assert r["perc_target_cov"] == 100.0


def test_regions_redundant_variants_collapse():
    rows = [
        make_hsp("c", "pUC19_ori", 100, 700, 1, 600, 600, bitscore=1100),
        make_hsp("c", "pBR322_ori", 102, 698, 1, 597, 600, bitscore=1050),
    ]
    loci = loci_from(rows, {"pUC19_ori": 600, "pBR322_ori": 600})
    regions = ez.resolve_contig_regions(loci)
    assert len(regions) == 1
    r = regions[0]
    assert r["aligned_target"] == "pUC19_ori"
    assert r["other_targets"] == ["pBR322_ori"]
    assert r["n_overlapping_targets"] == 2
    assert r["perc_target_cov"] == 100.0


def test_regions_two_loci_same_target_stay_separate():
    rows = [
        make_hsp("c", "IS", 100, 400, 1, 301, 301),
        make_hsp("c", "IS", 5000, 5300, 1, 301, 301),
    ]
    loci = loci_from(rows, {"IS": 301})
    assert len(ez.resolve_contig_regions(loci)) == 2


def test_regions_overlap_frac_boundary_inclusive():
    # exactly 0.5 reciprocal overlap must JOIN (the test is >=).
    # two distinct targets at one locus: spans 1..200 (200 bp) and 101..400,
    # overlap 100 bp / shorter 200 = 0.5. both must clear coverage to appear,
    # so give each full subject coverage of its own target.
    rows = [
        make_hsp("c", "a", 1, 200, 1, 200, 200, bitscore=500),
        make_hsp("c", "b", 101, 400, 1, 300, 300, bitscore=600),
    ]
    loci = loci_from(rows, {"a": 200, "b": 300}, min_perc_cov=80)
    regions = ez.resolve_contig_regions(loci, overlap_frac=0.5)
    assert len(regions) == 1
    assert regions[0]["n_overlapping_targets"] == 2


def test_regions_bitscore_tie_broken_by_pident_length():
    # equal bitscore; tie-break is pident * alignment-length. tie_long has the
    # longer alignment (601 vs 381), so 99.0*601 beats 99.9*381 despite lower
    # pident. spans overlap enough (>=50% of shorter) to cluster.
    rows = [
        make_hsp("c", "tie_short", 100, 480, 1, 380, 400, pident=99.9, bitscore=1000),
        make_hsp("c", "tie_long",  100, 700, 1, 601, 601, pident=99.0, bitscore=1000),
    ]
    loci = loci_from(rows, {"tie_short": 400, "tie_long": 601}, min_perc_cov=80)
    regions = ez.resolve_contig_regions(loci)
    assert len(regions) == 1
    assert regions[0]["aligned_target"] == "tie_long"


# =========================================================================== #
# gen_region_calls_table
# =========================================================================== #

def _multi_contig_loci():
    rows = [
        make_hsp("cA", "pUC19_ori", 100, 700, 1, 600, 600, bitscore=1100),
        make_hsp("cA", "pBR322_ori", 102, 698, 1, 597, 600, bitscore=1050),
        make_hsp("cA", "kanR", 3000, 3800, 1, 800, 800, bitscore=1400),
        make_hsp("cB", "IS", 100, 400, 1, 301, 301),
        make_hsp("cB", "IS", 5000, 5300, 1, 301, 301),
    ]
    return loci_from(rows, {"pUC19_ori": 600, "pBR322_ori": 600,
                            "kanR": 800, "IS": 301})


def test_region_calls_columns_and_content():
    loci = _multi_contig_loci()
    region_df, num_regions = ez.gen_region_calls_table(loci, overlap_frac=0.5)
    assert list(region_df.columns) == REGION_COLS
    assert num_regions == {"cA": 2, "cB": 2}
    assert len(region_df) == 4
    assert region_df["contig"].tolist() == ["cA", "cA", "cB", "cB"]


def test_region_calls_perc_target_cov_is_winner_coverage():
    loci = _multi_contig_loci()
    region_df, _ = ez.gen_region_calls_table(loci)
    ori = region_df[(region_df["contig"] == "cA") &
                    (region_df["region_start"] == 100)].iloc[0]
    assert ori["aligned_target"] == "pUC19_ori"
    assert ori["perc_target_cov"] == 100.0
    assert ori["other_targets"] == "pBR322_ori"


def test_region_calls_other_targets_none_for_singletons():
    loci = _multi_contig_loci()
    region_df, _ = ez.gen_region_calls_table(loci)
    kan = region_df[region_df["aligned_target"] == "kanR"].iloc[0]
    assert kan["other_targets"] == "none"


def test_region_calls_empty():
    region_df, num_regions = ez.gen_region_calls_table(empty_loci())
    assert num_regions == {}
    assert len(region_df) == 0
    assert list(region_df.columns) == REGION_COLS


# =========================================================================== #
# gen_contig_summary_table  (from loci)
# =========================================================================== #

def _contig_lengths():
    return {"cA": 15000, "cB": 15000}


def test_contig_summary_contig_column_named():
    cd = ez.gen_contig_summary_table(_multi_contig_loci(), _contig_lengths())
    assert "contig" in cd.columns
    assert "index" not in cd.columns


def test_contig_summary_columns_without_regions():
    cd = ez.gen_contig_summary_table(_multi_contig_loci(), _contig_lengths())
    assert list(cd.columns) == [
        "contig", "length", "bases_aligned_to_targets",
        "perc_contig_aligned_to_targets", "num_unique_target_hits", "num_total_target_hits",
    ]


def test_contig_summary_with_regions_adds_column():
    loci = _multi_contig_loci()
    region_df, num_regions = ez.gen_region_calls_table(loci)
    cd = ez.gen_contig_summary_table(loci, _contig_lengths(), num_regions)
    assert list(cd.columns) == [
        "contig", "length", "bases_aligned_to_targets",
        "perc_contig_aligned_to_targets", "num_regions",
        "num_unique_target_hits", "num_total_target_hits",
    ]
    assert dict(zip(cd["contig"], cd["num_regions"])) == num_regions


def test_contig_summary_num_total_hits_counts_loci():
    cd = ez.gen_contig_summary_table(_multi_contig_loci(), _contig_lengths())
    by_contig = dict(zip(cd["contig"], cd["num_total_target_hits"]))
    assert by_contig["cA"] == 3
    assert by_contig["cB"] == 2


def test_contig_summary_num_unique_can_exceed_num_regions():
    loci = _multi_contig_loci()
    region_df, num_regions = ez.gen_region_calls_table(loci)
    cd = ez.gen_contig_summary_table(loci, _contig_lengths(), num_regions)
    ca = cd[cd["contig"] == "cA"].iloc[0]
    assert ca["num_unique_target_hits"] == 3
    assert ca["num_regions"] == 2


def test_contig_summary_integer_columns_stay_integer():
    # regression: a contig present in contig_lengths but with NO passing loci
    # must not float-ify the integer columns of the contigs that do have loci.
    loci = _multi_contig_loci()
    lengths = dict(_contig_lengths())
    lengths["cEMPTY"] = 50000  # extra contig with no loci -> the NaN trigger
    cd = ez.gen_contig_summary_table(loci, lengths)
    for col in ("length", "bases_aligned_to_targets", "num_unique_target_hits",
                "num_total_target_hits"):
        assert str(cd[col].dtype) == "int64", f"{col} should be int, got {cd[col].dtype}"
    assert "cEMPTY" not in set(cd["contig"])


def test_contig_summary_empty():
    cd_plain = ez.gen_contig_summary_table(empty_loci(), {})
    assert "num_regions" not in cd_plain.columns
    assert len(cd_plain) == 0

    cd_regions = ez.gen_contig_summary_table(empty_loci(), {}, {})
    assert "num_regions" in cd_regions.columns
    assert len(cd_regions) == 0


# =========================================================================== #
# filter_blast_results : pident hard filter
# =========================================================================== #

def test_filter_blast_results_drops_below_pident():
    df = pd.DataFrame([
        make_hsp("c", "t", 100, 700, 1, 600, 600, pident=99.0),
        make_hsp("c", "t", 100, 700, 1, 600, 600, pident=70.0),
    ])
    out = ez.filter_blast_results(df, min_perc_id=80)
    assert len(out) == 1
    assert out.iloc[0]["pident"] == 99.0


# =========================================================================== #
# update_assembly_summary_table : counts of passing loci
# =========================================================================== #

def test_summary_counts_passing_loci_per_target():
    rows = [
        make_hsp("c1", "yopE", 100, 350, 1, 250, 300),
        make_hsp("c1", "yopE", 5000, 5250, 1, 250, 300),
        make_hsp("c1", "yopK", 2000, 2230, 1, 230, 250),
        make_hsp("c1", "bigT", 100, 10200, 1, 10000, 12000),  # ~83% cov
    ]
    loci = loci_from(rows, {"yopE": 300, "yopK": 250, "bigT": 12000})
    out = ez.update_assembly_summary_table(
        loci, {"yopE": 300, "yopK": 250, "bigT": 12000}, "asmX", pd.DataFrame())
    assert out.loc["asmX", "yopE"] == 2
    assert out.loc["asmX", "yopK"] == 1
    assert out.loc["asmX", "bigT"] == 1


def test_summary_undetected_target_is_zero():
    rows = [make_hsp("c1", "yopE", 100, 350, 1, 250, 300)]
    loci = loci_from(rows, {"yopE": 300, "yopK": 250})
    out = ez.update_assembly_summary_table(
        loci, {"yopE": 300, "yopK": 250}, "asmY", pd.DataFrame())
    assert out.loc["asmY", "yopE"] == 1
    assert out.loc["asmY", "yopK"] == 0


def test_summary_fragmented_hsps_count_once():
    rows = [
        make_hsp("c1", "yopE", 100, 250, 1, 150, 300),
        make_hsp("c1", "yopE", 350, 500, 151, 300, 300),
    ]
    loci = loci_from(rows, {"yopE": 300})
    out = ez.update_assembly_summary_table(
        loci, {"yopE": 300}, "asmX", pd.DataFrame())
    assert out.loc["asmX", "yopE"] == 1


def test_summary_empty_loci_all_zero():
    out = ez.update_assembly_summary_table(
        empty_loci(), {"yopE": 300, "yopK": 250}, "asmZ", pd.DataFrame())
    assert out.loc["asmZ", "yopE"] == 0
    assert out.loc["asmZ", "yopK"] == 0


def test_summary_count_matches_contig_num_total_hits():
    # headline consistency: each target's summary count == sum of its
    # num_total_hits across contigs, since both come from the loci.
    rows = [
        make_hsp("c1", "tgt", 100, 350, 1, 250, 300),
        make_hsp("c1", "tgt", 5000, 5250, 1, 250, 300),
        make_hsp("c2", "tgt", 100, 350, 1, 250, 300),
    ]
    loci = loci_from(rows, {"tgt": 300})
    out = ez.update_assembly_summary_table(
        loci, {"tgt": 300}, "asmX", pd.DataFrame())
    summary_count = out.loc["asmX", "tgt"]

    contig_df = ez.gen_contig_summary_table(loci, {"c1": 15000, "c2": 15000})
    contig_total = contig_df["num_total_target_hits"].sum()

    assert summary_count == contig_total == 3


# =========================================================================== #
# filter_undetected_assembly_targets
# =========================================================================== #

def test_filter_undetected_drops_all_zero_columns():
    summary = pd.DataFrame(
        {"yopE": [2, 0], "yopK": [0, 0], "bigT": [1, 0]},
        index=["asmA", "asmB"])
    out = ez.filter_undetected_assembly_targets(summary)
    assert "yopK" not in out.columns
    assert "yopE" in out.columns
    assert "bigT" in out.columns
