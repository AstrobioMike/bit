from bit.modules.general import get_package_path
from bit.tests.utils import run_cli
from pathlib import Path
import shutil
import pandas as pd
from bit.modules import ez_screen as ez


test_targets_fasta = get_package_path("tests/data/ez-screen-targets.fasta")
test_assembly_fasta = get_package_path("tests/data/ez-screen-assembly.fasta")

def test_ez_screen_assembly(tmp_path):
    out_prefix = tmp_path / "ez-screen"
    cmd = [
        "bit-ez-screen",
        "assembly",
        "-a", str(test_assembly_fasta),
        "-t", str(test_targets_fasta),
        "-o", str(out_prefix),
        "--filter-if-not-detected"
    ]

    run_cli(cmd)

    summary_tsv = Path(f"{out_prefix}-assembly-summary.tsv")

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
        "bit-ez-screen",
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


def test_long_targets_merge_and_detect():
    # Two overlapping HSPs on a long target (length 12k) for sseqid "LT"
    df = pd.DataFrame({
        "sseqid":["LT","LT"],
        "sstart":[100,180],
        "send":[200,260],
        "pident":[99.0, 99.0],
    })
    targets = {"LT": 12000}  # long
    out = ez.gen_long_targets_assembly_results_dict(df, targets, long_target_length_cutoff=10000, min_perc_cov=0.5)
    # merged coverage = 100..260 => 161 bp; 0.5% of 12k = 60bp -> DETECTED
    assert out["LT"] == "DETECTED"


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
