import subprocess
import sys

from bit.modules.MAG_abundance import write_MAG_tables
from bit.modules.general import get_package_path

HEADER = "ref\tlength\tnum_contigs\tnum_bases_detected\tdetection (%)\tdetection_at_10x (%)\tmean_coverage\tmedian_coverage\tmean_pid\tmedian_pid\tnum_mapped_reads\n"


def _cov_stats(path, rows):
    path.write_text(HEADER + "".join(f"MAGs/{m}.fasta\t1000\t2\t{d * 10:.0f}\t{d:.2f}\t0.00\t{c:.2f}\t{c:.2f}\t99.00\t99.50\t{n:,}\n"
                                     for m, d, c, n in rows))
    return str(path)


def _mapping_info(path, pairs):
    path.write_text(f"{pairs} reads; of these:\n")
    return str(path)


def test_MAG_by_sample_tables(tmp_path):
    cov = [_cov_stats(tmp_path / "s1.tsv", [("S1-MAG-1", 98.5, 12.25, 2000), ("S2-MAG-3", 0, 0, 0)]),
           _cov_stats(tmp_path / "s2.tsv", [("S1-MAG-1", 40, 0.5, 100), ("S2-MAG-3", 100, 30, 1500)])]
    info = [_mapping_info(tmp_path / "s1.txt", 5000), _mapping_info(tmp_path / "s2.txt", 1000)]
    outs = [tmp_path / f"{w}.tsv" for w in ("det", "pct", "cov")]
    write_MAG_tables(["s1", "s2"], cov, info, True, *outs)

    assert outs[0].read_text() == "MAG_ID\ts1\ts2\nS1-MAG-1\t98.50\t40.00\nS2-MAG-3\t0.00\t100.00\n"
    # reads: mates counted separately (10,000 and 2,000)
    assert outs[1].read_text() == "MAG_ID\ts1\ts2\nS1-MAG-1\t20.00\t5.00\nS2-MAG-3\t0.00\t75.00\n"
    assert outs[2].read_text() == "MAG_ID\ts1\ts2\nS1-MAG-1\t12.25\t0.50\nS2-MAG-3\t0.00\t30.00\n"


def test_MAG_by_sample_tables_without_MAGs(tmp_path):
    (tmp_path / "s1.tsv").touch()
    (tmp_path / "s1.txt").touch()
    outs = [tmp_path / f"{w}.tsv" for w in ("det", "pct", "cov")]
    write_MAG_tables(["s1"], [str(tmp_path / "s1.tsv")], [str(tmp_path / "s1.txt")], True, *outs)
    assert all(o.read_text() == "There were no MAGs recovered.\n" for o in outs)


def test_drep_clusters_script(tmp_path):
    tables = tmp_path / "drep" / "data_tables"
    tables.mkdir(parents=True)
    (tables / "Cdb.csv").write_text("genome,location,length,primary_cluster,secondary_cluster\n"
                                    "S1-MAG-1.fasta,x,1,1,1_1\nS2-MAG-1.fasta,x,1,1,1_1\nS3-MAG-1.fasta,x,1,1,1_1\n"
                                    "S1-MAG-2.fasta,x,1,2,2_0\n")
    (tables / "Wdb.csv").write_text("genome,cluster,score\nS2-MAG-1.fasta,1_1,95.1\nS1-MAG-2.fasta,2_0,82.1\n")
    script = get_package_path("smk/scripts/drep-clusters.py")
    out = subprocess.run([sys.executable, str(script), str(tmp_path / "drep")], capture_output=True, text=True, check=True)
    assert out.stdout == "S1-MAG-2\tS1-MAG-2\nS2-MAG-1\tS2-MAG-1\nS1-MAG-1\tS2-MAG-1\nS3-MAG-1\tS2-MAG-1\n"
