import gzip
import random

import pytest # type: ignore

from bit._native import scan_foldbacks
from bit.modules.filter_foldbacks import classify_pair, filter_foldbacks, output_paths
from bit.tests.utils import run_cli


_COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_COMP)[::-1]


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def template():
    return rand_seq(3000, random.Random(7))


def pair_from(fragment, name):
    r1, r2 = fragment[:150], rc(fragment[-150:])
    return (f"{name}/1", r1, "I" * len(r1)), (f"{name}/2", r2, "I" * len(r2))


def write_fastq(path, recs):
    with open(path, "w") as f:
        for name, seq, qual in recs:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_names(path):
    with gzip.open(path, "rt") as f:
        return [l[1:].strip() for i, l in enumerate(f) if i % 4 == 0]


def read_seqs(path):
    with gzip.open(path, "rt") as f:
        return [l.strip() for i, l in enumerate(f) if i % 4 == 1]


def test_normal_pair_has_no_fold(template):
    r1, r2 = pair_from(template[100:500], "n")
    scan = scan_foldbacks([r1[1]], [r2[1]], 17, 20)[0]
    assert classify_pair(scan) == (False, False)


def test_short_insert_overlap_is_not_a_fold(template):
    # mates overlapping each other (a 200 bp fragment) share sequence only reverse complemented
    r1, r2 = pair_from(template[100:300], "s")
    scan = scan_foldbacks([r1[1]], [r2[1]], 17, 20)[0]
    assert scan[6] == 0 and scan[7] > 20
    assert classify_pair(scan) == (False, False)


def test_fold_between_mates(template):
    # 300 bp of the template, then back over 250 bp of it on the other strand
    frag = template[1000:1300] + rc(template[1050:1300])
    r1, r2 = pair_from(frag, "b")
    scan = scan_foldbacks([r1[1]], [r2[1]], 17, 20)[0]
    in_read, between = classify_pair(scan)
    assert between and not in_read


def test_fold_in_read_cut_to_longer_arm(template):
    # read 1 runs 110 bp into the template then turns back 30 bp short of where it turned
    frag = template[500:610] + rc(template[530:580]) + template[2000:2300]
    r1, r2 = pair_from(frag, "r")
    s1, e1, f1, *_ = scan_foldbacks([r1[1]], [r2[1]], 17, 20)[0]
    assert f1 == 1 and s1 == 0 and 75 <= e1 <= 85     # the first arm, to where its pairing ends


def test_single_end_input(template):
    reads = [template[100:250], template[500:610] + rc(template[530:570])]
    out = scan_foldbacks(reads, None, 17, 20)
    assert out[0][2] == 0 and out[1][2] == 1
    assert all(row[3:] == (0, 0, 0, 0, 0) for row in out)


def test_mismatched_batches_rejected():
    with pytest.raises(ValueError):
        scan_foldbacks(["ACGT" * 40], ["ACGT" * 40, "ACGT" * 40], 17, 20)


def test_filter_foldbacks_writes_pairs_and_singles(tmp_path, template):
    pairs = [pair_from(template[i:i + 400], f"n{i}") for i in range(0, 2400, 300)]
    pairs.append(pair_from(template[1000:1300] + rc(template[1050:1300]), "fold"))
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    write_fastq(r1, [p[0] for p in pairs])
    write_fastq(r2, [p[1] for p in pairs])
    out = tmp_path / "out"
    out.mkdir()
    counts, paths = filter_foldbacks(str(r1), str(r2), str(out), show_progress=False)
    assert counts.reads_or_pairs_in == len(pairs)
    assert counts.foldback_pairs == 1 and counts.pairs_out == len(pairs) - 1
    assert "fold/1" not in read_names(paths["r1"]) and len(read_names(paths["r1"])) == len(pairs) - 1
    assert read_names(paths["singles"]) == ["fold/1", "fold/2"]
    assert open(paths["summary"]).readline() == "measure\tvalue\n"


def test_filter_foldbacks_single_end_cuts_reads(tmp_path, template):
    folded = template[500:610] + rc(template[530:570])
    write_fastq(tmp_path / "r.fq", [("a", template[:150], "I" * 150), ("b", folded, "I" * len(folded))])
    out = tmp_path / "out"
    out.mkdir()
    counts, paths = filter_foldbacks(str(tmp_path / "r.fq"), None, str(out), show_progress=False)
    seqs = read_seqs(paths["reads"])
    assert counts.reads_cut_at_a_fold == 1 and seqs[0] == template[:150]
    assert 70 <= len(seqs[1]) <= 85 and folded.startswith(seqs[1])


def test_output_paths_prefix():
    p = output_paths("d", "s1-", paired=True, fastq=False)
    assert p["r1"].endswith("d/s1-filtered_R1.fasta.gz") and p["singles"].endswith("s1-foldback-mates.fasta.gz")


def test_cli(tmp_path, template):
    pairs = [pair_from(template[i:i + 400], f"n{i}") for i in range(0, 2400, 300)]
    write_fastq(tmp_path / "r1.fq", [p[0] for p in pairs])
    write_fastq(tmp_path / "r2.fq", [p[1] for p in pairs])
    out = tmp_path / "cli"
    result = run_cli(["bit", "filter-foldbacks", "-i", str(tmp_path / "r1.fq"),
                      "-I", str(tmp_path / "r2.fq"), "-o", str(out), "-j", "1"])
    assert "pairs scanned" in result.stdout
    # no fold-back pairs: no mates file, and no line pointing to one
    assert "Mates of fold-back pairs" not in result.stdout
    assert sorted(p.name for p in out.iterdir()) == ["command-execution-info.txt",
                                                     "filtered_R1.fastq.gz", "filtered_R2.fastq.gz",
                                                     "summary.tsv"]


def test_no_foldback_mates_file_without_foldbacks(tmp_path, template):
    pairs = [pair_from(template[i:i + 400], f"n{i}") for i in range(0, 2400, 300)]
    write_fastq(tmp_path / "r1.fq", [p[0] for p in pairs])
    write_fastq(tmp_path / "r2.fq", [p[1] for p in pairs])
    out = tmp_path / "out"
    out.mkdir()
    counts, paths = filter_foldbacks(str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq"), str(out),
                                     show_progress=False)
    assert counts.foldback_pairs == 0 and paths["singles"] is None
    assert not (out / "foldback-mates.fastq.gz").exists()


# long reads (--hifi)

@pytest.fixture
def long_template():
    return rand_seq(20000, random.Random(11))


def test_trim_long_folds_cuts_to_longer_arm(long_template):
    from bit._native import trim_long_folds
    t = long_template
    # 3 kb out to 3,000, back from 2,800 (200 bp unpaired middle) to 1,500
    folded = t[:3000] + rc(t[1500:2800])
    plain = t[5000:17000]
    (s1, e1, f1, _, m1), (s2, e2, f2, _, m2) = trim_long_folds([folded, plain])
    assert (s1, f1) == (0, 1) and 2750 <= e1 <= 2800 and 150 <= m1 <= 250
    assert (s2, e2, f2, m2) == (0, len(plain), 0, -1)


def test_trim_long_folds_leaves_inverted_repeat_in_long_flanks(long_template):
    from bit._native import trim_long_folds
    t = long_template
    # a genomic inverted repeat: 300 bp arm, 150 bp spacer, the arm reverse complemented, in 3 kb flanks
    read = t[:3450] + rc(t[3000:3300]) + t[6000:9000]
    assert trim_long_folds([read])[0][2] == 0


def test_filter_foldbacks_hifi(tmp_path, long_template):
    t = long_template
    folded = t[:3000] + rc(t[1500:2800])
    short_left = t[10000:10120] + rc(t[9000:10080])   # longer arm on the right, kept
    recs = [("plain", t[4000:9000]), ("fold", folded), ("right", short_left)]
    write_fastq(tmp_path / "hifi.fq", [(n, s, "I" * len(s)) for n, s in recs])
    out = tmp_path / "out"
    out.mkdir()
    counts, paths = filter_foldbacks(str(tmp_path / "hifi.fq"), None, str(out), show_progress=False,
                                     hifi=True)
    seqs = dict(zip(read_names(paths["reads"]), read_seqs(paths["reads"])))
    assert counts.reads_or_pairs_in == 3 and counts.reads_cut_at_a_fold == 2
    assert seqs["plain"] == t[4000:9000]
    assert folded.startswith(seqs["fold"]) and 2750 <= len(seqs["fold"]) <= 2800
    assert short_left.endswith(seqs["right"]) and 1000 <= len(seqs["right"]) <= 1080
    summary = dict(l.rstrip("\n").split("\t") for l in open(paths["summary"]))
    assert summary["folds_cut"] == "2" and "unpaired_middle_median" in summary


def test_filter_foldbacks_hifi_drops_short_pieces(tmp_path, long_template):
    t = long_template
    # the kept arm (~950 bp) is under min_length, so the read goes
    folded = t[:1000] + rc(t[700:950])
    write_fastq(tmp_path / "hifi.fq", [("fold", folded, "I" * len(folded))])
    out = tmp_path / "out"
    out.mkdir()
    counts, paths = filter_foldbacks(str(tmp_path / "hifi.fq"), None, str(out), show_progress=False,
                                     hifi=True, min_length=2000)
    assert counts.reads_dropped_short == 1 and counts.reads_out == 0


def test_filter_foldbacks_hifi_refuses_pairs(tmp_path):
    with pytest.raises(ValueError):
        filter_foldbacks("r1.fq", "r2.fq", str(tmp_path), show_progress=False, hifi=True)


def test_cli_hifi(tmp_path, long_template):
    t = long_template
    folded = t[:3000] + rc(t[1500:2800])
    write_fastq(tmp_path / "hifi.fq", [("fold", folded, "I" * len(folded))])
    out = tmp_path / "cli"
    result = run_cli(["bit", "filter-foldbacks", "--hifi", "-i", str(tmp_path / "hifi.fq"),
                      "-o", str(out), "-j", "1"])
    assert "reads cut at a fold" in result.stdout and "unpaired middle" in result.stdout
    assert sorted(p.name for p in out.iterdir()) == ["command-execution-info.txt",
                                                     "filtered.fastq.gz", "summary.tsv"]


def test_cli_rejects_bad_kmer_size(tmp_path, template):
    import subprocess
    write_fastq(tmp_path / "r.fq", [("a", template[:150], "I" * 150)])
    result = subprocess.run(["bit", "filter-foldbacks", "-i", str(tmp_path / "r.fq"), "-o",
                             str(tmp_path / "o"), "-k", "4"], capture_output=True, text=True)
    assert result.returncode != 0 and "kmer-size" in result.stderr


def test_short_read_defaults():
    from bit.modules.filter_foldbacks import KMER_SIZE, MIN_SHARED
    # 16 shared 13-mers: a 28 bp stretch (see the module's notes)
    assert (KMER_SIZE, MIN_SHARED) == (13, 16)


def test_cli_detailed_help_shows_fine_tuning():
    import subprocess
    plain = subprocess.run(["bit", "filter-foldbacks", "-h"], capture_output=True, text=True)
    detailed = subprocess.run(["bit", "filter-foldbacks", "-H"], capture_output=True, text=True)
    out_plain = plain.stdout + plain.stderr
    out_detailed = detailed.stdout + detailed.stderr
    assert "--kmer-size" not in out_plain and "HOW IT WORKS" not in out_plain
    assert "--kmer-size" in out_detailed and "HOW IT WORKS" in out_detailed
