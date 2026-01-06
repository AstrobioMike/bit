import pytest
from bit.modules import input_parsing as m


class PrematureExit(Exception):
    """Raised by our monkeypatched notify_premature_exit to test failure paths."""

@pytest.fixture(autouse=True)
def patch_failure_hooks(monkeypatch):
    # capture messages sent to report_message; raise on premature exit
    messages = []

    def fake_report_message(msg, initial_indent="", subsequent_indent=""):
        messages.append(("report", msg, initial_indent, subsequent_indent))

    def fake_notify_premature_exit():
        raise PrematureExit("premature-exit")

    monkeypatch.setattr(m, "report_message", fake_report_message)
    monkeypatch.setattr(m, "notify_premature_exit", fake_notify_premature_exit)

    return messages


@pytest.mark.parametrize("fname,expected", [
    ("sample_R1_.fastq.gz", ("sample", "R1")),
    ("sample_R2_.fastq.gz", ("sample", "R2")),
    ("sample-R1-.fq.gz",    ("sample", "R1")),
    ("sample-R2-.fq.gz",    ("sample", "R2")),
    ("sample.R1..fastq.gz", ("sample", "R1")),
    ("sample.R2..fastq.gz", ("sample", "R2")),
    ("sample_1..fastq.gz",  ("sample", "R1")),
    ("sample_2..fastq.gz",  ("sample", "R2")),
    ("no_tag.fastq.gz",     (None, None)),
])
def test_parse_read_filename_variants(fname, expected):
    assert m.parse_read_filename(fname) == expected


def test_validate_extension_accepts_known_extensions(tmp_path):
    for ext in m.accepted_read_extensions:
        p = tmp_path / f"x{ext}"
        p.write_text("data")
        m.validate_extension(p)

def test_validate_extension_rejects_unknown_extension(tmp_path):
    p = tmp_path / "x.fq"  # not gzipped → reject
    p.write_text("data")
    with pytest.raises(PrematureExit):
        m.validate_extension(p)


def test_get_input_reads_dict_from_paths_happy_path(tmp_path):
    r1 = tmp_path / "samp_R1_.fastq.gz"
    r2 = tmp_path / "samp_R2_.fastq.gz"
    r1.write_text("r1")
    r2.write_text("r2")

    out = m.get_input_reads_dict_from_paths(r1, r2)
    assert list(out.keys()) == ["samp"]
    assert set(out["samp"].keys()) == {"R1", "R2"}
    assert out["samp"]["R1"] == str(r1.resolve())
    assert out["samp"]["R2"] == str(r2.resolve())


def test_get_input_reads_dict_from_paths_missing_designation_calls_exit(tmp_path):
    bad = tmp_path / "samp.fastq.gz"  # no R1/R2 tag
    bad.write_text("x")
    with pytest.raises(PrematureExit):
        m.get_input_reads_dict_from_paths(bad, None)


def test_get_input_reads_dict_from_paths_wrong_slot_calls_exit(tmp_path):
    # file says R2 but provided as R1 argument
    wrong = tmp_path / "samp_R2_.fastq.gz"
    wrong.write_text("x")
    with pytest.raises(PrematureExit):
        m.get_input_reads_dict_from_paths(wrong, None)


def test_get_input_reads_dict_from_dir_pairs_samples_and_ignores_noise(tmp_path):
    # A complete pair
    (tmp_path / "A_R1_.fastq.gz").write_text("a1")
    (tmp_path / "A_R2_.fastq.gz").write_text("a2")
    # B has extra unrelated files that should be ignored
    (tmp_path / "notes.txt").write_text("ignore me")
    (tmp_path / "weird.fq").write_text("ignore me")  # bad ext
    (tmp_path / "junk.fastq.gz").write_text("no R tag")  # will be parsed as None,None and skipped

    out = m.get_input_reads_dict_from_dir(tmp_path)
    assert list(out.keys()) == ["A"]
    assert set(out["A"].keys()) == {"R1", "R2"}


def test_get_input_reads_dict_from_dir_detects_incomplete_pairs_and_exits(tmp_path, patch_failure_hooks):
    # Only R1 present for B → should trigger error
    (tmp_path / "B_R1_.fastq.gz").write_text("b1")
    with pytest.raises(PrematureExit):
        m.get_input_reads_dict_from_dir(tmp_path)

    # Ensure a diagnostic report_message was sent
    reports = [x for x in patch_failure_hooks if x[0] == "report"]
    assert reports, "expected report_message() to be called"
    # The message mentions the input directory
    assert str(tmp_path) in reports[-1][1]


def test_get_input_reads_dict_from_dir_handles_multiple_samples(tmp_path):
    # C paired
    (tmp_path / "C_R1_.fastq.gz").write_text("c1")
    (tmp_path / "C_R2_.fastq.gz").write_text("c2")
    # D paired with different accepted tags
    (tmp_path / "D-R1-.fq.gz").write_text("d1")
    (tmp_path / "D-R2-.fq.gz").write_text("d2")

    out = m.get_input_reads_dict_from_dir(tmp_path)
    assert set(out.keys()) == {"C", "D"}
    assert set(out["C"].keys()) == {"R1", "R2"}
    assert set(out["D"].keys()) == {"R1", "R2"}


def test_get_input_reads_dict_from_dir_skips_files_with_no_designation(tmp_path):
    (tmp_path / "E.fastq.gz").write_text("no tag")
    out = m.get_input_reads_dict_from_dir(tmp_path)
    assert out == {}


def test_parse_read_filename_uses_basename_not_path(tmp_path):
    # Ensure directories in the path don’t confuse parsing
    p = tmp_path / "subdir"
    p.mkdir()
    f = p / "sample.R1..fastq.gz"
    f.write_text("x")
    # give full path; function uses Path(...).name internally
    samp, which = m.parse_read_filename(str(f))
    assert samp == "sample" and which == "R1"
