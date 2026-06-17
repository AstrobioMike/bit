import os
from unittest import mock
import pandas as pd # type: ignore
import pytest # type: ignore

import bit.modules.go.go as mod
from bit.modules.go.go import (
    add_percentages,
    count_terms,
    build_go_df,
    write_go_tables,
    combine_summaries,
    resolve_obo_path,
    check_go_dbs,
    get_term_info,
    summarize_annotations,
    slim_terms,
    _term_row,
    _print_related_terms,
)


# ───────────────────────── fakes ─────────────────────────

class FakeGOTerm:
    def __init__(self, namespace, depth, name, is_obsolete=False,
                 parents=None, children=None):
        self.namespace = namespace
        self.depth = depth
        self.name = name
        self.is_obsolete = is_obsolete
        self._parents = set(parents or [])
        self._children = set(children or [])

    def get_all_parents(self):
        return self._parents

    def get_all_children(self):
        return self._children


def make_fake_dag():
    # mix of namespaces, depths, one obsolete
    return {
        "GO:0008150": FakeGOTerm("biological_process", 0, "biological_process",
                                 children={"GO:0000001"}),
        "GO:0003674": FakeGOTerm("molecular_function", 0, "molecular_function"),
        "GO:0000001": FakeGOTerm("biological_process", 5, "deep BP term",
                                 parents={"GO:0008150"}),
        "GO:0009999": FakeGOTerm("biological_process", 3, "obsolete term",
                                 is_obsolete=True),
    }


class Args:
    """Minimal attribute bag for the arg-taking functions."""
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ───────────────────────── add_percentages ─────────────────────────

def test_add_percentages_basic():
    df = pd.DataFrame({"term_counts": [1, 3]})
    out = add_percentages(df)
    assert list(out["term_perc_of_annotated"]) == [25.0, 75.0]
    # original not mutated (function copies)
    assert "term_perc_of_annotated" not in df.columns


# ───────────────────────── build_go_df ─────────────────────────

def test_build_go_df_splits_obsolete_and_sorts():
    go = make_fake_dag()
    df, obsolete = build_go_df(go)

    assert obsolete == {"GO:0009999"}
    # obsolete dropped; 3 live terms remain
    assert len(df) == 3
    assert "GO:0009999" not in set(df["GO_term"])
    # when obsolete present, the is_obsolete column is dropped
    assert list(df.columns) == ["GO_term", "namespace", "depth", "name"]
    # sorted by depth ascending
    assert list(df["depth"]) == sorted(df["depth"])


def test_build_go_df_no_obsolete_keeps_is_obsolete_column():
    go = {
        "GO:1": FakeGOTerm("biological_process", 1, "a"),
        "GO:2": FakeGOTerm("molecular_function", 2, "b"),
    }
    df, obsolete = build_go_df(go)
    assert obsolete == set()
    # NOTE: with no obsolete terms the is_obsolete column is retained — see note below
    assert "is_obsolete" in df.columns


# ───────────────────────── count_terms ─────────────────────────

def test_count_terms_counts_and_skips(tmp_path):
    annots = tmp_path / "annots.tsv"
    annots.write_text(
        "# header comment\n"
        "gene1\tGO:0000001;GO:0000002\n"
        "gene2\tGO:0000001\n"
        "no_tab_line_should_skip\n"
    )
    counts = {"GO:0000001": 0, "GO:0000002": 0}
    count_terms(str(annots), counts, set(), "go_basic")
    assert counts == {"GO:0000001": 2, "GO:0000002": 1}


def test_count_terms_skips_obsolete(tmp_path):
    annots = tmp_path / "annots.tsv"
    annots.write_text("gene1\tGO:0000001;GO:0000099\n")
    counts = {"GO:0000001": 0}
    # GO:0000099 is obsolete -> skipped silently, no KeyError / no exit
    count_terms(str(annots), counts, {"GO:0000099"}, "go_basic")
    assert counts == {"GO:0000001": 1}


def test_count_terms_unknown_term_exits(tmp_path, capsys):
    annots = tmp_path / "annots.tsv"
    annots.write_text("gene1\tGO:0001234\n")
    counts = {"GO:0000001": 0}
    with pytest.raises(SystemExit) as e:
        count_terms(str(annots), counts, set(), "go_basic")
    assert e.value.code == 1
    assert "not present in the obo reference" in capsys.readouterr().out


def test_count_terms_strips_spaces(tmp_path):
    annots = tmp_path / "annots.tsv"
    annots.write_text("gene1\tGO:0000001; GO:0000002\n")   # space after ;
    counts = {"GO:0000001": 0, "GO:0000002": 0}
    count_terms(str(annots), counts, set(), "go_basic")
    assert counts["GO:0000002"] == 1


# ───────────────────────── resolve_obo_path ─────────────────────────

def test_resolve_obo_path_named_slims(monkeypatch):
    monkeypatch.setattr(mod, "check_go_dbs", lambda: "/godir")
    assert resolve_obo_path("goslim_metagenomics") == os.path.join("/godir", "goslim_metagenomics.obo")
    assert resolve_obo_path("go_basic") == os.path.join("/godir", "go-basic.obo")


def test_resolve_obo_path_passthrough(monkeypatch):
    monkeypatch.setattr(mod, "check_go_dbs", lambda: "/godir")
    assert resolve_obo_path("/my/custom.obo") == "/my/custom.obo"


# ───────────────────────── check_go_dbs ─────────────────────────

def test_check_go_dbs_returns_env_dir(monkeypatch):
    monkeypatch.setattr(mod, "get_go_data", lambda quiet=True: None)
    monkeypatch.setenv("GO_DB_DIR", "/some/go/dir")
    assert check_go_dbs() == "/some/go/dir"


# ───────────────────────── write_go_tables ─────────────────────────

def _go_df_with_counts(counts):
    df = pd.DataFrame({
        "GO_term": ["GO:1", "GO:2", "GO:3"],
        "namespace": ["biological_process", "molecular_function", "cellular_component"],
        "depth": [1, 2, 1],
        "name": ["a", "b", "c"],
        "term_counts": counts,
    })
    return add_percentages(df)


def test_write_go_tables_combined_only(tmp_path):
    df = _go_df_with_counts([5, 0, 2])
    args = Args(output_prefix=str(tmp_path / "go"), keep_zeros=False, by_namespace=False)
    write_go_tables(df, args)

    combined = tmp_path / "go.tsv"
    assert combined.is_file()
    out = pd.read_csv(combined, sep="\t")
    # zero-count row filtered out (keep_zeros False)
    assert set(out["GO_term"]) == {"GO:1", "GO:3"}


def test_write_go_tables_keep_zeros(tmp_path):
    df = _go_df_with_counts([5, 0, 2])
    args = Args(output_prefix=str(tmp_path / "go"), keep_zeros=True, by_namespace=False)
    write_go_tables(df, args)
    out = pd.read_csv(tmp_path / "go.tsv", sep="\t")
    assert len(out) == 3


def test_write_go_tables_by_namespace(tmp_path):
    df = _go_df_with_counts([5, 1, 2])
    args = Args(output_prefix=str(tmp_path / "go"), keep_zeros=False, by_namespace=True)
    write_go_tables(df, args)
    assert (tmp_path / "go-BP.tsv").is_file()
    assert (tmp_path / "go-MF.tsv").is_file()
    assert (tmp_path / "go-CC.tsv").is_file()


def test_write_go_tables_all_zero_exits(tmp_path, capsys):
    df = _go_df_with_counts([0, 0, 0])
    args = Args(output_prefix=str(tmp_path / "go"), keep_zeros=False, by_namespace=False)
    with pytest.raises(SystemExit) as e:
        write_go_tables(df, args)
    assert e.value.code == 1
    assert "no counts to any terms" in capsys.readouterr().out


# ───────────────────────── combine_summaries ─────────────────────────

def _write_summary(path, rows):
    pd.DataFrame(rows, columns=["GO_term", "namespace", "depth", "name",
                                "term_counts", "term_perc_of_annotated"]
                 ).to_csv(path, sep="\t", index=False)


def test_combine_summaries_auto_sample_names(tmp_path):
    s1 = tmp_path / "sampleA.tsv"
    s2 = tmp_path / "sampleB.tsv"
    _write_summary(s1, [["GO:1", "biological_process", 1, "a", 5, 50.0],
                        ["GO:2", "molecular_function", 2, "b", 5, 50.0]])
    _write_summary(s2, [["GO:1", "biological_process", 1, "a", 3, 75.0],
                        ["GO:3", "cellular_component", 1, "c", 1, 25.0]])

    out = tmp_path / "combined.tsv"
    args = Args(sample_names=[], input_files=[str(s1), str(s2)], output_file=str(out))
    combine_summaries(args)

    assert out.is_file()
    df = pd.read_csv(out, sep="\t")
    # sample-named count columns from filename stems
    assert "sampleA_counts" in df.columns
    assert "sampleB_counts" in df.columns
    # union of GO terms across both files, NaNs filled with 0
    assert set(df["GO_term"]) == {"GO:1", "GO:2", "GO:3"}
    assert not df.isnull().values.any()


def test_combine_summaries_name_count_mismatch_exits(tmp_path, capsys):
    s1 = tmp_path / "a.tsv"
    _write_summary(s1, [["GO:1", "biological_process", 1, "a", 1, 100.0]])
    args = Args(sample_names=["only_one", "too_many"],
                input_files=[str(s1)], output_file=str(tmp_path / "o.tsv"))
    with pytest.raises(SystemExit) as e:
        combine_summaries(args)
    assert e.value.code == 1
    assert "doesn't match the number" in capsys.readouterr().out


# ───────────────────────── _term_row / _print_related_terms ─────────────────────────

def test_term_row():
    go = make_fake_dag()
    assert _term_row(go, "GO:0000001") == ["GO:0000001", "biological_process", 5, "deep BP term"]


def test_print_related_terms_empty(capsys):
    go = make_fake_dag()
    _print_related_terms(go, set(), "Parent terms", "GO:0008150")
    assert "no parent terms for GO:0008150" in capsys.readouterr().out


def test_print_related_terms_with_terms(capsys):
    go = make_fake_dag()
    _print_related_terms(go, {"GO:0000001"}, "Child terms", "GO:0008150")
    out = capsys.readouterr().out
    assert "Child terms info:" in out
    assert "GO:0000001" in out


# ───────────────────────── get_term_info ─────────────────────────

def test_get_term_info_adds_go_prefix_and_prints(monkeypatch, capsys):
    go = make_fake_dag()
    monkeypatch.setattr(mod, "resolve_obo_path", lambda x: "/fake.obo")
    monkeypatch.setattr(mod, "load_obo", lambda p, **kw: go)

    args = Args(GO_obo_file="go_basic", GO_term="0000001", parents_only=False)
    get_term_info(args)
    out = capsys.readouterr().out
    assert "Input GO term info:" in out
    assert "Parent terms info:" in out                 # GO:0000001 has a parent
    assert "no child terms for GO:0000001" in out       # GO:0000001 has no children


def test_get_term_info_missing_term_exits(monkeypatch):
    go = make_fake_dag()
    monkeypatch.setattr(mod, "resolve_obo_path", lambda x: "/fake.obo")
    monkeypatch.setattr(mod, "load_obo", lambda p, **kw: go)

    args = Args(GO_obo_file="go_basic", GO_term="9990000", parents_only=False)
    with pytest.raises(SystemExit):
        get_term_info(args)


# ───────────────────────── summarize_annotations (end-to-end) ─────────────────────────

def test_summarize_annotations_end_to_end(tmp_path, monkeypatch):
    go = make_fake_dag()
    monkeypatch.setattr(mod, "resolve_obo_path", lambda x: "/fake.obo")
    monkeypatch.setattr(mod, "load_obo", lambda p, **kw: go)

    annots = tmp_path / "annots.tsv"
    # GO:0008150 x2, GO:0003674 x1, plus obsolete GO:0009999 (skipped)
    annots.write_text(
        "gene1\tGO:0008150;GO:0003674\n"
        "gene2\tGO:0008150;GO:0009999\n"
    )

    args = Args(GO_obo_file="go_basic", input_GO_annotations=str(annots),
                output_prefix=str(tmp_path / "summary"),
                keep_zeros=False, by_namespace=False)

    summarize_annotations(args)

    out = pd.read_csv(tmp_path / "summary.tsv", sep="\t")
    counts = dict(zip(out["GO_term"], out["term_counts"]))
    assert counts["GO:0008150"] == 2
    assert counts["GO:0003674"] == 1
    assert "GO:0009999" not in counts          # obsolete, never counted


# ───────────────────────── slim_terms ─────────────────────────

def test_slim_terms_builds_command(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "resolve_obo_path",
                        lambda x: "/resolved/" + os.path.basename(str(x)))
    out_file = tmp_path / "slim_out.tsv"
    args = Args(initial_GO_obo_file="go_basic", slimmed_GO_obo_file="goslim_metagenomics",
                input_GO_annotations="annots.tsv", mode="more",
                output_file=str(out_file))

    with mock.patch("subprocess.run") as run:
        slim_terms(args)

    assert run.called
    cmd = run.call_args[0][0]
    assert cmd[0] == "map_to_slim.py"
    assert "--association_file" in cmd and "annots.tsv" in cmd
    assert "--slim_out" in cmd and "more" in cmd
    # subprocess stdout was redirected to the opened output file
    assert run.call_args.kwargs["stdout"] is not None
