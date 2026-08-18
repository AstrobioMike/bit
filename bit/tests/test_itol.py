import pytest # type: ignore

from bit.modules.itol import (
    color_map,
    resolve_color,
    shape_map,
    read_targets,
    binary_dataset,
    colorstrip,
    style_dataset,
    text_dataset,
)


class Args:
    def __init__(self, **kw):
        self.__dict__.update(kw)


@pytest.fixture
def targets_file(tmp_path):
    f = tmp_path / "targets.txt"
    # blank/whitespace lines exercise the .strip() handling
    f.write_text("genomeA\ngenomeB\n  genomeC  \n")
    return f


@pytest.fixture
def messy_targets_file(tmp_path):
    """
    The shape real hand-made lists tend to arrive in: a blank line in the middle,
    a repeat, and trailing whitespace.
    """
    f = tmp_path / "messy.txt"
    f.write_text("genomeA\n\ngenomeB\ngenomeA\n  genomeC  \n\n")
    return f


def _lines(path):
    return path.read_text().splitlines()


def _data_rows(path):
    lines = _lines(path)
    return [l for l in lines[lines.index("DATA") + 1:] if l.strip()]


# ───────────────────────── read_targets ─────────────────────────

def test_read_targets_strips_blanks_and_duplicates(messy_targets_file):
    """
    A blank line used to become a row with an empty first field, which iToL drops
    silently along with the annotation on it; a repeated ID used to become a
    duplicate annotation row.
    """
    assert read_targets(str(messy_targets_file)) == ["genomeA", "genomeB", "genomeC"]


def test_read_targets_preserves_input_order(tmp_path):
    f = tmp_path / "ordered.txt"
    f.write_text("zeta\nalpha\nmu\n")
    assert read_targets(str(f)) == ["zeta", "alpha", "mu"]


# ───────────────────────── binary_dataset ─────────────────────────

def test_binary_dataset(targets_file, tmp_path):
    out = tmp_path / "binary.txt"
    args = Args(color="red", shape="star", dataset_label="my label",
                height=1.5, wanted_genomes=str(targets_file), output_file=str(out))
    binary_dataset(args)

    text = out.read_text()
    assert text.startswith("DATASET_BINARY\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tmy label\n" in text
    assert f"COLOR\t{color_map['red']}\n" in text
    assert "FIELD_SHAPES\t3\n" in text          # star -> "3"
    assert f"FIELD_COLORS\t{color_map['red']}\n" in text
    assert "HEIGHT_FACTOR\t1.5\n" in text


@pytest.mark.parametrize("shape", sorted(shape_map))
def test_binary_dataset_data_values_are_presence_flags_not_shape_codes(
        shape, targets_file, tmp_path):
    """
    The DATASET_BINARY template allows only 1, 0 or -1 under DATA -- the symbol is
    chosen once via FIELD_SHAPES. Writing the shape code per row rendered fine for
    `square` (code 1, which is also a valid presence flag) and silently produced
    nothing for every other shape, so this is parametrized over all of them.
    """
    out = tmp_path / f"binary-{shape}.txt"
    args = Args(color="blue", shape=shape, dataset_label="d", height=1,
                wanted_genomes=str(targets_file), output_file=str(out))
    binary_dataset(args)

    assert f"FIELD_SHAPES\t{shape_map[shape]}\n" in out.read_text()
    assert _data_rows(out) == ["genomeA\t1", "genomeB\t1", "genomeC\t1"]


def test_binary_dataset_bad_color_raises(targets_file, tmp_path):
    args = Args(color="chartreuse", shape="square", dataset_label="x",
                height=1, wanted_genomes=str(targets_file),
                output_file=str(tmp_path / "o.txt"))
    with pytest.raises(ValueError):
        binary_dataset(args)


def test_binary_dataset_bad_shape_raises(targets_file, tmp_path):
    args = Args(color="blue", shape="hexagon", dataset_label="x",
                height=1, wanted_genomes=str(targets_file),
                output_file=str(tmp_path / "o.txt"))
    with pytest.raises(KeyError):
        binary_dataset(args)


# ───────────────────────── colorstrip ─────────────────────────

def test_colorstrip_branches_on(targets_file, tmp_path):
    out = tmp_path / "strip.txt"
    args = Args(color="green", dataset_label="mystrip", strip_width=25,
                color_branches_too=True,
                wanted_genomes=str(targets_file), output_file=str(out))
    colorstrip(args)

    text = out.read_text()
    assert text.startswith("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
    assert "COLOR_BRANCHES\t1\n" in text
    assert "STRIP_WIDTH\t25\n" in text
    assert "BORDER_COLOR\t#999999\n" in text
    # each DATA row: "<target>\t<color>\t<label>"
    assert f"genomeA\t{color_map['green']}\tmystrip" in text
    assert text.count(f"\t{color_map['green']}\tmystrip") == 3


def test_colorstrip_branches_off(targets_file, tmp_path):
    out = tmp_path / "strip.txt"
    args = Args(color="green", dataset_label="s", strip_width=25,
                color_branches_too=False,
                wanted_genomes=str(targets_file), output_file=str(out))
    colorstrip(args)
    assert "COLOR_BRANCHES\t0\n" in out.read_text()


# ───────────────────────── style_dataset ─────────────────────────

def _style_args(targets_file, out, **over):
    base = dict(color="blue", dataset_label="data", what_to_color="branches",
                apply_to="node", line_weight=3,
                wanted_genomes=str(targets_file), output_file=str(out))
    base.update(over)
    return Args(**base)


def test_style_dataset_branches(targets_file, tmp_path):
    out = tmp_path / "style.txt"
    style_dataset(_style_args(targets_file, out, dataset_label="my gene"))

    text = out.read_text()
    assert text.startswith("DATASET_STYLE\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tmy gene\n" in text
    col = color_map["blue"]
    assert _data_rows(out) == [f"{g}\tbranch\tnode\t{col}\t3\tnormal"
                               for g in ("genomeA", "genomeB", "genomeC")]


def test_style_dataset_both_emits_branch_and_label_rows(targets_file, tmp_path):
    out = tmp_path / "style.txt"
    style_dataset(_style_args(targets_file, out, what_to_color="both"))

    rows = _data_rows(out)
    assert len(rows) == 6
    assert sum("\tbranch\t" in r for r in rows) == 3
    assert sum("\tlabel\t" in r for r in rows) == 3


def test_style_dataset_labels_only(targets_file, tmp_path):
    out = tmp_path / "style.txt"
    style_dataset(_style_args(targets_file, out, what_to_color="labels"))

    rows = _data_rows(out)
    assert all("\tlabel\t" in r for r in rows)
    assert not any("\tbranch\t" in r for r in rows)


def test_style_dataset_applies_to_clade(targets_file, tmp_path):
    """`clade` propagates the style to every descendant, `node` doesn't."""
    out = tmp_path / "style.txt"
    style_dataset(_style_args(targets_file, out, apply_to="clade"))
    assert all(r.split("\t")[2] == "clade" for r in _data_rows(out))


def test_style_dataset_bad_color_raises(targets_file, tmp_path):
    with pytest.raises(ValueError):
        style_dataset(_style_args(targets_file, tmp_path / "o.txt",
                                  color="chartreuse"))


# ───────────────────────── text_dataset ─────────────────────────

def test_text_dataset(targets_file, tmp_path):
    out = tmp_path / "text.txt"
    args = Args(color="blue", text_to_add="hello", dataset_label="data",
                wanted_genomes=str(targets_file), output_file=str(out))
    text_dataset(args)

    text = out.read_text()
    assert text.startswith("DATASET_TEXT\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tdata\n" in text
    col = color_map["blue"]
    # per-target: "<target>\thello\t-1\t<col>\tnormal\t1\t0"
    assert f"genomeA\thello\t-1\t{col}\tnormal\t1\t0" in text
    data_lines = [l for l in _lines(out) if "\thello\t-1\t" in l]
    assert len(data_lines) == 3


def test_text_dataset_honors_dataset_label(targets_file, tmp_path):
    out = tmp_path / "text.txt"
    args = Args(color="blue", text_to_add="hi", dataset_label="my notes",
                wanted_genomes=str(targets_file), output_file=str(out))
    text_dataset(args)
    assert "DATASET_LABEL\tmy notes\n" in out.read_text()



# ───────────────────────── color resolution ─────────────────────────

@pytest.mark.parametrize("name", sorted(color_map))
def test_resolve_color_accepts_every_named_color(name):
    assert resolve_color(name) == color_map[name]


@pytest.mark.parametrize("hex_code", ["#0000ff", "#abc"])
def test_resolve_color_accepts_hex(hex_code):
    assert resolve_color(hex_code) == hex_code


@pytest.mark.parametrize("bad", ["chartreuse", "0000ff", "#00ff"])
def test_resolve_color_rejects_anything_else(bad):
    with pytest.raises(ValueError):
        resolve_color(bad)
