import pytest # type: ignore

from bit.modules.itol import (
    color_map,
    binary_dataset,
    colorstrip,
    itol_map,
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


def _lines(path):
    return path.read_text().splitlines()


# ───────────────────────── binary_dataset ─────────────────────────

def test_binary_dataset(targets_file, tmp_path):
    out = tmp_path / "binary.txt"
    args = Args(color="red", shape="star", dataset_label="my label",
                height=1.5, input_file=str(targets_file), output_file=str(out))
    binary_dataset(args)

    text = out.read_text()
    assert text.startswith("DATASET_BINARY\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tmy label\n" in text
    assert f"COLOR\t{color_map['red']}\n" in text
    assert "FIELD_SHAPES\t3\n" in text          # star -> "3"
    assert f"FIELD_COLORS\t{color_map['red']}\n" in text
    assert "HEIGHT_FACTOR\t1.5\n" in text

    # rows after the DATA marker are the per-target rows
    lines = _lines(out)
    data_rows = lines[lines.index("DATA") + 1:]
    assert data_rows == ["genomeA\t3", "genomeB\t3", "genomeC\t3"]


def test_binary_dataset_bad_color_raises(targets_file, tmp_path):
    args = Args(color="chartreuse", shape="square", dataset_label="x",
                height=1, input_file=str(targets_file),
                output_file=str(tmp_path / "o.txt"))
    with pytest.raises(KeyError):
        binary_dataset(args)


def test_binary_dataset_bad_shape_raises(targets_file, tmp_path):
    args = Args(color="blue", shape="hexagon", dataset_label="x",
                height=1, input_file=str(targets_file),
                output_file=str(tmp_path / "o.txt"))
    with pytest.raises(KeyError):
        binary_dataset(args)


# ───────────────────────── colorstrip ─────────────────────────

def test_colorstrip_branches_on(targets_file, tmp_path):
    out = tmp_path / "strip.txt"
    args = Args(color="green", label="mystrip", width=25,
                color_branches_too=True,
                input_file=str(targets_file), output_file=str(out))
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
    args = Args(color="green", label="s", width=25,
                color_branches_too=False,
                input_file=str(targets_file), output_file=str(out))
    colorstrip(args)
    assert "COLOR_BRANCHES\t0\n" in out.read_text()


# ───────────────────────── itol_map ─────────────────────────

def test_itol_map_both(targets_file, tmp_path):
    out = tmp_path / "map.txt"
    args = Args(color="purple", what_to_color="both", line_weight=4,
                input_file=str(targets_file), output_file=str(out))
    itol_map(args)

    text = out.read_text()
    assert text.startswith("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
    col = color_map["purple"]
    # both label and branch rows for each target
    assert f"genomeA\tlabel\t{col}\tbold" in text
    assert f"genomeA\tbranch\t{col}\tnormal\t4" in text
    assert text.count("\tlabel\t") == 3
    assert text.count("\tbranch\t") == 3


def test_itol_map_labels_only(targets_file, tmp_path):
    out = tmp_path / "map.txt"
    args = Args(color="purple", what_to_color="labels", line_weight=4,
                input_file=str(targets_file), output_file=str(out))
    itol_map(args)
    text = out.read_text()
    assert "\tlabel\t" in text
    assert "\tbranch\t" not in text


def test_itol_map_branches_only(targets_file, tmp_path):
    out = tmp_path / "map.txt"
    args = Args(color="black", what_to_color="branches", line_weight=2,
                input_file=str(targets_file), output_file=str(out))
    itol_map(args)
    text = out.read_text()
    assert "\tbranch\t" in text
    assert "\tlabel\t" not in text


# ───────────────────────── text_dataset ─────────────────────────

def test_text_dataset(targets_file, tmp_path):
    out = tmp_path / "text.txt"
    args = Args(color="blue", text_to_add="hello",
                input_file=str(targets_file), output_file=str(out))
    text_dataset(args)

    text = out.read_text()
    assert text.startswith("DATASET_TEXT\nSEPARATOR TAB\n")
    assert "DATASET_LABEL\tdata\n" in text
    col = color_map["blue"]
    # per-target: "<target>\thello\t-1\t<col>\tnormal\t1\t0"
    assert f"genomeA\thello\t-1\t{col}\tnormal\t1\t0" in text
    data_lines = [l for l in _lines(out) if "\thello\t-1\t" in l]
    assert len(data_lines) == 3