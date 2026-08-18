"""
Writers for the iToL annotation files `bit itol` produces
"""

SEPARATOR = "\t"

color_map = {
    "blue": "#434da7",
    "green": "#48a743",
    "red": "#c01820",
    "purple": "#512f9c",
    "black": "#000000",
    "orange": "#e07b28",
}

# FIELD_SHAPES codes, per the DATASET_BINARY template
shape_map = {
    "square": "1",
    "circle": "2",
    "star": "3",
    "rtriangle": "4",
    "ltriangle": "5",
    "check": "6",
}

BINARY_PRESENT = "1"


def resolve_color(color):
    """Accept either a named color or a literal hex code."""
    if color in color_map:
        return color_map[color]

    if color.startswith("#") and len(color) in (4, 7):
        return color

    raise ValueError(
        f'"{color}" was passed to `-c`, but that\'s neither one of the named colors '
        f"({', '.join(color_map)}) nor a hex code like \"#0000ff\".")


def read_targets(path):
    """
    Read the single-column input file into an ordered, de-duplicated list
    """
    seen = set()
    targets = []

    with open(path, "r") as target_genomes:
        for line in target_genomes:
            target = line.strip()

            if not target or target in seen:
                continue

            seen.add(target)
            targets.append(target)

    return targets


def binary_dataset(args):
    """
    DATASET_BINARY: one filled symbol per target genome
    """
    col = resolve_color(args.color)
    shape = shape_map[args.shape]

    targets = read_targets(args.wanted_genomes)

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_BINARY\nSEPARATOR TAB\n\n")
        out_file.write(f"DATASET_LABEL\t{args.dataset_label}\n\n")
        out_file.write(f"COLOR\t{col}\n\n")
        out_file.write(f"FIELD_LABELS\t{args.dataset_label}\n\n")
        out_file.write(f"FIELD_SHAPES\t{shape}\n\n")
        out_file.write(f"FIELD_COLORS\t{col}\n\n")
        out_file.write(f"HEIGHT_FACTOR\t{args.height}\n\n")

        out_file.write("DATA\n")

        for target in targets:
            out_file.write(f"{target}\t{BINARY_PRESENT}\n")


def colorstrip(args):
    """
    DATASET_COLORSTRIP: a colored band alongside each target genome
    """
    col = resolve_color(args.color)

    targets = read_targets(args.wanted_genomes)

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n\n")
        out_file.write(f"DATASET_LABEL\t{args.dataset_label}\n")
        out_file.write(f"COLOR\t{col}\n\n")

        out_file.write(f"COLOR_BRANCHES\t{'1' if args.color_branches_too else '0'}\n\n")

        out_file.write(f"STRIP_WIDTH\t{args.strip_width}\n\n")
        out_file.write("BORDER_WIDTH\t1\n")
        out_file.write("BORDER_COLOR\t#999999\n\n")

        out_file.write("DATA\n")

        for target in targets:
            out_file.write(f"{target}\t{col}\t{args.dataset_label}\n")


def style_dataset(args):
    """
    DATASET_STYLE: replaces my old `bit itol map`
    """
    col = resolve_color(args.color)

    targets = read_targets(args.wanted_genomes)

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_STYLE\nSEPARATOR TAB\n\n")
        out_file.write(f"DATASET_LABEL\t{args.dataset_label}\n\n")
        out_file.write(f"COLOR\t{col}\n\n")

        out_file.write("DATA\n")

        for target in targets:
            if args.what_to_color in ("both", "branches"):
                out_file.write(
                    f"{target}\tbranch\tnode\t{col}\t{args.line_weight}\tnormal\n")
            if args.what_to_color in ("both", "labels"):
                out_file.write(
                    f"{target}\tlabel\tnode\t{col}\t1\tbold\n")


def text_dataset(args):
    """
    DATASET_TEXT: attach a text label to each target genome
    """
    col = resolve_color(args.color)

    targets = read_targets(args.wanted_genomes)

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_TEXT\nSEPARATOR TAB\n\n")
        out_file.write(f"DATASET_LABEL\t{args.dataset_label}\n\n")
        out_file.write(f"COLOR\t{col}\n\n")

        out_file.write("DATA\n")

        for target in targets:
            out_file.write(f"{target}\t{args.text_to_add}\t-1\t{col}\tnormal\t1\t0\n")
