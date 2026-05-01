color_map = {
    "blue": "#434da7",
    "green": "#48a743",
    "red": "#c01820",
    "purple": "#512f9c",
    "black": "#000000",
}

def binary_dataset(args):

    col = color_map[args.color]

    shape_map = {
        "square": "1",
        "circle": "2",
        "star": "3",
        "rtriangle": "4",
        "ltriangle": "5",
        "check": "6"
    }

    shape = shape_map[args.shape]

    target_list = []

    with open(args.input_file, "r") as target_genomes:
        for genome in target_genomes:
            target_list.append(genome.strip())

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_BINARY\nSEPARATOR TAB\n\n")
        out_file.write("DATASET_LABEL" + "\t" + str(args.dataset_label) + "\n\n")
        out_file.write("COLOR\t" + str(col) + "\n\n")
        out_file.write("FIELD_LABELS\tf1\n\n")
        out_file.write("FIELD_SHAPES\t" + str(shape) + "\n\n")
        out_file.write("FIELD_COLORS\t" + str(col) + "\n\n")
        out_file.write("HEIGHT_FACTOR\t" + str(args.height) + "\n\n")

        out_file.write("DATA\n")

        for target in target_list:
            out_file.write(str(target) + "\t" + str(shape) + "\n")


def colorstrip(args):

    col = color_map[args.color]

    target_list = []

    with open(args.input_file, "r") as target_genomes:
        for genome in target_genomes:
            target_list.append(genome.strip())

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_COLORSTRIP" + "\n" + "SEPARATOR TAB" + "\n\n" + "DATASET_LABEL" + "\t" + str(args.label) + "\n" + "COLOR" + "\t" + str(col) + "\n\n")

        if args.color_branches_too:
            out_file.write("COLOR_BRANCHES\t1\n\n")
        else:
            out_file.write("COLOR_BRANCHES\t0\n\n")

        out_file.write("STRIP_WIDTH" + "\t" + str(args.width) + "\n\n")
        out_file.write("BORDER_WIDTH" + "\t" + "1" + "\n")
        out_file.write("BORDER_COLOR" + "\t" + "#999999" + "\n\n")
        out_file.write("DATA\n\n")

        for target in target_list:
            out_file.write(str(target) + "\t" + str(col) + "\t" + str(args.label) + "\n")


def itol_map(args):

    col = color_map[args.color]

    target_list = []

    with open(args.input_file, "r") as target_genomes:
        for genome in target_genomes:
            target_list.append(genome.strip())

    with open(args.output_file, "w") as out_file:

        out_file.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n\n")

        if args.what_to_color in ["both", "labels"]:
            for target in target_list:
                out_file.write(str(target) + "\tlabel\t" + str(col) + "\tbold\n")

        if args.what_to_color in ["both", "branches"]:
            for target in target_list:
                out_file.write(str(target) + "\tbranch\t" + str(col) + "\tnormal\t" + str(args.line_weight) + "\n")


def text_dataset(args):

    col = color_map[args.color]

    target_list = []

    with open(args.input_file, "r") as target_genomes:
        for genome in target_genomes:
            target_list.append(genome.strip())

    with open(args.output_file, "w") as out_file:

        out_file.write("DATASET_TEXT\nSEPARATOR TAB\n\n")

        # setting DATASET_LABEL
        out_file.write("DATASET_LABEL\tdata\n\n")

        # setting dataset main color
        out_file.write("COLOR\t" + str(col) + "\n\n")

        # writing lines for each labels
        out_file.write("DATA\n")

        for target in target_list:
            out_file.write(str(target) + "\t" + str(args.text_to_add) + "\t" + "-1" + "\t" + str(col) + "\t" + "normal" + "\t" + "1" + "\t" + "0" + "\n")

