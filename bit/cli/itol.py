import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg

def build_parser(parent_subparsers=None):

    desc = """
        This program helps generate various Interacitve Tree of Life (iToL) files that can be dropped onto a
        tree on the website for visualization/annotation. See itol.embl.de/help.cgi for information on the different types.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "itol",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    add_help(parser)

    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "-g",
            "--wanted-genomes",
            help = "single-column file with genomes to decorate (must match the IDs in the tree file)",
            metavar = "<FILE>",
            required = True
        )

    def add_common_optional_arguments(group):
        group.add_argument(
            "-o",
            "--output-file",
            help = 'Name of output file for iToL (default: "itol.txt")',
            metavar = "<FILE>",
            default = "itol.txt"
        )

        group.add_argument(
            "-c",
            "--color",
            metavar="<STR>",
            help=('Color to use -- either a name (blue, green, red, purple, black, orange) '
                  'or a hex code like "#0000ff" (default: "blue")'),
            default="blue"
        )

    # `DATASET_LABEL` in the written file -- the dataset's name in iToL's Datasets
    # panel. Not added to `map`: TREE_COLORS has no such field.
    def add_dataset_label_argument(group):
        group.add_argument(
            "-d",
            "--dataset-label",
            metavar="<STR>",
            help='Label of the dataset (default: "data")',
            default="data"
        )

    ### subcommand cli for generating an iToL binary-dataset file ###
    binary_desc = """
        This subcommand creates a standard iToL binary-dataset file.
        """

    binary_parser = subparsers.add_parser(
        "binary-dataset",
        help="Create an iToL binary-dataset file",
        description=binary_desc,
        epilog="Ex. usage: `bit itol binary-dataset -g genomes.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    binary_required = binary_parser.add_argument_group("Required Parameters")
    binary_optional = binary_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(binary_required)
    add_common_optional_arguments(binary_optional)

    add_dataset_label_argument(binary_optional)

    binary_optional.add_argument(
        "-s",
        "--shape",
        help='Shape to add (default: "square)',
        choices=["square", "circle", "star", "rtriangle", "ltriangle", "check"],
        default="square"
    )

    binary_optional.add_argument(
        "-H",
        "--height-factor",
        metavar="<NUM>",
        help='Increase or decrease symbol size. Values below 1 will decrease the standard size, above 1 will increase it (default: "1")',
        default="1",
        dest="height"
    )

    add_help(binary_optional)

    add_version_arg(binary_optional)

    binary_parser.set_defaults(func="binary_dataset")


    ### subcommand cli for generating an iToL colorstrip file ###
    colorstrip_desc = """
        This subcommand creates a standard iToL colorstrip file.
        """

    colorstrip_parser = subparsers.add_parser(
        "colorstrip",
        help="Create an iToL colorstrip file",
        description=colorstrip_desc,
        epilog="Ex. usage: `bit itol colorstrip -g genomes.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    colorstrip_required = colorstrip_parser.add_argument_group("Required Parameters")
    colorstrip_optional = colorstrip_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(colorstrip_required)
    add_common_optional_arguments(colorstrip_optional)

    add_dataset_label_argument(colorstrip_optional)

    colorstrip_optional.add_argument(
        "-W",
        "--strip-width",
        metavar="<NUM>",
        help='Width of the colorstrip (default: "25")',
        default="25"
    )

    colorstrip_optional.add_argument(
        "--color-branches-too",
        help="Add this flag if wanting to color branches also",
        action="store_true"
    )

    add_help(colorstrip_optional)

    add_version_arg(colorstrip_optional)

    colorstrip_parser.set_defaults(func="colorstrip")


    ### subcommand cli for generating an iToL-style dataset file ###
    style_desc = """
        This subcommand creates a standard iToL style-dataset file for coloring labels and/or branches.
        """

    style_parser = subparsers.add_parser(
        "style",
        help="Create an iToL style-dataset file for coloring labels and/or branches",
        description=style_desc,
        epilog="Ex. usage: `bit itol style -g genomes.txt -d my-gene`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    style_required = style_parser.add_argument_group("Required Parameters")
    style_optional = style_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(style_required)
    add_common_optional_arguments(style_optional)

    add_dataset_label_argument(style_optional)

    style_optional.add_argument(
        "--what-to-color",
        help='What to color (default: "branches")',
        choices=["branches", "labels", "both"],
        default="branches"
    )

    style_optional.add_argument(
        "-a",
        "--apply-to",
        help=('Whether to style just the matched node, or the whole clade beneath it '
              '(default: "node")'),
        choices=["node", "clade"],
        default="node",
        dest="apply_to"
    )

    style_optional.add_argument(
        "-l",
        "--line-weight",
        metavar="<NUM>",
        help='Line weight if coloring branches (default: "3")',
        default=3
    )

    add_help(style_optional)

    add_version_arg(style_optional)

    style_parser.set_defaults(func="style_dataset")


    ### subcommand cli for generating an iToL text dataset file ###
    text_desc = """
        This subcommand creates a standard iToL text-dataset file.
        """

    text_parser = subparsers.add_parser(
        "text-dataset",
        help="Create an iToL text-dataset file",
        description=text_desc,
        epilog="Ex. usage: `bit itol text-dataset -g genomes.txt -t 'my note'`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    text_required = text_parser.add_argument_group("Required Parameters")
    text_optional = text_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(text_required)

    text_required.add_argument(
        "-t",
        "--text-to-add",
        metavar="<STR>",
        help='Text to add to the target genomes',
        required=True,
    )

    add_common_optional_arguments(text_optional)

    add_dataset_label_argument(text_optional)

    add_help(text_optional)

    add_version_arg(text_optional)

    text_parser.set_defaults(func="text_dataset")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so approriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in ("-v", "--version"):
            from bit.modules.general import report_version
            report_version()
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            parser.subparsers.choices[cmd].print_help(sys.stderr)
            sys.exit(0)

        print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)


    args = parser.parse_args()

    from bit.modules.general import check_files_are_found
    check_files_are_found([args.wanted_genomes])

    from bit.modules.itol import (binary_dataset, colorstrip, style_dataset,
                                  text_dataset)

    args = preflight_checks(args)

    func_map = {
        "binary_dataset": binary_dataset,
        "colorstrip": colorstrip,
        "style_dataset": style_dataset,
        "text_dataset": text_dataset
    }

    func = func_map[args.func]

    func(args)


def preflight_checks(args):

    from bit.modules.general import report_message, notify_premature_exit

    from bit.modules.itol import resolve_color

    try:
        resolve_color(args.color)
    except ValueError as e:
        report_message(str(e))
        notify_premature_exit()

    if args.func == "binary_dataset":
        try:
            args.height = float(args.height)
        except ValueError:
            report_message("The value passed to `--height-factor' must be a number.")
            notify_premature_exit()

    elif args.func == "colorstrip":
        try:
            args.strip_width = int(args.strip_width)
        except ValueError:
            report_message("The value passed to `--strip-width' must be an integer.")
            notify_premature_exit()

    elif args.func == "style_dataset":
        try:
            args.line_weight = float(args.line_weight)
        except ValueError:
            report_message("The value passed to `--line-weight' must be a number.")
            notify_premature_exit()

    return args
