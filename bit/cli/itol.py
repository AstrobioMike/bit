import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help

def build_parser():

    desc = """
        This program helps generate various Interacitve Tree of Life (iToL) files that can be dropped onto a
        tree on the website for visualization/annotation. See itol.embl.de/help.cgi for information on the different types.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "-i",
            "--input-file",
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
            help='Color to use (default: "blue")',
            choices=["blue", "green", "red", "purple", "black"],
            default="blue"
        )

     ### subcommand cli for generating iToL binary dataset file ###

    ### subcommand cli for generating iToL binary dataset file ###
    binary_desc = """
        This subcommand creates a standard iToL binary-dataset file.
        """

    binary_parser = subparsers.add_parser(
        "binary-dataset",
        help="Create an iToL binary-dataset file",
        description=binary_desc,
        epilog="Ex. usage: `bit-itol binary-dataset -i genomes.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    binary_required = binary_parser.add_argument_group("Required Parameters")
    binary_optional = binary_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(binary_required)
    add_common_optional_arguments(binary_optional)

    binary_optional.add_argument(
        "-d",
        "--dataset-label",
        metavar="<STR>",
        help='Label of the dataset (default: "data")',
        default="data"
    )

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

    binary_parser.set_defaults(func="binary_dataset")


    ### subcommand cli for generating iToL colorstrip file ###
    colorstrip_desc = """
        This subcommand creates a standard iToL colorstrip file.
        """

    colorstrip_parser = subparsers.add_parser(
        "colorstrip",
        help="Create an iToL colorstrip file",
        description=colorstrip_desc,
        epilog="Ex. usage: `bit-itol colorstrip -i genomes.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    colorstrip_required = colorstrip_parser.add_argument_group("Required Parameters")
    colorstrip_optional = colorstrip_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(colorstrip_required)
    add_common_optional_arguments(colorstrip_optional)

    colorstrip_optional.add_argument(
        "-l",
        "--label",
        metavar="<STR>",
        help='Label used in the legend table (default: "label1")',
        default="label1"
    )

    colorstrip_optional.add_argument(
        "-w",
        "--width",
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

    colorstrip_parser.set_defaults(func="colorstrip")


    ### subcommand cli for generating iToL map file ###
    map_desc = """
        This subcommand creates a standard iToL-map file for coloring labels and/or branches.
        """

    map_parser = subparsers.add_parser(
        "map",
        help="Create an iToL map file for coloring labels and/or branches",
        description=map_desc,
        epilog="Ex. usage: `bit-itol map -i genomes.txt`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    map_required = map_parser.add_argument_group("Required Parameters")
    map_optional = map_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(map_required)
    add_common_optional_arguments(map_optional)

    map_optional.add_argument(
        "-w",
        "--what-to-color",
        help='What to color (default: "both")',
        choices=["branches", "labels", "both"],
        default="both"
    )

    map_optional.add_argument(
        "-l",
        "--line-weight",
        metavar="<NUM>",
        help='Line weight if coloring branches (default: "2")',
        default=2
    )

    add_help(map_optional)

    map_parser.set_defaults(func="itol_map")


    ### subcommand cli for generating iToL text dataset file ###
    text_desc = """
        This subcommand creates a standard iToL text-dataset file.
        """

    text_parser = subparsers.add_parser(
        "text-dataset",
        help="Create an iToL text-dataset file",
        description=text_desc,
        epilog="Ex. usage: `bit-itol text-dataset -i genomes.txt`",
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

    add_help(text_optional)

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

        if cmd in parser.subparsers.choices:
            parser.subparsers.choices[cmd].print_help(sys.stderr)
            sys.exit(0)

        print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)


    args = parser.parse_args()

    from bit.modules.general import check_files_are_found
    check_files_are_found([args.input_file])

    from bit.modules.itol import binary_dataset, colorstrip, itol_map, text_dataset

    args = preflight_checks(args)

    func_map = {
        "binary_dataset": binary_dataset,
        "colorstrip": colorstrip,
        "itol_map": itol_map,
        "text_dataset": text_dataset
    }

    func = func_map[args.func]

    func(args)


def preflight_checks(args):

    from bit.modules.general import report_message, notify_premature_exit

    if args.func == "binary_dataset":
        try:
            args.height = float(args.height)
        except ValueError:
            report_message("The value passed to `--height-factor' must be a number.")
            notify_premature_exit()

    elif args.func == "colorstrip":
        try:
            args.width = int(args.width)
        except ValueError:
            report_message("The value passed to `--width' must be an integer.")
            notify_premature_exit()

    elif args.func == "itol_map":
        try:
            args.line_weight = float(args.line_weight)
        except ValueError:
            report_message("The value passed to `--line-weight' must be a number.")
            notify_premature_exit()

    return args
