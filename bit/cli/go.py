import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg


def build_parser():

    desc = """
        This program has helpers for working Gene Ontology (GO) annotations. See subcommand-specific
        help menus for more info.
        """

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
            "-i",
            "--input-GO-annotations",
            help="Input GO annotations file (2-column, tab-delimited, where the first column holds gene IDs, and the second column holds GO terms delimited with a semi-colon; no header).",
            metavar="<FILE>",
            required=True
        )

    def add_common_optional_arguments(group):

        group.add_argument(
            "-g",
            "--GO-obo-file",
            metavar="<FILE>",
            help='GO obo file to use (e.g. from: geneontology.org/docs/download-ontology/). By default will \
                use "go-basic.obo". "goslim_metagenomics.obo" is also a pre-packaged option (enter `-g goslim_metagenomics` to specify it). Or \
                a different obo-formatted file can be specified here.',
            default="go_basic"
        )



    ### subcommand cli for getting GO term info ###
    get_term_info_desc = """
        This subcommand gets information for an individual GO term.
        """

    get_term_info_parser = subparsers.add_parser(
        "get-term-info",
        help="Get information on individual GO terms",
        description=get_term_info_desc,
        epilog="Ex. usage: `bit-go get-term-info GO:0004386`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_term_info_required = get_term_info_parser.add_argument_group("Required Parameters")
    get_term_info_optional = get_term_info_parser.add_argument_group("Optional Parameters")

    get_term_info_required.add_argument(
        "GO_term",
        help="GO term to get information on (e.g., GO:0004386)",
        metavar="<GO-TERM>"
    )

    add_common_optional_arguments(get_term_info_optional)

    get_term_info_optional.add_argument(
        "--parents-only",
        help="Add this flag to report parents only, and no children.",
        action="store_true")


    add_help(get_term_info_optional)

    add_version_arg(get_term_info_optional)

    get_term_info_parser.set_defaults(func="get_term_info")


    ### subcommand cli for summarizing GO annotations ###
    summarize_annotations_desc = """
        This subcommand takes an GO-annotations tab-delimited file. By default, it returns one
        tab-delimited summary table with counts and percentages for each GO term. If the `--by-namespace` \
        flag is added, it will also return separate summaries for each GO namespace.
        """

    summarize_annotations_parser = subparsers.add_parser(
        "summarize-annotations",
        help="Summarize GO annotations",
        description=summarize_annotations_desc,
        epilog="Ex. usage: `bit-go summarize-annotations -i GO-annotations.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    summarize_annotations_required = summarize_annotations_parser.add_argument_group("Required Parameters")
    summarize_annotations_optional = summarize_annotations_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(summarize_annotations_required)

    add_common_optional_arguments(summarize_annotations_optional)

    summarize_annotations_optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help='Prefix for output summary file(s) (default: "GO-summary")',
        default="GO-summary"
    )
    summarize_annotations_optional.add_argument(
        "--keep-zeros",
        help="By default the program will remove any rows with for GO terms with 0 counts to them. Add this flag to keep them.",
        action="store_true"
    )

    summarize_annotations_optional.add_argument(
        "--by-namespace",
        help="Add this flag to return separate summaries for each GO namespace.",
        action="store_true"
    )

    add_help(summarize_annotations_optional)

    add_version_arg(summarize_annotations_optional)

    summarize_annotations_parser.set_defaults(func="summarize_annotations")


    ### subcommand cli for combining GO-annotation summaries ###
    combine_summaries_desc = """
        This subcommand takes multiple GO summary tables produced by `bit-go summarize-annotations` and combines them into a single table.
        """

    combine_summaries_parser = subparsers.add_parser(
        "combine-summaries",
        help="Combine GO summary tables",
        description=combine_summaries_desc,
        epilog="Ex. usage: `bit-go combine-summaries -i sample-1-GO-summary.tsv sample-2-GO-summary.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    combine_summaries_required = combine_summaries_parser.add_argument_group("Required Parameters")
    combine_summaries_optional = combine_summaries_parser.add_argument_group("Optional Parameters")

    combine_summaries_required.add_argument(
        "-i",
        "--input-files",
        metavar="<FILE(s)>",
        nargs="+",
        type=str,
        help="Space-delimited list of `bit-go summarize-annotations` output files",
        required=True
    )

    combine_summaries_optional.add_argument(
        "-n",
        "--sample-names",
        metavar="<NAME(s)>",
        nargs="+",
        help='Sample names provided as a space-delimited list, be sure it matches the order of the input files (by default will use basename of input files up to last period)',
        default=[]
    )

    combine_summaries_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output combined summaries (default: "combined-GO-summaries.tsv")',
        default="combined-GO-summaries.tsv"
    )

    add_help(combine_summaries_optional)

    add_version_arg(combine_summaries_optional)

    combine_summaries_parser.set_defaults(func="combine_summaries")


    ### subcommand cli for slimming down GO terms ###
    slim_terms_desc = """
        This subcommand wraps the goatools `map_to_slim.py` program (github.com/tanghaibao/Goatools#map-go-terms-to-goslim-terms).
        See there for more details, and if you use it in your work, be sure to properly cite them :)
        https://www.nature.com/articles/s41598-018-28948-z. It is included here to streamline integration with
        with the GO databases stored with `bit` and programs like `bit-go summarize-annotations`.
        """

    slim_terms_parser = subparsers.add_parser(
        "slim-terms",
        help="Slim down GO annotations to a specified slim obo",
        description=slim_terms_desc,
        epilog="Ex. usage: `bit-go slim-terms -i GO-annotations.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    slim_terms_required = slim_terms_parser.add_argument_group("Required Parameters")
    slim_terms_optional = slim_terms_parser.add_argument_group("Optional Parameters")

    add_common_required_arguments(slim_terms_required)

    slim_terms_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Name for output slimmed annotation file. (default: "GO-slimmed.tsv").',
        default="GO-slimmed.tsv"
    )

    slim_terms_optional.add_argument(
        "-g",
        "--initial-GO-obo-file",
        metavar="<FILE>",
        help='Initial GO obo file holding relationships of all terms used to perform the annotation (e.g. from: geneontology.org/docs/download-ontology/). By default \
             this program will use "go-basic.obo" that is stored with `bit`. Or a different obo-formatted file can be specified here.',
        default="go_basic"
    )

    slim_terms_optional.add_argument(
        "-s",
        "--slimmed-GO-obo-file",
        metavar="<FILE>",
        help='Slimmed GO obo file holding relationships to collapse GO terms (e.g. from: geneontology.org/docs/download-ontology/#subsets;). By default will \
             use "goslim_metagenomics.obo" that is stored with `bit`. Or a different obo-formatted file can be specified here.',
        default="goslim_metagenomics"
    )

    slim_terms_optional.add_argument(
        "-m",
        "--mode",
        help='Set if the slimmer should return only direct ancestors, or all ancestors. Default setting is to return all.',
        choices=["direct", "all"],
        default="all"
    )

    add_help(slim_terms_optional)

    add_version_arg(slim_terms_optional)

    slim_terms_parser.set_defaults(func="slim_terms")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv)==1: # pragma: no cover
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
            if sys.stdin.isatty():
                parser.subparsers.choices[cmd].print_help(sys.stderr)
                sys.exit(0)
            # else: stdin is being piped, fall through to parse_args()
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    from bit.modules.go.go import (get_term_info, summarize_annotations,
                                   combine_summaries, slim_terms)

    func_map = {
        "get_term_info": get_term_info,
        "summarize_annotations": summarize_annotations,
        "combine_summaries": combine_summaries,
        "slim_terms": slim_terms
    }

    func = func_map[args.func]

    func(args)
