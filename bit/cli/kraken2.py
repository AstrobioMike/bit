#!/usr/bin/env python
import os
import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg


def build_parser():

    desc = """
        This program provides subcommands for working with kraken2 (or bracken) output. See subcommand-specific
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

    ####################################################
    ### subcommand cli for generating taxonomy plots ###
    ####################################################
    tax_plots_desc = """
        This subcommand takes a kraken2 report as input and generates bar plots for the most abundant taxa,
        or those above a specified percent threshold, at each rank.
        """

    tax_plots_parser = subparsers.add_parser(
        "tax-plots",
        help="Generate taxonomy bar plots from a kraken2 report",
        description=tax_plots_desc,
        epilog="Ex. usage: `bit-kraken2 tax-plots -i kraken2.report`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    tax_plots_required = tax_plots_parser.add_argument_group("Required Parameters")
    tax_plots_optional = tax_plots_parser.add_argument_group("Optional Parameters")

    tax_plots_required.add_argument(
        "-i",
        "--input-kraken2-report",
        metavar="<FILE>",
        help="Input kraken2 report",
        required=True,
        dest="input_kraken"
    )

    tax_plots_optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help='Prefix for output plots (default: "")',
        default=""
    )

    tax_plots_optional.add_argument(
        "-m",
        "--maximum-taxa-per-plot",
        metavar="<INT>",
        help='Maximum number of unique taxa in a plot; the rest are grouped into "Other" (default: 10)',
        type=int,
        default=10,
        dest="max_taxa"
    )

    tax_plots_optional.add_argument(
        "-p",
        "--minimum-percent-threshold",
        metavar="<FLOAT>",
        help='Minimum percent of reads a taxon must have to be shown; overrides --maximum-taxa-per-plot when > 0 '
             '(default: 0.0)',
        type=float,
        default=0.0,
        dest="min_percent"
    )

    tax_plots_optional.add_argument(
        "--no-annots",
        metavar="",
        help="Don't add annotations (total reads, %% unclassified, %% stuck at root) as overlain text on plots",
    )

    add_help(tax_plots_optional)

    add_version_arg(tax_plots_optional)

    tax_plots_parser.set_defaults(func="tax_plots")


    #########################################################
    ### subcommand cli for generating a tax summary table ###
    #########################################################
    tax_summary_desc = """
        This subcommand takes one or more kraken2 (or bracken) reports and generates a combined table with
        full-lineage information, counts, and percentages. The lineage info is constructed directly from the
        report files. Completely unclassified things will be reported as "Unclassified" at each rank. If you
        see NA for all ranks, it's likely from something like "cellular organisms" or "root" — check the taxid
        in the original report to find out.
        """

    tax_summary_parser = subparsers.add_parser(
        "tax-summary",
        help="Generate a taxonomy summary table from kraken2/bracken report(s)",
        description=tax_summary_desc,
        epilog="Ex. usage: `bit-kraken2 tax-summary -i kraken2.report`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    tax_summary_required = tax_summary_parser.add_argument_group("Required Parameters")
    tax_summary_optional = tax_summary_parser.add_argument_group("Optional Parameters")

    tax_summary_required.add_argument(
        "-i",
        "--input-reports",
        metavar="<FILE(s)>",
        nargs="+",
        help="Input kraken2 or bracken report file(s) (space-delimited if multiple)",
        required=True,
    )

    tax_summary_optional.add_argument(
        "-s",
        "--sample-names",
        metavar="<NAME(s)>",
        nargs="+",
        help="Wanted sample name(s) (space-delimited if multiple); must match the order of the input files "
             "(by default uses the basename of each input file up to the last period)",
        default=None
    )

    tax_summary_optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help="Output TSV file (default: tax-summary.tsv)",
        default="tax-summary.tsv",
    )

    add_help(tax_summary_optional)

    add_version_arg(tax_summary_optional)

    tax_summary_parser.set_defaults(func="tax_summary")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

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
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    func_map = {
        "tax_summary": run_tax_summary,
        "tax_plots": run_tax_plots,
    }

    func_map[args.func](args)


def run_tax_summary(args):

    from bit.modules.general import check_files_are_found
    from bit.modules.kraken2_tax_summary import kraken2_tax_summary

    check_files_are_found(args.input_reports)

    input_map = check_and_setup_file_to_sample_map(args.input_reports, args.sample_names)

    kraken2_tax_summary(input_map, args.output_tsv)



def check_and_setup_file_to_sample_map(input_reports, sample_names):

    from bit.modules.general import report_message, notify_premature_exit

    input_map = {}

    if sample_names is None or len(sample_names) == 0:

        for file in input_reports:
            curr_sample = os.path.basename(file).rsplit('.', 1)[0]
            if curr_sample in input_map.values():
                report_message('\n    It seems the file "' + file + '" would have an overlapping basename with another input report. That\'s not gonna fly.\n')
                report_message("\n    You can provide custom sample names with the `-s/--sample-names` flag to get around this.\n", color="none")
                notify_premature_exit()
            input_map[file] = curr_sample


    if sample_names is not None and len(sample_names) != len(input_reports):
        report_message("\n    It seems the number of provided sample names doesn't match the number of provided input files :(")
        report_message("\n    Check usage with `bit-kraken2-tax-summary -h`.\n")
        notify_premature_exit()

    if sample_names is not None:
        for file, name in zip(input_reports, sample_names):
            if name in input_map.values():
                report_message('\n    It seems the name "' + name + '" was given more than once in the provided sample names. That\'s not gonna fly.\n')
                report_message("\n    Please make sure all sample names provided with the `-s/--sample-names` flag are unique.\n", color="none")
                notify_premature_exit()
            input_map[file] = name

    return input_map


def run_tax_plots(args):

    from bit.modules.general import check_files_are_found
    from bit.modules.kraken2_tax_plots import gen_kraken2_tax_plots

    check_files_are_found([args.input_kraken])

    gen_kraken2_tax_plots(args.input_kraken, args.output_prefix, args.max_taxa, args.min_percent, args.no_annots)
