import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            wrap_help,
                            add_common_snakemake_arguments,
                            reconstruct_invocation,
                            add_version_arg)


def build_parser(parent_subparsers=None, show_fine=False):
    """ builds the ez-screen parser.

        show_fine controls whether the assembly subcommand's 'Fine-tuning
        Parameters' group is included. By default (show_fine=False) it is
        omitted so the standard `-h` menu stays uncluttered; the group is added
        only when building the detailed help (triggered by -H/--show-detailed-help,
        handled in main()). The fine-tuning options still WORK when passed either
        way -- when omitted here they're added in a hidden pass so argparse still
        recognizes and defaults them (see _add_fine_tuning_arguments). """

    desc = """
        This program helps detect target-genes/regions present in assemblies or reads. See subcommand-specific
        help menus for more info.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "ez-screen",
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

    def add_common_optional_arguments(group):
        group.add_argument(
            "-o",
            "--output-prefix",
            help = 'Output prefix (default: "ez-screen")',
            metavar = "<STR>",
            default = "ez-screen", type = str
        )
        group.add_argument(
            "-M",
            "--min-perc-cov",
            help = 'Minimum percent coverage required of a target (default: 80)',
            metavar = "<INT>",
            default = 80,
            type = float
        )

    ### subcommand cli for assembly screening ###
    assembly_description = ("This subcommand searches input assemblies for given target sequences. "
                            "It generates summary tables of identified targets overall, and per contig, "
                            "a region-based output that cleans up redundancies of multiple targets hitting "
                            "the same loci, and it extracts \"islands\" of contiguous region stretches in fasta format.")

    assembly_parser = subparsers.add_parser(
        "assembly",
        help="Run BLAST/DIAMOND-based screening of targets in assemblies",
        description=assembly_description,
        epilog="Ex. usage: `bit ez-screen assembly -a assembly.fasta -t targets.fasta` "
               "OR `bit ez-screen assembly -a assembly.fasta -t nt-targets-blastdb-prefix prot-targets.dmnd`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    assembly_required = assembly_parser.add_argument_group('Required Parameters')
    assembly_general = assembly_parser.add_argument_group('General Parameters')

    assembly_required.add_argument(
        "-a",
        "--assemblies",
        help = "Assembly files in fasta format",
        metavar = "<FILE>",
        required = True,
        nargs = '+')


    assembly_required.add_argument(
        "-t",
        "--targets",
        help = "Nucleotide targets (fasta or BLAST db) and/or protein targets (fasta or DIAMOND db) to search for",
        metavar = "<FILE/PREFIX>",
        required = True,
        nargs = '+')

    add_common_optional_arguments(assembly_general)

    assembly_general.add_argument(
        "--min-nt-perc-id",
        help = 'Minimum percent ID for a nucleotide target (default: 80)',
        metavar = "<INT>",
        default = 80,
        type = float
    )

    assembly_general.add_argument(
        "--min-aa-perc-id",
        help = 'Minimum percent ID for a protein target (default: 70)',
        metavar = "<INT>",
        default = 70,
        type = float
    )

    assembly_general.add_argument(
        "-n",
        "--num-threads",
        help = "Number of threads to use during database searches (default: 5)",
        metavar = "<INT>",
        default = 5,
        type = int
    )

    assembly_general.add_argument(
        "-h",
        "--help",
        action = "help",
        help = wrap_help("Show basic help message and exit")
    )

    assembly_general.add_argument(
        "-H",
        "--detailed-help",
        help = "Show detailed help, including fine-tuning parameters",
        action = "store_true")

    if show_fine:
        # visible fine-tuning group, with real help strings, for the detailed menu
        assembly_fine = assembly_parser.add_argument_group('Fine-tuning Parameters')
        _add_fine_tuning_arguments(assembly_fine, hidden=False)
    else:
        # standard menu: omit the group from help entirely, but still register
        # the same arguments (hidden) so they parse and get their defaults
        _add_fine_tuning_arguments(assembly_parser, hidden=True)

    add_version_arg(assembly_general)

    assembly_parser.set_defaults(func="run_assembly")

    ### subcommand cli for reads screening ###
    reads_description = ("This subcommand maps input reads (currently paired-end only; bwa mem) against a set of target genes or regions (in "
                            "nucleotide format). It generates mapping info in the form of BAM files and summary tables "
                            "that report for each target: how many reads mapped to it; what their average identity was; "
                            "and what proportion of the target was covered by the reads. Results will only be returned for "
                            "targets for which the reads mapping to them surpass the tunable minimum target coverage and average percent-identity thresholds.")
    reads_parser = subparsers.add_parser(
        "reads",
        help="Run mapping-based screening of reads against targets",
        description=reads_description,
        epilog="Ex. usage: `bit ez-screen reads -t targets.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    reads_required = reads_parser.add_argument_group('Required Parameters')
    reads_optional = reads_parser.add_argument_group('Optional Parameters')
    reads_snakemake = reads_parser.add_argument_group('Snakemake Parameters')

    reads_required.add_argument(
        "-t",
        "--targets",
        help = "Targets you want to search for, e.g. genes/regions (as nucleotide fasta)",
        metavar = "<FILE>",
        required = True
    )

    reads_optional.add_argument(
        "-r",
        "--reads-dir",
        help = "Directory holding the input gzipped reads (default: current directory)",
        metavar = "<DIR>",
        action = "store",
        default = ".",
        type = str
    )

    reads_optional.add_argument(
        "-m",
        "--min-perc-id",
        help = 'Minimum percent ID for a hit to be counted (default: 80)',
        metavar = "<INT>",
        default = 80,
        type = float
    )

    add_common_optional_arguments(reads_optional)

    add_help(reads_optional)

    add_version_arg(reads_optional)

    add_common_snakemake_arguments(reads_snakemake)

    reads_parser.set_defaults(func="run_reads")

    return parser


def _add_fine_tuning_arguments(target, hidden=False):
    """ adds the assembly subcommand's fine-tuning arguments to `target` (either
        an argument group, for the visible detailed menu, or the assembly parser
        itself with help suppressed, for the standard menu).

        hidden=True suppresses each from the help display while keeping them
        fully functional, so `--island-gap 3000` etc. still parse on the standard
        menu even though they aren't listed there. """

    def h(text):
        return argparse.SUPPRESS if hidden else text

    target.add_argument(
        "--hit-merge-gap",
        help=h("When counting hits per contig, multiple alignments (HSPs) of the same "
               "target separated by up to this many bp are counted as one hit/locus "
               "(default: 200)"),
        metavar="<INT>",
        default=200,
        type=int)

    target.add_argument(
        "--no-region-resolution",
        help=h("By default, hits from different targets at the same contig "
               "locus are collapsed into a single region. Add this flag to disable this"),
        dest="resolve_regions",
        action="store_false")

    target.add_argument(
        "--region-overlap-frac",
        help=h("Two hits at the same contig locus are treated as the same region when they "
               "overlap by at least this fraction of the shorter hit's span (default: 0.5)"),
        metavar="<FLOAT>",
        default=0.5,
        type=float)

    target.add_argument(
        "--no-island-extraction",
        help=h("By default, a dense cluster of target regions (an 'island') will be extracted "
               "into its own fasta based on below parameters. Add this flag to disable this"),
        dest="island_extraction",
        action="store_false")

    target.add_argument(
        "--island-gap",
        help=h("Resolved regions within this many bp of each other on a contig are "
               "chained into one island (default: 2500)"),
        metavar="<INT>",
        default=2500,
        type=int)

    target.add_argument(
        "--island-contig-floor",
        help=h("Islands are only extracted from contigs at least this many bps long (default: 20000)"),
        metavar="<INT>",
        default=20000,
        type=int)

    target.add_argument(
        "--island-contig-ratio",
        help=h("An island is only extracted if its contig is at least this many times "
               "longer than said island (default: 3)"),
        metavar="<FLOAT>",
        default=3.0,
        type=float)

    target.add_argument(
        "--island-min-span",
        help=h("Minimum bps of an island for it to be extracted (default: 2000)"),
        metavar="<INT>",
        default=2000,
        type=int)

    target.add_argument(
        "--island-min-regions",
        help=h("Minimum number of resolved regions in an island for it to be extracted "
               "(default: 3)"),
        metavar="<INT>",
        default=3,
        type=int)

    target.add_argument(
        "--island-buffer",
        help=h("When cutting out an island, include this many bp of flanking contig "
               "on each side (default: 500)"),
        metavar="<INT>",
        default=500,
        type=int)

    target.add_argument(
        "--assemblies-as-rows",
        help=h('Add this flag to have the summary table have assemblies as rows and targets as columns'),
        action="store_true")

    target.add_argument(
        "--report-all-targets",
        help=h("Add this flag to report all targets in the output summary table (even those with 0 hits detected)"),
        action="store_true")


def main():

    # detailed help: -H/--show-detailed-help anywhere on an assembly invocation
    # prints the full menu (including fine-tuning params) and exits, mirroring how
    # -h works regardless of position and without requiring the otherwise-required
    # -a/-t. handled before normal parsing so an incomplete command line still
    # shows help rather than erroring on missing required args.
    argv = sys.argv[1:]
    if "assembly" in argv and ("-H" in argv or "--show-detailed-help" in argv):
        detailed_parser = build_parser(show_fine=True)
        detailed_parser.subparsers.choices["assembly"].print_help(sys.stderr)
        sys.exit(0)

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv)==1:
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

    from bit.modules.ez_screen import run_assembly, run_reads

    func_map = {
        "run_assembly": run_assembly,
        "run_reads": run_reads
    }

    func = func_map[args.func]

    # reconstructing the full command-line invocation for logging
    full_cmd_executed = reconstruct_invocation(parser, args)
    func(args, full_cmd_executed)


if __name__ == "__main__":
    main()
