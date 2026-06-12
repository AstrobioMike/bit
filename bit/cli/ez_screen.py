import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_common_snakemake_arguments,
                            reconstruct_invocation,
                            add_version_arg)


def build_parser(parent_subparsers=None):

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
            help = 'Minimum percent coverage of a target to be counted (default: 80)',
            metavar = "<INT>",
            default = 80,
            type = float
        )

    ### subcommand cli for assembly screening ###
    assembly_description = ("This subcommand searches input assemblies for a set of target genes or regions. "
                            "It generates summary tables of identified targets overall, and per contig, and "
                            "a region-based output that cleans up redundancies of multiple targets hitting "
                            "the same loci.")

    assembly_parser = subparsers.add_parser(
        "assembly",
        help="Run BLAST/DIAMOND-based screening of targets in assemblies",
        description=assembly_description,
        epilog="Ex. usage: `bit ez-screen assembly -a assembly.fasta -t targets.fasta` "
               "OR `bit ez-screen assembly -a assembly.fasta -t nt-targets.fasta prot-targets.faa`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    assembly_required = assembly_parser.add_argument_group('Required Parameters')
    assembly_optional = assembly_parser.add_argument_group('Optional Parameters')

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

    add_common_optional_arguments(assembly_optional)

    assembly_optional.add_argument(
        "--min-nt-perc-id",
        help = 'Minimum percent ID for a nucleotide hit to be counted (blastn; default: 80)',
        metavar = "<INT>",
        default = 80,
        type = float
    )

    assembly_optional.add_argument(
        "--min-aa-perc-id",
        help = 'Minimum percent ID for a protein hit to be counted (DIAMOND blastx; default: 70)',
        metavar = "<INT>",
        default = 70,
        type = float
    )

    assembly_optional.add_argument(
        "-n",
        "--num-threads",
        help = "Number of threads to use during search (other than nt fasta input targets; default: 5)",
        metavar = "<INT>",
        default = 5,
        type = int
    )

    assembly_optional.add_argument(
        "--hit-merge-gap",
        help="When counting hits per contig, multiple alignments (HSPs) of the same "
             "target separated by up to this many bp are counted as one hit "
             "(one locus), so fragmented alignments don't inflate 'num_total_hits' "
             "(default: 200)",
        metavar="<INT>",
        default=200,
        type=int)

    assembly_optional.add_argument(
        "--dont-resolve-regions",
        help="By default, hits from different targets at the same contig "
             "locus are collapsed into a single called region, keeping the best hit "
             "and reporting things in a '-region-calls.tsv' file. Add this flag to "
             "disable that",
        dest="resolve_regions",
        action="store_false")

    assembly_optional.add_argument(
        "--region-overlap-frac",
        help="Two hits at the same contig locus are treated as the same region when they "
             "overlap by at least this fraction of the shorter hit's span (default: 0.5)",
        metavar="<FLOAT>",
        default=0.5,
        type=float)

    assembly_optional.add_argument(
        "--assemblies-as-rows",
        help = 'Add this flag to have the summary table have assemblies as rows and targets as columns',
        action = "store_true")

    assembly_optional.add_argument(
        "--report-all-targets",
        help = "Add this flag to report all targets in the output summary table (even those with 0 hits detected)",
        action = "store_true")

    add_help(assembly_optional)

    add_version_arg(assembly_optional)

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

def main():

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
