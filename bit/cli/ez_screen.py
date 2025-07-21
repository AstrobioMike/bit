import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_common_snakemake_arguments,
                            reconstruct_invocation)
from bit.modules.ez_screen import run_assembly, run_reads


def build_parser():
    parser = argparse.ArgumentParser(
        description="This program helps detect target genes/regions present in assemblies or reads",
        epilog="Use `bit-ez-screen assembly -h` or `bit-ez-screen reads -h` to see subcommand-specific help.",
        formatter_class=CustomRichHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "-t", "--targets",
            help = "Targets you want to search for, e.g. genes/regions (as nucleotide fasta)",
            metavar = "<FILE>", required = True
        )

    def add_common_optional_arguments(group):
        group.add_argument(
            "-o", "--output-prefix",
            help = 'Output prefix (default: "ez-screen")',
            metavar = "<STR>", default = "ez-screen", type = str
        )
        group.add_argument(
            "-m", "--min-perc-id",
            help = 'Minimum percent ID for a hit to be counted (default: 80)',
            metavar = "<INT>", default = 80, type = float
        )
        group.add_argument(
            "-M", "--min-perc-cov",
            help = 'Minimum percent coverage of a target to be counted (default: 80)',
            metavar = "<INT>", default = 80, type = float
        )

    ### subcommand cli for assembly screening ###
    assembly_description = ("This subcommand takes a set of target genes or regions to find (in nucleotide format) "
                            "and searches for them in input assemblies with blastn. It generates the general BLAST outputs "
                            "as well as simplified summary tables that report how many times each target was "
                            "found in each input assembly (if the target was < 10,000 bases), or if the target was "
                            "detected at all (if the target was >= 10,000 bases), based on tunable minimum target "
                            "coverage and percent-identity thresholds.")
    assembly_parser = subparsers.add_parser(
        "assembly",
        help="Run BLAST-based screening of targets in assemblies",
        description=assembly_description,
        epilog="Ex. usage: `bit-ez-screen assembly -a assembly.fasta -t targets.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    assembly_required = assembly_parser.add_argument_group('REQUIRED PARAMETERS')
    assembly_optional = assembly_parser.add_argument_group('OPTIONAL PARAMETERS')

    assembly_required.add_argument("-a", "--assemblies", help = "Assembly files in fasta format", metavar = "<FILE>", required = True,
                        nargs = '+')

    add_common_required_arguments(assembly_required)
    add_common_optional_arguments(assembly_optional)


    assembly_optional.add_argument("-T", "--transpose-output-tsv", help = 'Set this flag if we want to have the output table have targets as rows rather than columns.', action = "store_true")

    assembly_optional.add_argument("-f", "--filter-if-not-detected",
                        help = "By default, all targets are included in the output table, even if they weren't detected \
                                in any input assemblies. Add this flag if you'd like to filter them out of the final output table.",
                        action = "store_true")

    add_help(assembly_optional)

    assembly_parser.set_defaults(func=run_assembly)

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
        epilog="Ex. usage: `bit-ez-screen reads -s sample-1 -t targets.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    reads_required = reads_parser.add_argument_group('REQUIRED PARAMETERS')
    reads_optional = reads_parser.add_argument_group('OPTIONAL PARAMETERS')
    reads_snakemake = reads_parser.add_argument_group('SNAKEMAKE PARAMETERS')

    add_common_required_arguments(reads_required)

    reads_optional.add_argument("-r", "--reads-dir", help = "Directory holding the input gzipped reads (default: current directory)",
                        metavar = "<DIR>", action = "store", default = ".", type = str)

    add_common_optional_arguments(reads_optional)

    add_help(reads_optional)

    add_common_snakemake_arguments(reads_snakemake)

    reads_parser.set_defaults(func=run_reads)

    return parser

def main():
    parser = build_parser()
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so approriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]
        if cmd in parser.subparsers.choices:
            parser.subparsers.choices[cmd].print_help(sys.stderr)
            sys.exit(0)
        # unknown command: fall back to global help and error
        print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # reconstructing the full command-line invocation for logging
    full_cmd_executed = reconstruct_invocation(parser, args)
    args.func(args, full_cmd_executed)


if __name__ == "__main__":
    main()
