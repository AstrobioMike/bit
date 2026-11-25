import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_seed,
                            add_version_arg)


def build_parser(parent_subparsers=None):

    desc = """
        This program generates reads in FASTQ format from one or multiple input FASTA files,
        without simulating sequencing errors. The genomes can first be mutated (see --mutation-rate and others),
        to introduce some variation.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "gen-reads",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit gen-reads -i genome.fasta`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    general = parser.add_argument_group("General Parameters")
    paired = parser.add_argument_group("Paired-end Parameters")
    long = parser.add_argument_group("Long-read Parameters")

    required.add_argument(
        "-i",
        "--input-fastas",
        metavar="<FILES>",
        nargs='+',
        help="List of input FASTA files",
        required=True,
    )
    general.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help="Prefix for output FASTQ files (default: generated-reads)",
        default="generated-reads",
    )
    general.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=10,
        help="Number of parallel processes to use for read generation when multiple input fastas are provided (default: 10)",
    )
    general.add_argument(
         "-t",
         "--type",
         help="Type of reads to generate (default: 'paired-end')",
         choices=["paired-end", "single-end", "long"],
         default="paired-end",
    )
    general.add_argument(
        "-n",
        "--num-reads",
        metavar="<INT>",
        help="Number of total reads to generate ([bold]NOT[/bold] read-pairs; default: 1000000)",
        type=int,
        default=1000000,
    )
    general.add_argument(
        "-r",
        "--read-length",
        metavar="<INT>",
        help="Length of each read (default: 150 for paired-end and single-end; 5000 for long reads)",
        type=int,
        default=None,
    )
    general.add_argument(
        "-c",
        "--coverage",
        metavar="<INT or FILE>",
        help="Wanted coverage(s). This can be just an integer if wanting the same coverage for \
              all input fasta files, or a coverages file (tab-delimited) specifying wanted \
              coverages per input fasta (1st column should be path to input fasta; 2nd column \
              should be wanted coverage for each input). This overrides the --num-reads option. \
              Cannot be combined with the --proportions-file option.",
    )
    general.add_argument(
        "-p",
        "--proportions-file",
        metavar="<FILE>",
        help="Proportions file (tab-delimited) specifying read-proportions per input fasta \
              (1st column should be path to input fasta; 2nd column should be wanted \
              proportion for each input with the total summing to one). If not provided, \
              equal proportions are assumed. Cannot be combined with the --coverage option.",
    )
    general.add_argument(
        "-C",
        "--circularize",
        help="Treat input contigs as circular and allow fragments that span the end-to-start \
              boundary to be generated. (You probably don't want this if your inputs include \
              linear genomes or if they're split across multiple contigs.)",
        action="store_true",
        default=False,
    )
    general.add_argument(
        "--include-Ns",
        help="By default, regions with Ns are skipped. Provide this flag to generate reads that might include Ns.",
        action="store_true",
        default=False,
    )

    general.add_argument(
        "--per-read-tsv",
        help="Add this flag to write a per-read provenance tsv; when set, read headers are kept minimal, "
             "otherwise the same info is in each read's header",
        action="store_true",
        default=False,
    )

    general.add_argument(
        "--ground-truth-assembly",
        help="Add this flag to write out a ground-truth assembly fasta (genomes may be fragmented if reads are generated with low coverage)",
        action="store_true",
        default=False,
    )

    add_seed(general)
    add_help(general)
    add_version_arg(general)

    mutation = parser.add_argument_group("Mutation Parameters")
    mutation.add_argument(
        "--mutation-rate",
        metavar="<FLOAT>",
        help="Mutate each input genome at this rate (0-1, e.g. 0.01 for about 99%% identity) before generating reads "
             "from it; the mutated genomes are written to '<output prefix>-mutated-genomes/', "
             "and the per-read sources and ground-truth assembly refer to them (default: 0, no mutation)",
        type=float,
        default=0.0,
    )
    mutation.add_argument(
        "--ti-tv-ratio",
        metavar="<FLOAT>",
        help="Transition/transversion ratio of the substitutions (default: 1.0)",
        type=float,
        default=1.0,
    )
    mutation.add_argument(
        "--indel-rate",
        metavar="<FLOAT>",
        help="Share of the mutations that are indels (0-1; e.g. 0.2 is 1 in 5; insertions and deletions equally "
             "likely), which shift positions relative to the input genome (default: 0.0)",
        type=float,
        default=0.0,
    )

    paired.add_argument(
        "-f",
        "--fragment-size",
        metavar="<INT>",
        help="Size of the fragment from which paired-end reads are generated (default: 500; "
            "these may be shorter than the specified size when shorter contigs are present)",
        type=int,
        default=500,
    )
    paired.add_argument(
        "-F",
        "--fragment-size-range",
        metavar="<INT>",
        help="Percent range (+/-) around --fragment-size when using paired-end mode (default: 10). E.g., 10 with "
             "--fragment-size 500 generates fragments from 450 to 550",
        type=int,
        default=10,
    )

    long.add_argument(
        "-L",
        "--long-read-length-range",
        metavar="<INT>",
        help="Percent range (+/-) around --read-length when using long-read mode (default: 50). E.g., 50 with "
             "--read-length 5000 generates reads from 2500 to 7500",
        type=int,
        default=50,
    )

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    for name in ("mutation_rate", "indel_rate"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 1")

    if args.type == "long" and not args.read_length:
        args.read_length = 5000
    elif not args.read_length:
        args.read_length = 150

    from bit.modules.gen_reads import generate_reads
    generate_reads(args)
