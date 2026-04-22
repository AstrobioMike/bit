import sys
import argparse
from bit.modules.gen_reads import generate_reads
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_seed)


def main():

    desc = """
        This script generates perfect (no error model) reads in FASTQ format from one or
        multiple input FASTA files. See `bit-mutate-seqs` if wanting to introduce variation
        to a fasta prior to read-generation. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-gen-reads -i genome.fasta`",
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
        help="Prefix for output FASTQ files (default: perfect-reads)",
        default="perfect-reads",
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

    add_seed(general)

    add_help(general)

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

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.type == "long" and not args.read_length:
        args.read_length = 5000
    elif not args.read_length:
        args.read_length = 150

    generate_reads(args)
