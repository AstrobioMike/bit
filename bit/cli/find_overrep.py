import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg,
                            reconstruct_invocation)
from bit.modules.find_overrep import (DEFAULT_K,
                                      DEFAULT_MIN_PERC,
                                      DEFAULT_MIN_FOLD,
                                      DEFAULT_MAX_DIFF,
                                      DEFAULT_MAX_LEN,
                                      DEFAULT_MAX_READS_FOR_DISCOVERY)


def build_parser(parent_subparsers=None):

    desc = """
        This program finds and counts overrepresented sequences in reads without needing to know
        what they are ahead of time, allowing for variation or errors. It's designed 
        for short reads (<= ~300 bp), and currently probably not suitable for long reads.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "find-overrep-seqs",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit find-overrep-seqs -i R1.fastq.gz -I R2.fastq.gz`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-reads",
        metavar="<FILE>",
        required=True,
        help="Input reads (read 1 if paired), fastq or fasta, gzipped or not",
    )

    optional.add_argument(
        "-I",
        "--input-reads-2",
        metavar="<FILE>",
        help="Input read 2 file if paired, fastq or fasta, gzipped or not",
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="overrep-seqs",
        help='Directory for output files (default: "overrep-seqs")',
    )

    optional.add_argument(
        "-O",
        "--output-prefix",
        metavar="<STR>",
        default="",
        help=("String to be prepended to output files (including separator if wanted, "
              "e.g., 'sample-1-'; default: '')"),
    )

    optional.add_argument(
        "-m",
        "--min-perc",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MIN_PERC,
        help=(f"Minimum percent of reads a sequence must be in to be reported (default: {DEFAULT_MIN_PERC}; "
              "those at a consistent place in reads are reported down to a tenth of this)"),
    )

    optional.add_argument(
        "--min-fold",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MIN_FOLD,
        help=("Minimum fold over the count expected from the reads' composition for a k-mer to seed "
              f"a sequence (default: {DEFAULT_MIN_FOLD})"),
    )

    optional.add_argument(
        "-d",
        "--max-diff",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MAX_DIFF,
        help=("Maximum edits (substitutions + indels) allowed when counting sequences in reads, as a "
              f"fraction of the sequence's length (default: {DEFAULT_MAX_DIFF}, e.g., 2 for a 20-25 bp primer)"),
    )

    optional.add_argument(
        "-k",
        "--kmer-size",
        metavar="<INT>",
        type=int,
        default=DEFAULT_K,
        help=f"K-mer size for finding overrepresented sequences, 7-12 (default: {DEFAULT_K})",
    )

    optional.add_argument(
        "--max-len",
        metavar="<INT>",
        type=int,
        default=DEFAULT_MAX_LEN,
        help=f"Maximum length a found sequence is extended to (default: {DEFAULT_MAX_LEN})",
    )

    optional.add_argument(
        "-a",
        "--adapters",
        metavar="<FILE>",
        help=("Fasta of adapter sequences to trim read-through of, as it appears in reads (i.e., the "
              "start of the adapter right after the insert), in addition to the built-in Illumina ones"),
    )

    optional.add_argument(
        "--no-builtin-adapters",
        action="store_true",
        help=("Don't trim standard Illumina adapters before searching for overrepresented sequences (adapters are still inferred from "
              "paired reads, and any given with --adapters are used)"),
    )

    optional.add_argument(
        "-M",
        "--max-reads-for-discovery",
        metavar="<INT>",
        type=int,
        default=DEFAULT_MAX_READS_FOR_DISCOVERY,
        help=("Number of reads (per file) used for finding the overrepresented sequences; counting them "
              f"is always done across all reads (default: {DEFAULT_MAX_READS_FOR_DISCOVERY:,}; set 0 to use all)"),
    )

    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=5,
        help="Number of threads to use (default: 5)",
    )

    add_force(optional)

    add_help(optional)

    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if not 0 < args.min_perc <= 100:
        parser.error("--min-perc must be greater than 0 and at most 100")
    if args.min_fold < 1:
        parser.error("--min-fold must be 1 or greater")
    if not 0 <= args.max_diff < 0.5:
        parser.error("--max-diff must be at least 0 and less than 0.5")
    if not 7 <= args.kmer_size <= 12:
        parser.error("--kmer-size must be between 7 and 12")
    if args.max_len < args.kmer_size:
        parser.error("--max-len must be at least the k-mer size")
    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater")
    if args.max_reads_for_discovery < 0:
        parser.error("--max-reads-for-discovery must be 0 or greater")

    from bit.modules.general import check_files_are_found, notify_premature_exit
    from bit.modules.find_overrep import setup_output_dir, write_outputs, summary_lines

    read_paths = [args.input_reads] + ([args.input_reads_2] if args.input_reads_2 else [])
    check_files_are_found(read_paths + ([args.adapters] if args.adapters else []))
    setup_output_dir(args.output_dir, args.output_prefix, args.force_overwrite,
                     reconstruct_invocation(parser, args))

    try:
        result = run(args, read_paths)
    except (ValueError, OSError) as e:
        print(f"\n    {e}")
        notify_premature_exit()

    paths = write_outputs(result, args.output_dir, args.output_prefix)

    print()
    for line in summary_lines(result):
        print(f"    {line}")
    print()
    print(f"    Overrepresented sequences written to: '{paths['seqs']}' (and '{paths['fasta']}')")
    print(f"    Where they start in reads written to: '{paths['starts']}'")
    print(f"    Adapter read-through counts written to: '{paths['adapters']}'")
    print()


def run(args, read_paths):

    from bit.modules.find_overrep import find_overrep, load_adapters

    return find_overrep(
        read_paths,
        k=args.kmer_size,
        min_perc=args.min_perc,
        min_fold=args.min_fold,
        max_diff=args.max_diff,
        max_len=args.max_len,
        adapters=load_adapters(args.adapters, use_builtin=not args.no_builtin_adapters),
        jobs=args.jobs,
        max_reads_for_discovery=args.max_reads_for_discovery,
        verbose=True,
    )
