import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            ParagraphRichHelpFormatter,
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


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        This program finds and counts overrepresented sequences in reads while allowing for variation or errors.
        It's designed for short reads (up to a few hundred bases) and currently isn't suitable for longer reads.
        Run with `-H` to see the detailed help menu, including how the program works and fine-tuning parameters.
        """

    detailed_desc = desc + """

        HOW IT WORKS: Reads are read through five times. (1) All k-mers (default k=11) are counted, and
        for paired reads, adapters not among the built-in or given ones are inferred. (2) K-mers are counted
        again with adapters trimmed off, and those in at least --min-perc percent of reads and at least
        --min-fold times what the reads' composition predicts, as well as any that are in at least 1% of reads,
        become candidate seed k-mers. (3) Rarer k-mers (down to a tenth of --min-perc) that tend to start at the
        same spot in the reads become candidates too. (4) K-mer candidates are extended into full sequences, first from
        the k-mer counts, and then using the reads that hold them. (5) These full sequences (allowing up to --max-diff
        edits per base) are counted across all reads, along with where they start in the reads, how often
        they're reverse-complemented, and which searched sequences may hold other sequences (e.g., primer dimers).
        Passes 1-4 use up to --max-reads-for-discovery reads per file, while pass 5 always uses all of them.
        """

    description = detailed_desc if show_fine else desc
    formatter = ParagraphRichHelpFormatter if show_fine else CustomRichHelpFormatter

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "find-overrep-seqs",
            description=description,
            formatter_class=formatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=description,
            epilog="Ex. usage: `bit find-overrep-seqs -i R1.fastq.gz -I R2.fastq.gz`",
            formatter_class=formatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    general = parser.add_argument_group("General Parameters")

    required.add_argument(
        "-i",
        "--input-reads",
        metavar="<FILE>",
        required=True,
        help="Input reads (read 1 if paired), fastq or fasta, gzipped or not",
    )

    general.add_argument(
        "-I",
        "--input-reads-2",
        metavar="<FILE>",
        help="Input read 2 file if paired, fastq or fasta, gzipped or not",
    )

    general.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="overrep-seqs",
        help='Directory for output files (default: "overrep-seqs")',
    )

    general.add_argument(
        "-O",
        "--output-prefix",
        metavar="<STR>",
        default="",
        help=("String to be prepended to output files (including separator if wanted, "
              "e.g., 'sample-1-'; default: '')"),
    )

    general.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=5,
        help="Number of threads to use (default: 5)",
    )

    add_force(general)

    general.add_argument(
        "-H",
        "--detailed-help",
        action="store_true",
        help="Show detailed help, including how it works and fine-tuning parameters",
    )

    add_help(general)

    add_version_arg(general)

    if show_fine:
        # visible fine-tuning group, with real help strings, for the detailed menu
        fine = parser.add_argument_group("Fine-tuning Parameters")
        _add_fine_tuning_arguments(fine, hidden=False)
    else:
        # standard menu: leave the group out of the help entirely, but still register the same
        # arguments (hidden) so they parse and get their defaults
        _add_fine_tuning_arguments(parser, hidden=True)

    return parser


def _add_fine_tuning_arguments(target, hidden=False):

    def h(text):
        return argparse.SUPPRESS if hidden else text

    target.add_argument(
        "-m",
        "--min-perc",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MIN_PERC,
        help=h(f"Minimum percent of reads a sequence must be in to be reported (default: {DEFAULT_MIN_PERC}; "
               "those at a consistent place in reads are reported down to a tenth of this)"),
    )

    target.add_argument(
        "--min-fold",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MIN_FOLD,
        help=h("Minimum fold over the count expected from the reads' composition for a k-mer to seed "
               f"a sequence (default: {DEFAULT_MIN_FOLD})"),
    )

    target.add_argument(
        "-d",
        "--max-diff",
        metavar="<FLOAT>",
        type=float,
        default=DEFAULT_MAX_DIFF,
        help=h("Maximum edits (substitutions + indels) allowed when counting sequences in reads, as a "
               f"fraction of the sequence's length (default: {DEFAULT_MAX_DIFF}, e.g., 2 for a 20-25 bp primer)"),
    )

    target.add_argument(
        "-k",
        "--kmer-size",
        metavar="<INT>",
        type=int,
        default=DEFAULT_K,
        help=h(f"K-mer size for finding overrepresented sequences, 7-12 (default: {DEFAULT_K}; "
               "smaller finds shorter sequences but with more chance matches, larger the reverse)"),
    )

    target.add_argument(
        "--max-len",
        metavar="<INT>",
        type=int,
        default=DEFAULT_MAX_LEN,
        help=h(f"Maximum length a found sequence is extended to (default: {DEFAULT_MAX_LEN})"),
    )

    target.add_argument(
        "-a",
        "--adapters",
        metavar="<FILE>",
        help=h("Fasta of adapter sequences to trim read-through of, as it appears in reads (i.e., the "
               "start of the adapter right after the insert), in addition to the built-in Illumina ones"),
    )

    target.add_argument(
        "--no-builtin-adapter-trim",
        action="store_true",
        help=h("Don't trim standard Illumina adapters before searching for overrepresented sequences "
               "(adapters are still inferred from paired reads, and any given with --adapters are used)"),
    )

    target.add_argument(
        "-M",
        "--max-reads-for-discovery",
        metavar="<INT>",
        type=int,
        default=DEFAULT_MAX_READS_FOR_DISCOVERY,
        help=h("Number of reads (per file) used for finding the overrepresented sequences; counting them "
               f"is always done across all reads (default: {DEFAULT_MAX_READS_FOR_DISCOVERY:,}; set 0 to use all)"),
    )


def main():

    argv = sys.argv[1:]
    if "-H" in argv or "--detailed-help" in argv:
        build_parser(show_fine=True).print_help(sys.stderr)
        sys.exit(0)

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
        adapters=load_adapters(args.adapters, use_builtin=not args.no_builtin_adapter_trim),
        jobs=args.jobs,
        max_reads_for_discovery=args.max_reads_for_discovery,
        verbose=True,
    )
