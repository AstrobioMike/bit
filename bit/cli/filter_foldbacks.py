import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            ParagraphRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg,
                            reconstruct_invocation)
from bit.modules.filter_foldbacks import (KMER_SIZE,
                                          MIN_SHARED,
                                          MIN_LENGTH,
                                          LONG_MIN_LENGTH)


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        This program removes fold-back artifacts of PTA/MDA in short Illumina reads (default) or long reads
        (--hifi)). Run with `-H` to see the detailed help menu, including how the
        program works and fine-tuning parameters.
        """



    detailed_desc = desc + """

        The most common artifact of PTA (and frequent in MDA) is a fold-back, where partway through being
        copied, the molecule folds back on itself and the polymerase carries on copying the strand it just
        made. This results in the output read holding a stretch of sequence followed by its own reverse complement.
        This program cuts reads with such a fold down to their longer arm. With paired reads, pairs whose
        fragment folded back between the mates have both mates written as single reads (because they don't
        form a proper pair orientation-wise).

        HOW IT WORKS: With the default of short reads, each read, or each pair, is checked two ways with
        k-mers (default 13-mers). (1) A fold within a read: when there are at least k-mers (--min-shared, default 16)
        whose reverse complement is also in the same read. The read is cut to its longer arm, ending where
        that arm's paired stretch ends, so the join and anything unpaired beside it go; it's then checked again,
        as a read can fold more than once. (2) With pairs, a fold between the mates: read 2 holding at least --min-shared
        of read 1's k-mers in the same orientation. Here the unpaired sequence is trimmed from both, and the remaining
        pieces are sent as individual reads to the output (not the filtered_R1/R2).

        OUTPUTS: Input pairs with neither signature detected go as they are to filtered_R1/R2. Both mates of a
        pair with either found go to foldback-mates as single reads (again, because they aren't proper pairs).
        Cut pieces under --min-length are dropped.

        LONG READS (--hifi): For each read, about 1 in 4 21-mers are used as anchors. It is deterined to have fold
        if 5 or more anchors that appear both forward and reverse-complemented agree on where the read "turns", and those
        pairs comprise at least a quarter of the anchors on the read's shorter side of the turn (which helps delineate
        amplification artifacts from legit inverted repeats). A folded read is cut to its longer arm and checked again.
        --kmer-size and --min-shared are for short reads and aren't used in this mode.
        """

    description = detailed_desc if show_fine else desc
    formatter = ParagraphRichHelpFormatter if show_fine else CustomRichHelpFormatter
    epilog = ("Ex. usage: `bit filter-foldbacks -i R1.fastq.gz -I R2.fastq.gz -o filtered`, "
              "or `bit filter-foldbacks --hifi -i hifi-reads.fastq.gz -o filtered`")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "filter-foldbacks",
            description=description,
            formatter_class=formatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=description,
            epilog=epilog,
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
        help=("Input read 2 file if paired, fastq or fasta, gzipped or not"),
    )

    general.add_argument(
        "--hifi",
        action="store_true",
        help=("Add this flag if the input are pacbio hifi reads (run with -H for detailed help)"),
    )

    general.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="filter-foldbacks",
        help='Directory for output files (default: "filter-foldbacks")',
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
        help="Show detailed help, including fine-tuning parameters and how the program works",
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
        "--min-length",
        metavar="<INT>",
        type=int,
        default=None,
        help=h(f"Minimum length of a read cut at a fold for it to be kept (default: {MIN_LENGTH}; "
               f"{LONG_MIN_LENGTH} with --hifi)"),
    )

    target.add_argument(
        "-k",
        "--kmer-size",
        metavar="<INT>",
        type=int,
        default=KMER_SIZE,
        help=h(f"K-mer size for finding folds in short reads, 5-31 (default: {KMER_SIZE}). Not used "
               "with --hifi"),
    )

    target.add_argument(
        "-s",
        "--min-shared",
        metavar="<INT>",
        type=int,
        default=MIN_SHARED,
        help=h(f"Shared k-mers needed to call a fold in short reads, within a read or between mates "
               f"(default: {MIN_SHARED}). Lower catches more fold-backs but also potentially more, legit short inverted "
               "repeats real genomes carry. Not used with --hifi"),
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

    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater")
    if args.min_length is not None and args.min_length < 1:
        parser.error("--min-length must be 1 or greater")
    if args.hifi and args.input_reads_2:
        parser.error("--hifi takes single long reads: no --input-reads-2")
    if not 5 <= args.kmer_size <= 31:
        parser.error("--kmer-size must be between 5 and 31")
    if args.min_shared < 1:
        parser.error("--min-shared must be 1 or greater")

    from bit.modules.general import check_files_are_found
    from bit.modules.filter_foldbacks import filter_foldbacks, setup_output_dir
    from bit._native import set_threads

    check_files_are_found([args.input_reads] + ([args.input_reads_2] if args.input_reads_2 else []))
    setup_output_dir(args.output_dir, args.output_prefix, args.force_overwrite,
                     reconstruct_invocation(parser, args))
    set_threads(args.jobs)

    counts, paths = filter_foldbacks(
        read_1=args.input_reads,
        read_2=args.input_reads_2,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        min_length=args.min_length,
        min_shared=args.min_shared,
        k=args.kmer_size,
        hifi=args.hifi,
    )

    print_summary(counts, paths, paired=args.input_reads_2 is not None, hifi=args.hifi)


def print_summary(counts, paths, paired, hifi=False):

    unit = "pairs" if paired else "reads"
    n = max(counts.reads_or_pairs_in, 1)
    print()
    print(f"        {counts.reads_or_pairs_in:,} {unit} scanned")
    if paired:
        print(f"        {counts.foldback_pairs:,} fold-back pairs ({100 * counts.foldback_pairs / n:.2f}%): "
              f"{counts.pairs_with_fold_in_a_read:,} with a fold in a read, "
              f"{counts.pairs_with_fold_between_mates:,} folded between the mates")
        print(f"        {counts.pairs_out:,} pairs kept as pairs; {counts.mates_written_as_singles:,} "
              f"mates of fold-back pairs written as single reads")
    else:
        print(f"        {counts.reads_cut_at_a_fold:,} reads cut at a fold "
              f"({100 * counts.reads_cut_at_a_fold / n:.2f}%)")
        q = counts.middle_quartiles() if hifi else None
        if q:
            print(f"        unpaired middle between the arms: median {q[1]:,} bp "
                  f"(quartiles {q[0]:,}-{q[2]:,})")
    if counts.reads_dropped_short:
        print(f"        {counts.reads_dropped_short:,} reads shorter than --min-length after cutting dropped")
    print()
    if paired:
        print(f"    Filtered pairs written to: '{paths['r1']}' and '{paths['r2']}'")
        if paths["singles"] is not None:
            print(f"    Mates of fold-back pairs written to: '{paths['singles']}'")
    else:
        print(f"    Filtered reads written to: '{paths['reads']}'")
    print(f"    Summary written to: '{paths['summary']}'")
    print()
