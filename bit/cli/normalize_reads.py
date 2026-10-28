import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            ParagraphRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg,
                            reconstruct_invocation)
from bit.modules.normalize_reads import (TARGET_DEPTH,
                                   KMER_SIZE,
                                   DENSITY,
                                   NOVEL_FRACTION,
                                   MIN_READS,
                                   ERROR_WEIGHT,
                                   ERROR_SHARE)


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        This program normalizes read depth (e.g., for assembly) while keeping reads that add sequence
        not yet covered, including from low-abundance members and minority variants. Run with `-H` to
        see the detailed help menu, including how the program works and fine-tuning parameters.
        """

    detailed_desc = desc + """

        HOW IT WORKS: Each read, or each pair (both mates together), is reduced to sampled k-mers (by
        default 25-mers, about 1 in 2 of them). Programs like bbnorm keep or drop a read by its median
        k-mer depth, so a read from a rare genome or a minority variant that is mostly shared, deep
        sequence is thinned like the rest. Here a read is judged by its least-covered part instead.

        The reads are read twice. The first pass counts how many reads hold each k-mer (in a fixed
        amount of memory, --memory). The second judges and writes: k-mers in fewer than --min-reads
        reads, or far rarer than the k-mers around them, are taken for sequencing errors and count for
        nothing. A read's depth is that of its 10% least-covered k-mers (--novel-fraction), and it's
        kept with probability --target-depth / that depth (always, when that's at or under the
        target). Reads with a few k-mers on a minor branch (k-mers in 3+ reads, at 5-50% of the read's
        typical depth: a minority variant) are kept up to half the target for that branch
        (--no-minor-branches turns this off). Which reads are kept depends only on the reads
        themselves, not on their order, and is the same run to run.

        OUTPUTS: normalized_R1/R2 (or normalized, if single-end), as the reads were in the input,
        gzipped; summary.tsv; and with --depth-table, depths.tsv.gz, each kept read's typical
        k-mer depth in the whole input (e.g., for abundances).
        """

    description = detailed_desc if show_fine else desc
    formatter = ParagraphRichHelpFormatter if show_fine else CustomRichHelpFormatter
    epilog = "Ex. usage: `bit normalize-reads -i R1.fastq.gz -I R2.fastq.gz -o normalized`"

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "normalize-reads",
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
        help="Input read 2 file if paired, fastq or fasta, gzipped or not",
    )

    general.add_argument(
        "-d",
        "--target-depth",
        metavar="<INT>",
        type=int,
        default=TARGET_DEPTH,
        help=f"Depth to cap reads at (default: {TARGET_DEPTH})",
    )

    general.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="normalize-reads",
        help='Directory for output files (default: "normalize-reads")',
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

    general.add_argument(
        "-m",
        "--memory",
        metavar="<FLOAT>",
        type=float,
        default=None,
        help=("GB of memory for counting k-mers (default: about the input's size on disk, "
              "0.25-8, at most half the computer's memory)"),
    )

    general.add_argument(
        "--depth-table",
        action="store_true",
        help="Also write each kept read's k-mer depth in the whole input (depths.tsv.gz)",
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
        fine = parser.add_argument_group("Fine-tuning Parameters")
        _add_fine_tuning_arguments(fine, hidden=False)
    else:
        _add_fine_tuning_arguments(parser, hidden=True)

    return parser


def _add_fine_tuning_arguments(target, hidden=False):

    def h(text):
        return argparse.SUPPRESS if hidden else text

    target.add_argument(
        "-k",
        "--kmer-size",
        metavar="<INT>",
        type=int,
        default=KMER_SIZE,
        help=h(f"K-mer size, 15-32 (default: {KMER_SIZE})"),
    )

    target.add_argument(
        "--density",
        metavar="<INT>",
        type=int,
        default=DENSITY,
        help=h(f"About 1 in this many k-mers are used (default: {DENSITY}; 1 uses every k-mer)"),
    )

    target.add_argument(
        "--novel-fraction",
        metavar="<FLOAT>",
        type=float,
        default=NOVEL_FRACTION,
        help=h(f"Share of a read's (or pair's) least-covered k-mers its depth is judged by "
               f"(default: {NOVEL_FRACTION})"),
    )

    target.add_argument(
        "--min-reads",
        metavar="<INT>",
        type=int,
        default=MIN_READS,
        help=h(f"Reads a k-mer needs to be in to count as sequence rather than an error "
               f"(default: {MIN_READS})"),
    )

    target.add_argument(
        "--error-weight",
        metavar="<FLOAT>",
        type=float,
        default=ERROR_WEIGHT,
        help=h(f"Where reads are thinned, keep ones holding a k-mer judged an error at this "
               f"times the usual rate, and error-free ones more often to make up for it "
               f"(0-1; default: {ERROR_WEIGHT:g}, which treats them alike)"),
    )

    target.add_argument(
        "--error-share",
        metavar="<FLOAT>",
        type=float,
        default=ERROR_SHARE,
        help=h(f"K-mers in under this share of the reads around them are taken for sequencing "
               f"errors (0-1; default: {ERROR_SHARE:g}). Raising it drops recurrent errors at "
               f"very high depth, and with them variants under that share"),
    )

    target.add_argument(
        "--strand-aware",
        action="store_true",
        help=h("Also count k-mers by the strand reads hold them on (50%% more memory), and take "
               "well-covered k-mers seen almost only on one strand, well under the depth around "
               "them, for systematic sequencing errors"),
    )

    target.add_argument(
        "--no-minor-branches",
        action="store_true",
        help=h("Don't keep reads for k-mers of a minor branch (minority variants)"),
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
    if args.target_depth < 1:
        parser.error("--target-depth must be 1 or greater")
    if not 15 <= args.kmer_size <= 32:
        parser.error("--kmer-size must be between 15 and 32")
    if args.density < 1:
        parser.error("--density must be 1 or greater")
    if not 0 < args.novel_fraction <= 1:
        parser.error("--novel-fraction must be greater than 0 and at most 1")
    if args.min_reads < 1:
        parser.error("--min-reads must be 1 or greater")
    if not 0 < args.error_share < 1:
        parser.error("--error-share must be between 0 and 1")
    if not 0 <= args.error_weight <= 1:
        parser.error("--error-weight must be between 0 and 1")
    if args.memory is not None and args.memory <= 0:
        parser.error("--memory must be greater than 0")

    from bit.modules.general import check_files_are_found
    from bit.modules.normalize_reads import normalize, setup_output_dir
    from bit._native import set_threads

    check_files_are_found([args.input_reads] + ([args.input_reads_2] if args.input_reads_2 else []))
    setup_output_dir(args.output_dir, args.output_prefix, args.force_overwrite,
                     reconstruct_invocation(parser, args))
    set_threads(args.jobs)

    counts, paths = normalize(
        read_1=args.input_reads,
        read_2=args.input_reads_2,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        target_depth=args.target_depth,
        k=args.kmer_size,
        density=args.density,
        novel_fraction=args.novel_fraction,
        minor_branches=not args.no_minor_branches,
        min_reads=args.min_reads,
        sketch_bytes=None if args.memory is None else int(args.memory * 2**30),
        error_weight=args.error_weight,
        error_share=args.error_share,
        strand_aware=args.strand_aware,
        depth_table=args.depth_table,
    )

    print_summary(counts, paths, paired=args.input_reads_2 is not None)


def print_summary(counts, paths, paired):

    unit = "pairs" if paired else "reads"
    n = max(counts.reads_or_pairs_in, 1)
    b = max(counts.bases_in, 1)
    print()
    print(f"        {counts.reads_or_pairs_in:,} {unit} in ({counts.bases_in:,} bases)")
    print(f"        {counts.reads_or_pairs_out:,} {unit} kept ({100 * counts.reads_or_pairs_out / n:.2f}%; "
          f"{100 * counts.bases_out / b:.2f}% of bases)")
    if counts.no_anchors:
        print(f"        {counts.no_anchors:,} {unit} with no usable k-mers (too short, or Ns) dropped")
    print(f"        median k-mer depth {counts.median_depth_in:,} in the input, about "
          f"{counts.median_depth_out:,} in what was kept")
    print(f"        {100 * counts.error_share:.2f}% of k-mers judged sequencing errors")
    print()
    if paired:
        print(f"    Normalized pairs written to: '{paths['r1']}' and '{paths['r2']}'")
    else:
        print(f"    Normalized reads written to: '{paths['reads']}'")
    if "depth_table" in paths:
        print(f"    Depths of kept reads written to: '{paths['depth_table']}'")
    print(f"    Summary written to: '{paths['summary']}'")
    if counts.sketch_too_full:
        from bit.modules.general import report_message
        report_message(f"The k-mer counting ended up {100 * counts.sketch_fill:.0f}% full, so counts of rare "
                       "k-mers can run high, which can drop reads that should have been kept. Giving it "
                       "more memory (--memory) would help.", initial_indent="    ", subsequent_indent="    ")
    print()
