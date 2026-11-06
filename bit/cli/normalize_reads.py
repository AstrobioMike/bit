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
        This program normalizes read depth (e.g., prior to assembly) and is, with default settings,
        geared towards retaining low-abundance members and minority variants. Run with `-H` to see
        the detailed help menu, including how the program works and fine-tuning parameters.
        """

    detailed_desc = desc + """

        HOW IT WORKS: One of my favorites, bbnorm, keeps or drops a read roughly based on its median
        k-mer depth. This means, by design, a read that holds some sequence from a rare genome or a
        minority variant alongside a majority of shared sequence with a more abundantly represented member of
        the community will likely be dropped (losing that rare sequence component). Here a read is judged
        by its 10th percentile k-mer depth instead, in an effort to retain more low-abundance variants
        if they are present. This helps to recover low-abundance variants that may otherwise be lost, but it
        can also lead to slightly more fragmented assemblies in some situations if the assembler is unable to
        disentangle the low-abundance variants from the more abundant sequence. The `--novel-fraction` parameter
        can be raised from the default to 0.5 to behave more like bbnorm while benefiting from the speed-up
        offered by this program.

        The reads are read twice, with standard Illumina adapters, adapter read-through (if paired-end),
        and poly-G tails trimmed off in both passes (--no-adapter-trimming turns this off). The first
        pass counts how many reads hold each k-mer, counting exact duplicates only once, so the errors
        they share don't pass for sequence seen in two reads. The second judges and writes: k-mers in
        fewer than --min-reads reads, or far rarer than the k-mers around them, are judged as sequencing
        errors and therefore don't contribute to the read's depth, as are well-covered k-mers seen almost
        only on one strand and well under the depth around them (likely systematic errors).
        A read's depth is that of its 10% least-covered k-mers (--novel-fraction), and it's kept
        with probability --target-depth / that depth. Where reads are thinned, reads holding an error are
        kept less often than the rest (--error-weight), and only one copy of an exact duplicate is written
        (reads at or under the target keep all their duplicate copies, which can help assemblers hold on to
        low-coverage sequence). Reads with a few k-mers on a minor branch (k-mers in 3+ reads, at 5-50% of the read's
        trusted k-mer depth: a minority variant) are kept up to half the target for that branch
        (--no-minor-branches turns this off).

        OUTPUTS: normalized_R1/R2 (or normalized, if single-end); summary.tsv; and with --depth-table,
        depths.tsv.gz, each kept read's typical k-mer depth in the whole input.
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
        help=f"K-mer depth to normalize to; the read depth that comes out is somewhat higher "
             f"(default: {TARGET_DEPTH})",
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
        help=("GB of memory for counting k-mers, plus half again for counting them by strand "
              "(default: 0.5 per GB of gzipped input, at least 0.25, and with the strand counts "
              "at most half the memory free)"),
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
        help=h(f"K-mer size, odd, 15-31 (default: {KMER_SIZE})"),
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
        help=h(f"Reads a k-mer needs to be in to be potentially trusted rather than assumed to be an error "
               f"(default: {MIN_READS}; 'potentially' because --error-share is at play also)"),
    )

    target.add_argument(
        "--error-share",
        metavar="<FLOAT>",
        type=float,
        default=ERROR_SHARE,
        help=h(f"A k-mer seen under this proportion as often as its neighboring k-mers is judged an error "
               f"(0-1; default: {ERROR_SHARE:g}). Raising it can help drop recurrent errors at "
               f"very high depth, at the cost of potentially losing low-abundance variants"),
    )

    target.add_argument(
        "--error-weight",
        metavar="<FLOAT>",
        type=float,
        default=ERROR_WEIGHT,
        help=h(f"Where reads are thinned, keep ones holding a k-mer judged an error at this "
               f"times the usual rate, thereby attempting to favor error-free reads being retained "
               f"(0-1, 1 meaning keep them at the same rate; default: {ERROR_WEIGHT:g})"),
    )

    target.add_argument(
        "--no-minor-branches",
        action="store_true",
        help=h("Don't keep reads for k-mers of a minor branch (minority variants)"),
    )

    target.add_argument(
        "--no-adapter-trimming",
        action="store_true",
        help=h("Add this flag to skip trimming adapters (built-in Illumina adapters, and past the insert where "
               "mates overlap) and poly-G tails. Trimming is fast and best left on: untrimmed adapters can "
               "wreak havoc on assembly"),
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
    if not (15 <= args.kmer_size <= 31 and args.kmer_size % 2 == 1):
        parser.error("--kmer-size must be odd and between 15 and 31")
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
        trim_adapters=not args.no_adapter_trimming,
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
    if counts.reads_adapter_trimmed or counts.reads_poly_g_trimmed:
        print(f"        {counts.reads_adapter_trimmed:,} reads trimmed of adapter read-through, "
              f"{counts.reads_poly_g_trimmed:,} of poly-G tails ({counts.bases_trimmed:,} bases)")
    if counts.too_short_after_trimming:
        print(f"        {counts.too_short_after_trimming:,} {unit} too short after trimming dropped")
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
        gb = counts.suggested_memory_gb
        at_least = "at least " if counts.sketch_fill > 0.99 else ""
        report_message(f"The k-mer counting ended up {100 * counts.sketch_fill:.0f}% full, so counts of rare "
                       "k-mers can run high: errors can pass as sequence (keeping reads that needn't be "
                       "kept), and rare sequence can look deeper than it is (dropping reads that should "
                       f"have been kept). Rerunning with {at_least}--memory {gb:g} would bring it to about "
                       f"half full (using about {1.5 * gb:g} GB, with the counts by strand).",
                       initial_indent="    ", subsequent_indent="    ")
    print()
