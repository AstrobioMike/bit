import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg,
                            reconstruct_invocation)


def build_parser(parent_subparsers=None):

    desc = """
        This program removes fold-back artifacts of whole-genome amplification (e.g., PTA or MDA)
        from short reads (the default) or long reads (`--hifi`). An amplified molecule that switched to the opposite strand nearby reads as
        a stretch and then its own reverse complement. A read with such a fold in it is cut to its
        longer arm, up to where that arm's paired stretch ends (dropping the join and any unpaired,
        possibly untemplated, bases beside it). With paired reads, a pair whose fragment folded back
        between the mates (read 2 holds read 1's sequence in the same orientation) is also caught.
        Pairs with neither are written as pairs; both mates of a fold-back pair are written as single
        reads, since each is still real sequence though the pairing isn't. With `--hifi` (long reads,
        e.g., PacBio HiFi), each read is cut at its folds the same way `shifia --pta` does: a long
        read holds whole fold-backs, so no pairs are needed. Works from the reads alone (no reference
        needed).
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "filter-foldbacks",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog=("Ex. usage: `bit filter-foldbacks -i R1.fastq.gz -I R2.fastq.gz -o filtered`, "
                    "or `bit filter-foldbacks --hifi -i hifi-reads.fastq.gz -o filtered`"),
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
        help=("Input read 2 file if paired, fastq or fasta, gzipped or not (paired reads catch "
              "more fold-backs: those between the mates as well as those within a read)"),
    )

    optional.add_argument(
        "--hifi",
        action="store_true",
        help=("Input reads are long single reads (e.g., PacBio HiFi): each is cut at its folds as "
              "`shifia --pta` does (21-mers, a fold's arms pairing at least a quarter of the read's "
              "shorter side, so inverted repeats in long flanks are left alone); --min-shared isn't "
              "used (default: short reads)"),
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="filter-foldbacks",
        help='Directory for output files (default: "filter-foldbacks")',
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
        "--min-length",
        metavar="<INT>",
        type=int,
        default=None,
        help="Minimum length of a read cut at a fold for it to be kept (default: 50; 100 with --hifi)",
    )

    optional.add_argument(
        "--min-shared",
        metavar="<INT>",
        type=int,
        default=20,
        help=("Shared k-mers (17-mers) needed to call a fold in short reads, within a read or "
              "between mates; 20 means at least 36 bp shared as one stretch, more than most short "
              "inverted repeats in genomes (default: 20)"),
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

    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater")
    if args.min_length is not None and args.min_length < 1:
        parser.error("--min-length must be 1 or greater")
    if args.hifi and args.input_reads_2:
        parser.error("--hifi takes single long reads: no --input-reads-2")
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
        print(f"    Mates of fold-back pairs written to: '{paths['singles']}'")
    else:
        print(f"    Filtered reads written to: '{paths['reads']}'")
    print(f"    Summary written to: '{paths['summary']}'")
    print()
