import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg)


def build_parser(parent_subparsers=None):

    desc = """
        This program assigns reads to the reference sequences they exactly match or nearly exactly match
        (within `--max-edits`). It can be useful for sorting reads among highly similar references,
        like a mix of near-identical plasmids or constructs, where mapping-based approaches might struggle
        with multimapping and low MAPQs that need to be annoying parsed. Add `--circular` for circular references
        like plasmids so reads spanning the origin are matched. Reads assigned to exactly one input reference
        are reported as "unique", and reads tied between more than one as "ambiguous".
        For paired-end input, a pair is assigned to references that both mates match (ranked by the mates'
        combined edit distance when edits are allowed).
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "assign-reads",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit assign-reads -r ref-1.fasta ref-2.fasta -1 reads.fastq.gz --circular`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-r",
        "--refs",
        metavar="<FILE(s)>",
        nargs="+",
        required=True,
        help="Reference fasta file(s); every sequence is treated as a separate reference",
    )

    required.add_argument(
        "-1",
        "--read-1",
        metavar="<FILE>",
        required=True,
        help="Input reads (read 1 if paired), fastq or fasta, gzipped or not",
    )

    optional.add_argument(
        "-2",
        "--read-2",
        metavar="<FILE>",
        help="Input read 2 file if paired",
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        default="assign-reads",
        help='Output-file prefix (default: "assign-reads")',
    )

    optional.add_argument(
        "-c",
        "--circular",
        action="store_true",
        help="Treat references as circular, so reads spanning the origin are matched",
    )

    optional.add_argument(
        "-e",
        "--max-edits",
        metavar="<INT>",
        type=int,
        default=0,
        help=("Maximum edit distance (substitutions + indels) allowed between a read and a reference "
              "(default: 0, exact matches only)"),
    )

    optional.add_argument(
        "--min-frac-of-ref",
        metavar="<FLOAT>",
        type=float,
        default=0.0,
        help=("Only count a match if the read is at least this fraction of the reference's length "
              "(applies to single-end input only; default: 0)"),
    )

    optional.add_argument(
        "--min-read-len",
        metavar="<INT>",
        type=int,
        default=0,
        help="Skip reads shorter than this (for paired-end, applies to each mate; default: 0)",
    )

    optional.add_argument(
        "--write-reads",
        action="store_true",
        help=("Write reads uniquely assigned to each reference to <prefix>-reads/ "
              "(one file per reference, or an R1/R2 pair of files if paired-end)"),
    )

    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=5,
        help="Number of parallel processes for assigning reads (default: 5)",
    )

    optional.add_argument(
        "-k",
        "--kmer-size",
        metavar="<INT>",
        type=int,
        default=31,
        help="k-mer size for the prefilter; must be <= read lengths for the prefilter to apply (default: 31)",
    )

    optional.add_argument(
        "-s",
        "--num-kmer-samples",
        metavar="<INT>",
        type=int,
        default=25,
        help=("Number of k-mers sampled across each read for the exact-match prefilter. More samples "
              "exclude non-matching references more often. With --max-edits, all non-overlapping k-mers "
              "of each read are used instead (default: 25)"),
    )

    optional.add_argument(
        "--batch-size",
        metavar="<INT>",
        type=int,
        default=500,
        help="Number of reads (or pairs) sent to each process at a time (default: 500)",
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
    if args.kmer_size < 1:
        parser.error("--kmer-size must be 1 or greater")
    if args.num_kmer_samples < 1:
        parser.error("--num-kmer-samples must be 1 or greater")
    if args.max_edits < 0:
        parser.error("--max-edits must be 0 or greater")
    if args.batch_size < 1:
        parser.error("--batch-size must be 1 or greater")
    if not 0 <= args.min_frac_of_ref <= 1:
        parser.error("--min-frac-of-ref must be between 0 and 1")
    if args.read_2 and args.min_frac_of_ref > 0:
        parser.error("--min-frac-of-ref only applies to single-end input")

    from bit.modules.general import check_files_are_found
    from bit.modules.assign_reads import assign_reads, check_outputs

    check_files_are_found(args.refs + [args.read_1] + ([args.read_2] if args.read_2 else []))
    check_outputs(args.output_prefix, args.write_reads, args.force_overwrite)

    summary = assign_reads(
        ref_paths=args.refs,
        read_1=args.read_1,
        read_2=args.read_2,
        output_prefix=args.output_prefix,
        circular=args.circular,
        max_edits=args.max_edits,
        k=args.kmer_size,
        num_kmer_samples=args.num_kmer_samples,
        min_read_len=args.min_read_len,
        min_frac_of_ref=args.min_frac_of_ref,
        write_reads=args.write_reads,
        jobs=args.jobs,
        batch_size=args.batch_size,
    )

    print_summary(summary)


def print_summary(summary):

    total = summary["total"]
    unit = summary["unit"]
    pct = lambda x: f"{100 * x / total:.2f}%" if total else "NA"

    print()
    print(f"        {'Total ' + unit + ':':<34}{total:,}")
    print(f"        {'Uniquely assigned:':<34}{summary['unique']:,} ({pct(summary['unique'])})")
    print(f"        {'Ambiguous (multiple refs):':<34}{summary['ambiguous']:,} ({pct(summary['ambiguous'])})")
    print()

    paths = summary["paths"]
    print(f"    Per-{unit[:-1]} assignments written to: '{paths['hits']}'")
    print(f"    Per-reference summary written to: '{paths['summary']}'")
    if summary["write_reads"]:
        print(f"    Uniquely assigned reads written to: '{paths['reads_dir']}/'")
    print()
