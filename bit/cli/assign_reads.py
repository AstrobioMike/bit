import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_force,
                            add_version_arg,
                            reconstruct_invocation)


def build_parser(parent_subparsers=None):

    desc = """
        This program assigns reads to the reference sequences they exactly match (or nearly exactly match,
        within `--max-edits`). It can be useful for finding the origins of reads from known references
        where mapping-based approaches might struggle with multi-mapping and low MAPQs that need to be annoyingly parsed.
        Note the `--per-seq` and `--circular` options. Reads assigned to exactly one input reference are reported
        as "unique", and reads tied between more than one as "ambiguous".
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
            epilog="Ex. usage: `bit assign-reads -r ref-1.fasta ref-2.fasta -i reads.fastq.gz --circular`",
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
        help=("Reference fasta file(s). By default, each file is one reference made up of all of its "
              "sequences, named by its file name (see --per-seq)"),
    )

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
        default="assign-reads",
        help='Directory for output files (default: "assign-reads")',
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
        "--per-seq",
        action="store_true",
        help=("Treat each sequence in the input reference fasta(s) as its own reference, named by its "
              "sequence name, rather than each file being one reference"),
    )

    optional.add_argument(
        "-c",
        "--circular",
        action="store_true",
        help="Treat reference sequences as circular, so reads spanning the origin are matched",
    )

    optional.add_argument(
        "-e",
        "--max-edits",
        metavar="<INT>",
        type=int,
        default=0,
        help=("Maximum edit distance (substitutions + indels) allowed between a read and a reference "
              "(default: 0, exact matches only, higher values will slow things down)"),
    )

    optional.add_argument(
        "--min-frac-of-seq",
        metavar="<FLOAT>",
        type=float,
        default=0.0,
        help=("Only count a match if the read is at least this fraction of the length of the "
              "sequence it matches (applies to single-end input only; default: 0)"),
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
    if args.max_edits < 0:
        parser.error("--max-edits must be 0 or greater")
    if not 0 <= args.min_frac_of_seq <= 1:
        parser.error("--min-frac-of-seq must be between 0 and 1")
    if args.input_reads_2:
        parser.error("--min-frac-of-seq only applies to single-end input")

    from bit.modules.general import check_files_are_found
    from bit.modules.assign_reads import assign_reads, setup_output_dir

    check_files_are_found(args.refs + [args.input_reads] + ([args.input_reads_2] if args.input_reads_2 else []))
    setup_output_dir(args.output_dir, args.output_prefix, args.force_overwrite,
                     reconstruct_invocation(parser, args))

    summary = assign_reads(
        ref_paths=args.refs,
        read_1=args.input_reads,
        read_2=args.input_reads_2,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        per_seq=args.per_seq,
        circular=args.circular,
        max_edits=args.max_edits,
        min_frac_of_seq=args.min_frac_of_seq,
        write_reads=args.write_reads,
        jobs=args.jobs,
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
    if summary["wrote_seq_summary"]:
        print(f"    Per-sequence summary written to: '{paths['seq_summary']}'")
    if summary["write_reads"]:
        print(f"    Uniquely assigned reads written to: '{paths['reads_dir']}/'")
    # print(f"    Command info written to: '{paths['log']}'")
    print()
