import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from bit.modules.mapped_read_stats import (get_mapped_reads_pids,
                                              get_summary_stats)


def build_parser(parent_subparsers=None):

    desc = """
        This program takes an input bam file and generates percent-identity and clipping information for
        mapped reads. By default, it just prints out some summary stats. Specify an output file if you also want it
        to write out per-read info. To get coverage and detection information, use `bit cov-stats` instead.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "mapped-read-stats",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit mapped-read-stats input.bam`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "input_bam",
        help="Input bam file",
        metavar="<FILE>",
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv if wanting per-read info written out',
    )

    optional.add_argument(
        "--include-non-primary",
        action="store_true",
        help=("Also include secondary and supplementary (non-primary) alignments in the "
              "percent-identity stats. Read-length and clipping stats stay primary-only. "
              "With this set, the percent-identity stats describe alignments rather than "
              "reads (compare 'Num alignments' against 'Num mapped reads' in the summary)."),
    )

    add_help(optional)

    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    store_reads = bool(args.output_tsv)

    ref_read_pids, pid_stats, pid_gc_stats, clip_agg = get_mapped_reads_pids(
        args.input_bam, include_non_primary=args.include_non_primary,
        store_read_pids=store_reads)

    print_summary_stats(pid_stats, pid_gc_stats, clip_agg, args.include_non_primary)

    if args.output_tsv:
        write_out_read_stats(args.output_tsv, ref_read_pids)


def print_summary_stats(pid_stats, pid_gc_stats=None, clip_agg=None, include_non_primary=False):

    summary_stats = get_summary_stats(pid_stats, pid_gc_stats, clip_agg, include_non_primary)

    def fmt(val):
        s = f"{val:.2f}"
        return s.rstrip("0").rstrip(".")

    print()

    for name, val in summary_stats:
        value = val if isinstance(val, str) else fmt(val)
        print(f"        {name:<{28}}{value}")

    print()


def write_out_read_stats(output_tsv, ref_read_pids):

    with open(output_tsv, "w") as output:

        output.write("contig\tread_ID\tread_length\tread_aligned_bases\tref_start\t"
                     "ref_end\tsoft_clipped\thard_clipped\tclipped_percent\t"
                     "percent_identity_gap_aware\tpercent_identity_gap_compressed\n")

        for refname, data_list in ref_read_pids.items():
            for (read_id, pid, gc_pid, aln, total_len, soft, hard, frac,
                 ref_start, ref_end) in data_list:
                output.write(f"{refname}\t{read_id}\t{total_len}\t{aln}\t{ref_start}\t"
                             f"{ref_end}\t{soft}\t{hard}\t{frac * 100:.2f}\t"
                             f"{pid:.2f}\t{gc_pid:.2f}\n")

    print(f"    Per-read stats were written to: '{output_tsv}'\n")
