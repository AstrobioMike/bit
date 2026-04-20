import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.get_mapped_reads_pid import (get_mapped_reads_pids,
                                              get_summary_stats)


def build_parser():

    desc = """
        This script takes an input bam file and generates percent-identity information for mapped reads
        based on edit distance (using the NM field) and total alignment length. By default,
        it just prints out some summary stats. Specify an output file if you also want it to write out the
        percent identities for each mapped read. [bold]TO ALSO GET[/bold] coverage and detection information, use
        `bit-cov-stats` instead. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-mapped-reads-pid input.bam`",
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
        help='Name of output tsv if wanting per-read percent identities written out',
    )

    optional.add_argument(
        "--include-non-primary",
        action="store_true",
        help="Also calculate percent identities for secondary and supplementary (non-primary) alignments",
    )


    add_help(optional)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    store_reads = bool(args.output_tsv)
    ref_read_pids, pid_stats = get_mapped_reads_pids(args.input_bam, include_non_primary=args.include_non_primary, store_read_pids=store_reads)

    print_summary_stats(pid_stats)

    if args.output_tsv:
        write_out_read_pids(args.output_tsv, ref_read_pids)


def print_summary_stats(pid_stats):

    mean, summary_stats = get_summary_stats(pid_stats)

    def fmt(val):
        s = f"{val:.2f}"
        return s.rstrip("0").rstrip(".")

    print(f"\n    Mean percent identity of mapped reads: {fmt(mean)}%\n")

    for name, val in summary_stats:
        value = val if isinstance(val, str) else fmt(val)

        print(f"        {name:<{20}}{value}")

    print()


def write_out_read_pids(output_tsv, ref_read_pids):

    with open(output_tsv, "w") as output:
        output.write("contig\tread_ID\tpercent_identity\n")
        for refname, data_list in ref_read_pids.items():
            for read_id, pid in data_list:
                output.write(f"{refname}\t{read_id}\t{pid:.2f}\n")

    print(f"    Per-read percent identities were written to: '{output_tsv}'\n")
