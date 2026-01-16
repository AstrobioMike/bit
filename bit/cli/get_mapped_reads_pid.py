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
        percent identities for each mapped read. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-mapped-reads-pid input.bam`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

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

    add_help(optional)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv)==1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    ref_read_pids, all_pids = get_mapped_reads_pids(args.input_bam)

    print_summary_stats(all_pids)

    if args.output_tsv:
        write_out_read_pids(args.output_tsv, ref_read_pids)


def print_summary_stats(all_pids):

    mean, summary_stats = get_summary_stats(all_pids)

    def fmt(val):
        s = f"{val:.2f}"
        return s.rstrip("0").rstrip(".")

    print(f"\n    Mean percent identity of mapped reads: {fmt(mean)}%\n")

    for name, val in summary_stats:
        value = str(val) if isinstance(val, int) else fmt(val)
        print(f"        {name:<{20}}{value}")

    print()


def write_out_read_pids(output_tsv, ref_read_pids):

    with open(output_tsv, "w") as output:
        output.write("contig\tread_ID\tpercent_identity\n")
        for refname, data_list in ref_read_pids.items():
            for read_id, pid in data_list:
                output.write(f"{refname}\t{read_id}\t{pid:.2f}\n")

    print(f"    Per-read percent identities were written to: '{output_tsv}'\n")
