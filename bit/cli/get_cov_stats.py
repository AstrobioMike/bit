import sys
import argparse
from bit.modules.get_cov_stats import get_cov_stats
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import report_message, notify_premature_exit, is_gzipped


def build_parser():

    desc = """
        This script generates whole-reference detection and coverage info
        for specified references given the input reference fasta(s) and either a bam file OR
        a mosdepth-produced per-base.bed.gz file. If provided a bam file, it will generate
        the required mosdepth files as well as also adding mean percent ID of mapped reads
        to each input reference in the output table. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-cov-stats -r reference.fasta -b mapping.bam` or `bit-get-cov-stats -r reference.fasta --bed per-base.bed.gz`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS (choose one of bam or bed input in addition to refs)")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-r",
        "--reference-fastas",
        metavar="<STR>",
        help='Path to reference fasta file(s)',
        nargs="+",
    )

    required.add_argument(
        "-b",
        "--bam",
        metavar="<PATH>",
        help="Path to bam file (which will allow reporting of mean percent ID per input reference) OR",
    )

    required.add_argument(
        "--bed",
        metavar="<PATH>",
        help="Path to mosdepth-produced per-base.bed.gz file (if already available)",
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help='Name of the output prefix (default: "coverage-stats")',
        default="coverage-stats",
    )

    add_help(optional)

    return parser


def main():
        parser = build_parser()

        if len(sys.argv) == 1:  # pragma: no cover
            parser.print_help(sys.stderr)
            sys.exit(0)

        args = parser.parse_args()

        check_required_inputs(args)

        get_cov_stats(args)


def check_required_inputs(args):
    if not args.reference_fastas:
        report_message("You must provide at least one reference fasta file via the -r/--reference-fastas parameter.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if not args.bam and not args.bed:
        report_message("You must provide either a bam file via the -b/--bam parameter or "
                       "a mosdepth-produced per-base.bed.gz file via the --bed parameter.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if args.bam and args.bed:
        report_message("You have provided both a bam file and a bed file. Please provide "
                       "solely the bed file IF it is a mosdepth-produced per-base.bed.gz file. "
                       "If you don't have that, then provide solely the bam file.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if args.bed and not is_gzipped(args.bed):
        report_message("The provided bed file does not appear to be gzipped. Please provide "
                       "a mosdepth-produced per-base.bed.gz file that is gzipped.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
