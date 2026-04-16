import sys
import argparse
from bit.modules.get_cov_stats import get_cov_stats
from bit.cli.common import CustomRichHelpFormatter, add_force, add_help
from bit.modules.general import report_message, notify_premature_exit, is_gzipped


def build_parser():

    desc = """
        This script generates whole-reference and contig-level detection and coverage info
        given the input reference fasta(s) and a bam file AND/OR a mosdepth-produced per-base.bed.gz file.
        If providing a bam file, it will also report mean and median percent ID values of mapped reads to
        each input reference and contig. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-get-cov-stats -r reference.fasta -b mapping.bam` or `bit-get-cov-stats -r reference.fasta --bed per-base.bed.gz`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters (choose one of `--reference-fastas` or `--reference-list` then `--bam` and/or `--bed`)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-r",
        "--reference-fastas",
        metavar="<STR>",
        help='Path to reference fasta file(s) (this takes precedence over -R/[cyan]--reference-list[/] if both provided) OR',
        nargs="+",
    )

    required.add_argument(
        "-R",
        "--reference-list",
        metavar="<FILE>",
        help="Path to a file containing reference fasta paths, one per line.\n"
             "This is an alternative to -r/[cyan]--reference-fastas[/].",
    )

    required.add_argument(
        "-b",
        "--bam",
        metavar="<PATH>",
        help="Path to bam file (which will allow reporting of mean percent ID per input reference) AND/OR",
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
        help='Name of the output prefix (default: "cov-stats")',
        default="cov-stats",
    )

    optional.add_argument(
        "--skip-per-contig",
        action="store_true",
        help="Specify this to skip reporting stats on a per-contig basis (which might save a lot of spacetime depending on the inputs)",
    )

    optional.add_argument(
        "--include-non-primary",
        action="store_true",
        help="Include secondary and supplementary alignments in the percent ID calculations (if being done), and run mosdepth with '--flag 1540' set",
    )

    optional.add_argument(
        "--skip-read-pids",
        action="store_true",
        help="Skip calculating mean and median percent ID values of mapped reads to each input reference and contig (can save spacetime)",
    )

    add_force(optional)
    add_help(optional)

    return parser


def main():
        parser = build_parser()

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(0)

        args = parser.parse_args()

        check_required_inputs(args)

        args = set_ref_input(args)

        get_cov_stats(args)


def set_ref_input(args):
    if args.reference_fastas and args.reference_list:
        report_message("You have provided both -r/--reference-fastas and -R/--reference-list parameters, please only provide one.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if getattr(args, "reference_list", None):
        ref_file = args.reference_list
        try:
            with open(ref_file, "r") as fh:
                lines = [l.strip() for l in fh if l.strip()]
            args.reference_fastas = lines
        except Exception:
            report_message(f"We were unable to read the speficied file: {ref_file}")
            notify_premature_exit()
    return args


def check_required_inputs(args):
    if not args.reference_fastas and not args.reference_list:
        report_message("You must provide at least one input reference fasta file via the -r/--reference-fastas or -R/--reference-list parameters.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if not args.bam and not args.bed:
        report_message("You must provide either a bam file via the -b/--bam parameter or "
                       "a mosdepth-produced per-base.bed.gz file via the --bed parameter.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
    if args.bed and not is_gzipped(args.bed):
        report_message("The provided bed file does not appear to be gzipped. Please provide "
                       "a mosdepth-produced per-base.bed.gz file that is gzipped.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()
