import sys
import argparse
from bit.modules.cov_analyzer import run_cov_analyzer
from bit.cli.common import CustomRichHelpFormatter, reconstruct_invocation

def main():

    desc = """
        This script identifies and pulls out regions of relatively higher and lower coverage when
        given a reference fasta and a bam file. It generates a bed file of the specified window-
        and step-size, utilizes mosdepth to get the coverage of those windows, then generates stats
        for those windows and pulls out regions with coverage above and below specified thresholds.
        It outputs a table of coverage stats for all windows, a table of merged adjacent windows
        ("regions"), and identified regions in fasta format. By default it looks for coverage variations
        from the global mean coverage, but you can tell it to use per-contig coverages by adding the
        `--per-contig` flag. Additionally, it is recommended to exclude contigs holding mitochondrial
        genomes or chloroplasts due to their generally very relative high coverage, unless you
        specifically want to investigate them too (if so, probably use `--per-contig` mode). For
        version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-cov-analyzer -r reference.fasta -b mapping.bam",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )
    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-r",
        "--reference-fasta",
        metavar="<FILE>",
        help='Reference fasta file',
        required=True,
    )
    required.add_argument(
        "-b",
        "--bam-file",
        metavar="<FILE>",
        help="Bam file",
        required=True,
    )
    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<STR>",
        help='Directory for the output files (default: "cov-analyzer/")',
        default="cov-analyzer",
    )
    optional.add_argument(
        "-H",
        "--high-threshold",
        metavar="<FLOAT>",
        help='High threshold for regions-of-interest (default: 5; e.g., 5 would\
            identify regions with 5-fold+ higher coverage than the mean of the whole ref)',
        type=float,
        default=5,
    )
    optional.add_argument(
        "-L",
        "--low-threshold",
        metavar="<FLOAT>",
        help='Low threshold for regions-of-interest (default: 5; e.g., 5 would\
            identify regions with 5-fold+ lower coverage than the mean of the whole ref)',
        type=float,
        default=5,
    )
    optional.add_argument(
        "-m",
        "--min-region-length",
        metavar="<INT>",
        type=int,
        help='Minimum region length (in bp) to be reported (default: 0; so reports all regions, even single windows)',
    )
    optional.add_argument(
        "-e",
        "--exclude-contigs",
        metavar="<STR>",
        nargs="+",
        help='Space-delimited list of contigs to exclude from analysis (e.g., "-e NC_000932.1 NC_037304.1"; default: None)',
        type=str,
        default=[],
    )
    optional.add_argument(
        "--per-contig",
        action="store_true",
        help='Run the analysis using per-contig mean-coverages for each contig to find high/low regions, instead of using the mean coverage of the whole reference (default: False)',
        default=False,
    )
    optional.add_argument(
        "-s",
        "--sliding-window-size",
        metavar="<INT>",
        help='Sliding window size (default: 50)',
        type=int,
        default=50,
    )
    optional.add_argument(
        "-S",
        "--step-size",
        metavar="<INT>",
        help='Step size for sliding window (default: 10)',
        type=int,
        default=10,
    )
    optional.add_argument(
        "-g",
        "--allowed-gap",
        metavar="<INT>",
        help='Number of bases allowed between qualifying windows (those with coverages above/below the specified thresholds) to still merge them into one contiguous region (default: 1000)',
        type=int,
        default=1000,
    )
    optional.add_argument(
        "-B",
        "--buffer",
        metavar="<INT>",
        help='Add this length to each side of a region of interest when pulled out as fasta (default: 100)',
        type=int,
        default=100,
    )
    optional.add_argument(
        "--no-window-stats",
        help='Add this flag to skip writing out individual window stats (saves spacetime)',
        action="store_true",
    )
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit",
    )

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    # reconstructing the full command-line invocation for logging
    full_cmd_executed = reconstruct_invocation(parser, args)

    run_cov_analyzer(
        reference_fasta = args.reference_fasta,
        bam_file = args.bam_file,
        output_dir = args.output_dir,
        high_threshold = args.high_threshold,
        low_threshold = args.low_threshold,
        min_region_length = args.min_region_length,
        exclude_contigs = args.exclude_contigs,
        per_contig = args.per_contig,
        sliding_window_size = args.sliding_window_size,
        step_size = args.step_size,
        allowed_gap = args.allowed_gap,
        buffer = args.buffer,
        no_window_stats = args.no_window_stats,
        log_file = f"{args.output_dir}/runlog.txt",
        full_cmd_executed = full_cmd_executed
    )
