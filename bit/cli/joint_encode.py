import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_seed
from bit.modules.general import check_files_are_found, report_message, notify_premature_exit

def build_parser():

    desc = """
        This program trims and joint encodes a 3di alignment and an amino-acid alignment (that is derived from the 3di alignment).
        It removes columns based on a gap threshold, and it swaps amino-acid characters into the base 3di alignment
        based on specified column-variability criteria. It produced the joint-encoded alignment and a partitions file.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-joint-encode -i 3di-alignment.fasta -a AA-alignment.fasta -o output-alignment.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group('Required Parameters')
    optional = parser.add_argument_group('Optional Parameters')

    required.add_argument(
        "-i",
        "--input-3di-alignment",
        metavar="<FILE>",
        help="Input fasta file containing 3di alignment",
        required=True,
    )
    required.add_argument(
        "-a",
        "--input-aa-alignment",
        metavar="<FILE>",
        help="Input fasta file containing amino-acid alignment (must be derived from the 3di alignment)",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar="<STR>",
        help="Output prefix (default: 'joint-encoded')",
        default="joint-encoded",
    )
    optional.add_argument(
        "-g",
        "--gap-threshold",
        metavar="<FLOAT>",
        type=float,
        help="Threshold for removing gap columns, e.g. 0.5 says if >= 50%% of the positions are gaps, that column will be removed (default: 0.5)",
        default=0.5,
    )
    optional.add_argument(
        "-s",
        "--swap-3di-variability-threshold",
        metavar="<FLOAT>",
        type=float,
        help="Columns in the 3di alignment with variability <= this value will have their 3di characters swapped to AA characters IF the AA column variability is higher (default: 0.10)",
        default=0.20,
    )
    optional.add_argument(
        "-t",
        "--trim-variability-upper-percentile",
        metavar="<FLOAT>",
        type=float,
        help="""
            Columns will be trimmed from the alignment if their variability in 3di space is above this percentile, OR
            if variability in 3di space is below the `-s` threshold AND their variability in AA space is above this percentile.
            (default: 95)
            """,
        default=95,
    )
    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    preflight_checks(args)

    from bit.modules.joint_encode import joint_encode

    joint_encode(args)


def preflight_checks(args):

    check_files_are_found([args.input_3di_alignment, args.input_aa_alignment])

    # check that alignments match in length

    if args.gap_threshold < 0 or args.gap_threshold > 1:
        report_message("`-g|--gap-threshold` must be between 0 and 1.")
        notify_premature_exit()

    if args.swap_3di_variability_threshold < 0 or args.swap_3di_variability_threshold > 1:
        report_message("`-s|--swap-3di-variability-threshold` must be between 0 and 1.")
        notify_premature_exit()

    if args.trim_variability_upper_percentile < 0 or args.trim_variability_upper_percentile > 100:
        report_message("`-t|--trim-variability-upper-percentile` must be between 0 and 100.")
        notify_premature_exit()
