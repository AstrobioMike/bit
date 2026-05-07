import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_seed
from bit.modules.general import check_files_are_found, report_message, notify_premature_exit

def build_parser():

    desc = """
        This program trims and joint encodes a 3di alignment and an amino-acid alignment (that is derived from the 3di alignment).
        It removes columns based on a gap threshold, and it swaps amino-acid characters into the base 3di alignment
        based on specified column-variability criteria (see `bit-calc-variability-is-msa -h` for details on "variability" here).
        It produces the joint-encoded alignment and a partitions file required for appropriate phylogenetics analysis.
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
        "--3di-lower-bound",
        metavar="<FLOAT>",
        type=float,
        help="""
            For a given column, if the variability in 3di space is <= this value,
            it will be swapped to the corresponding AA character IF the variability
            in AA space is > the 3di variability AND < the AA-upper-bound.
            If those conditions are not met, the 3di character will be retained. (default: 0.2)
            """,
        default=0.2,
        dest="threedi_lower_bound"
    )
    optional.add_argument(
        "--3di-upper-bound",
        metavar="<FLOAT>",
        type=float,
        help="""
            For a given column, if the variability in 3di space is >= this value,
            it will be swapped to the corresponding AA character IF the variability
            in AA space is < the AA-upper-bound. If those conditions are not met,
            the column will be trimmed. (default: 0.7)
            """,
        default=0.7,
        dest="threedi_upper_bound"
    )
    optional.add_argument(
        "--AA-upper-bound",
        metavar="<FLOAT>",
        type=float,
        help="""
            For a given column, when an amino-acid character is being considered to replace a 3di
            character, it must have a variability < this value in order to be swapped in. If those
            conditions are not met, the column will be trimmed. (default: 0.8)
            """,
        default=0.8,
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

    # need check that alignments match in length

    bounds_params = [
        ("gap_threshold", "-g|--gap-threshold"),
        ("threedi_lower_bound", "--3di-lower-bound"),
        ("threedi_upper_bound", "--3di-upper-bound"),
        ("AA_upper_bound", "--AA-upper-bound"),
    ]

    for attr, flag in bounds_params:
        val = getattr(args, attr)
        if val < 0 or val > 1:
            report_message(f"`{flag}` must be between 0 and 1.")
            notify_premature_exit()
