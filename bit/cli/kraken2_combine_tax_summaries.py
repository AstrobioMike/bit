#!/usr/bin/env python
import os
import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)
from bit.modules.general import report_message, notify_premature_exit


def build_parser():

    desc = """
        This combines outputs from the `bit-kraken2-tax-summary` program. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-kraken2-combine-tax-summaries -i tax-summary-1.tsv tax-summary-2.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-files",
        metavar="<FILE>",
        nargs="+",
        help="Spade-delimited list of `bit-kraken2-tax-summary` output files, can be provided with shell wildcards",
        required=True
    )

    optional.add_argument(
        "-n",
        "--sample-names",
        metavar="<NAME>",
        nargs="+",
        help='Space-delimited list of wanted sample names. Be sure it matches the order of the input files (by default will use basename of input files up to last period)',
        default=''
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<STR>",
        help='Output combined summaries tsv (default: "combined-kraken2-tax-summaries.tsv")',
        default="combined-kraken2-tax-summaries.tsv"
    )

    add_help(optional)

    return parser


def main():
    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    all_samples = preflight_checks(args)

    from bit.modules.kraken2_combine_tax_summaries import combine_kraken2_tax_summaries

    combine_kraken2_tax_summaries(all_samples, args.output_file)


def preflight_checks(args):

    all_samples = {}

    if len(args.sample_names) == 0:
        for file in args.input_files:
            curr_sample = os.path.basename(file).rsplit('.', 1)[0]

            if file in all_samples:
                report_message('\n    It seems the file "' + file + '" is trying to get in here twice.')
                report_message("\n    That's not gonna fly :(\n")
                notify_premature_exit()
                sys.exit(1)

            all_samples[file] = curr_sample

    else:

        # checking if sample names provided the length equals the number of input files
        if len(args.sample_names) != len(args.input_files):
            report_message("\n    It seems the number of provided sample names doesn't match the number of provided input files :(")
            report_message("\n    Check usage with `bit-kraken2-combine-tax-summaries -h`.\n")
            notify_premature_exit()
            sys.exit(0)

        # setting iterator
        i = 0

        for curr_sample in args.sample_names:

            all_samples[args.input_files[i]] = curr_sample
            i += 1

    return all_samples
