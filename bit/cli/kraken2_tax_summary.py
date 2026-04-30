#!/usr/bin/env python
import os
import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)
from bit.modules.general import check_files_are_found, report_message, notify_premature_exit

def build_parser():

    desc = """
        This script can take an individual or multiple kraken2 (or bracken) reports and generates a combined table with full-lineage information,
        counts, and percentages. The lineage info is constructed directly from the report files. Completely unclassified things will
        be reported as "Unclassified" at each rank. If you see NA for all ranks, it's likely from something like "cellular organisms"
        or "root". Check for the taxid in the initial kraken2/bracken report to find out. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-kraken2-tax-summary -i kraken2.report",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-reports",
        metavar="<FILE>",
        nargs="+",
        help="Input kraken2 or braken report file(s) (space-delimited if multiple)",
        required=True,
    )

    optional.add_argument(
        "-s",
        "--sample-names",
        metavar="<NAME>",
        nargs="+",
        help='Wanted sample name(s) (space-delimited if multiple). Be sure it matches the order of the input files (by default will use basename of input files up to last period)',
        default=None
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help="Output TSV file (default: tax-summary.tsv)",
        default="tax-summary.tsv",
    )

    add_help(optional)

    return parser

def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    input_map = preflight_checks(args.input_reports, args.sample_names)

    from bit.modules.kraken2_tax_summary import kraken2_tax_summary

    kraken2_tax_summary(input_map, args.output_tsv)


def preflight_checks(input_reports, sample_names):

    check_files_are_found(input_reports)

    input_map = check_and_setup_file_to_sample_map(input_reports, sample_names)

    return input_map

def check_and_setup_file_to_sample_map(input_reports, sample_names):

    input_map = {}

    if sample_names is None or len(sample_names) == 0:

        for file in input_reports:
            curr_sample = os.path.basename(file).rsplit('.', 1)[0]
            if curr_sample in input_map.values():
                report_message('\n    It seems the file "' + file + '" would have an overlapping basename with another input report. That\'s not gonna fly.\n')
                report_message("\n    You can provide custom sample names with the `-s/--sample-names` flag to get around this.\n", color="none")
                notify_premature_exit()
            input_map[file] = curr_sample


    if sample_names is not None and len(sample_names) != len(input_reports):
        report_message("\n    It seems the number of provided sample names doesn't match the number of provided input files :(")
        report_message("\n    Check usage with `bit-kraken2-tax-summary -h`.\n")
        notify_premature_exit()

    if sample_names is not None:
        for file, name in zip(input_reports, sample_names):
            if name in input_map.values():
                report_message('\n    It seems the name "' + name + '" was given more than once in the provided sample names. That\'s not gonna fly.\n')
                report_message("\n    Please make sure all sample names provided with the `-s/--sample-names` flag are unique.\n", color="none")
                notify_premature_exit()
            input_map[file] = name

    return input_map
