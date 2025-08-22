#!/usr/bin/env python
import sys
import argparse
from bit.modules.kraken2_to_taxon_summaries import kraken2_to_taxon_summaries
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)

def main():

    desc = """
        This script generates a full-lineage count table from a kraken2 report. Completely unclassified things will be reported as
        "Unclassified" at each rank. If you see NA for all ranks, it's likely from something like "cellular organisms"
        or "root". Check for the taxid in the initial kraken report to find out. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-kraken2-to-taxon-summaries -i kraken2.report -o output.tsv",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-report",
        metavar="<FILE>",
        help="Input kraken2 report file",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help="Output TSV file (default: output.tsv)",
        default="output.tsv",
    )

    add_help(optional)

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    kraken2_to_taxon_summaries(args.input_report, args.output_tsv)
