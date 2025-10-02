#!/usr/bin/env python
import sys
import argparse
from bit.modules.seqs import check_for_fastq_dup_headers
from bit.modules.general import report_message
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)


def main():

    desc = """
        This script reports any duplicate headers in a fastq file. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-check-for-fastq-dup-headers -i reads.fastq.gz`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group('REQUIRED PARAMETERS')
    optional = parser.add_argument_group('OPTIONAL PARAMETERS')

    required.add_argument(
        "-i",
        "--input-fastq",
        metavar="<FILE>",
        help="Input fastq file (can be gzipped or not)",
        required=True,
    )

    optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help="Output file to record duplicate headers (will just print to screen if not specified)",
    )

    add_help(optional)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)


    args = parser.parse_args()

    dup_keys, seq_count = check_for_fastq_dup_headers(args.input_fastq, args.output_file)

    if len(dup_keys) > 0:

        if args.output_file:
            message = f"{len(dup_keys)} duplicate header(s) found among {seq_count} input fastq entries. They were written to '{args.output_file}'.\n"
            report_message(message, color="yellow")
            print()
            with open(args.output_file, "w") as out:
                out.write("\n".join(dup_keys))
                out.write("\n")
        else:
            message = f"{len(dup_keys)} duplicate header(s) found among {seq_count} input fastq entries:\n\n"
            report_message(message)
            message = "\n     ".join(dup_keys) + "\n"
            print(f"\n     {message}")
            report_message("If you want to save these to a file, re-run this command with the -o option.")
            print()
    else:
        message = f"There were no duplicate headers detected among the {seq_count} input fastq entries :)"
        report_message(message, color="green")
        print("")
