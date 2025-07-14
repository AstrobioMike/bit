import sys
import argparse
from pathlib import Path
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            reconstruct_invocation)
# from bit.modules.assemble import run_assembly

def main():

    desc = """
        This program runs an assembly pipeline with optional QC and digital normalization
        (short-read and paired-end only currently). For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-assemble -1 R1.fastq.gz -2 R2.fastq.gz -c 8` or `bit-assemble -r reads-dir/ -c 8`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    optional.add_argument(
        "-1",
        "--read-1",
        metavar="<FILE>",
        help="Input read 1 file (gzipped fastq format; mutually exclusive with --reads-dir)",
        type=Path,
    )
    optional.add_argument(
        "-2",
        "--read-2",
        metavar="<FILE>",
        help="Input read 2 file (gzipped fastq format; mutally exclusive with --reads-dir)",
        type=Path,
    )
    optional.add_argument(
        "-r",
        "--reads-dir",
        metavar="<DIR>",
        help="Directory containing gzipped read files (default: current directory)",
        type=Path,
        default=Path("."),
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

    print(args)
