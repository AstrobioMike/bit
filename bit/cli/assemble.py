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

    inputs = parser.add_argument_group("INPUT READS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    inputs.add_argument(
        "-1",
        "--read-1",
        metavar="<FILE>",
        help="Input read 1 file (gzipped fastq format; incompatible with --reads-dir)",
        type=Path,
    )
    inputs.add_argument(
        "-2",
        "--read-2",
        metavar="<FILE>",
        help="Input read 2 file (gzipped fastq format; incompatible with --reads-dir)",
        type=Path,
    )
    inputs.add_argument(
        "-r",
        "--reads-dir",
        metavar="<DIR>",
        help="Directory containing gzipped read files (default: current directory; incompatible with -1 and -2)",
        type=Path,
        default=Path("."),
    )
    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help="Directory for the output files (default: 'assembly/')",
        type=Path,
        default=Path("assembly"),
    )
    optional.add_argument(
        "-a",
        "--assembler",
        choices=["megahit", "spades"],
        help="Assembler to use (default: megahit)",
        default="megahit",
    )
    optional.add_argument(
        "--isolate",
        action="store_true",
        help="Run assembly in '--isolate' mode rather than '--meta' (only applies to spades assembler)",
    )
    optional.add_argument(
        "--skip-fastp",
        action="store_true",
        help="Skip fastp quality trimming/filtering",
    )
    optional.add_argument(
        "--skip-bbnorm",
        action="store_true",
        help="Skip bbnorm digital normalization",
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
    full_cmd_executed = reconstruct_invocation(parser, args)

    print(args)
    print(full_cmd_executed)
