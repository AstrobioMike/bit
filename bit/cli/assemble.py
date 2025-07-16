import sys
import argparse
from pathlib import Path
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_common_snakemake_arguments,
                            reconstruct_invocation)
from bit.modules.general import report_message, notify_premature_exit
# from bit.modules.assemble import run_assembly


def build_parser():

    desc = """
        This program runs an assembly workflow with optional QC and digital normalization
        (short-read and paired-end only currently). For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-assemble -1 R1.fastq.gz -2 R2.fastq.gz` or `bit-assemble -r reads-dir/`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")
    megahit = parser.add_argument_group('MEGAHIT PARAMETERS')
    spades = parser.add_argument_group('SPADES PARAMETERS')
    snakemake = parser.add_argument_group('SNAKEMAKE PARAMETERS')

    required.add_argument(
        "-1",
        "--read-1",
        metavar="<FILE>",
        help="Input read 1 file (gzipped fastq format; good for when you have one sample; incompatible with --reads-dir)",
        type=Path,
    )
    required.add_argument(
        "-2",
        "--read-2",
        metavar="<FILE>",
        help="Input read 2 file (gzipped fastq format; good for when you have one sample; incompatible with --reads-dir)",
        type=Path,
    )
    required.add_argument(
        "-r",
        "--reads-dir",
        metavar="<DIR>",
        help="Directory containing gzipped read files (helpful when you have multiple samples; incompatible with -1 and/or -2)",
        type=Path,
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
        "-t",
        "--threads",
        metavar="<INT>",
        help="Number of threads/cpus to pass to assembler (default: 2; may be multiplied by number of snakemake jobs)",
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
    add_help(optional)

    megahit.add_argument(
        "-m",
        "--memory",
        metavar="<STR>",
        help="Max memory to use (default: '20e9'). Can be specified as \
            fraction of total memory (e.g., '0.5') or value in bytes (e.g., '20e9' \
            would be 20 GB)",
        default="30e9",
        type=str,
    )
    megahit.add_argument(
        "--min-contig-len",
        metavar="<INT>",
        help="Minimum contig length to report (default: 500)",
        default=500,
        type=int,
    )

    spades.add_argument(
        "--isolate",
        action="store_true",
        help="Run assembly in '--isolate' mode rather than '--meta'",
    )

    add_common_snakemake_arguments(snakemake)

    return parser


def main():
    parser = build_parser()
    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    # doing checks on inputs
    check_required_inputs(args)

    # reconstructing the full command-line invocation for logging
    full_cmd_executed = reconstruct_invocation(parser, args)

    print(args)
    print(full_cmd_executed)


def check_required_inputs(args):
    if args.reads_dir and (args.read_1 or args.read_2):
        report_message("We cannot accept a --reads-dir AND -1 and/or -2 inputs. "
                       "Please use one method or the other to specify inputs. -1 "
                       "and -2 are good when you have a single sample to do. The "
                       "--reads-dir is for when you have many in a directory.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()

    if not args.reads_dir and (not args.read_1 or not args.read_2):
        report_message("You need to specify the input reads either by pointing to "
                       "them specifically with the -1 and -2 parameters or pointing "
                       "to the directory holding them with the --reads-dir parameter.",
                       initial_indent = "    ", subsequent_indent = "    ")
        notify_premature_exit()

    if args.reads_dir:
        reads_dir = Path(args.reads_dir)
        if not reads_dir.is_dir():
            print(f"Error: The specified reads directory '{reads_dir}' does not exist.")
            sys.exit(1)
