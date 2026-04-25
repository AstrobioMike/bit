import sys
import argparse
from pathlib import Path
from rich_argparse import RawTextRichHelpFormatter # type: ignore
from bit.cli.common import (wrap_help,
                            wrap_multiline_help,
                            add_help,
                            add_common_snakemake_arguments,
                            reconstruct_invocation)
from bit.modules.general import report_message, notify_premature_exit

RawTextRichHelpFormatter.group_name_formatter = lambda name: "Usage" if name.lower() == "usage" else name


def build_parser():

    raw_desc = (
        "This program runs an assembly workflow with optional QC and digital normalization "
        "(short-read and paired-end only currently). For version info, run `bit-version`."
    )

    desc = wrap_help(raw_desc, 4)

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-assemble -1 R1.fastq.gz -2 R2.fastq.gz` or `bit-assemble -r reads-dir/`",
        formatter_class=RawTextRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters (choose one input method)")
    general = parser.add_argument_group("General Parameters")
    assembly = parser.add_argument_group('Assembly Parameters')
    snakemake = parser.add_argument_group('Snakemake Parameters')

    required.add_argument(
        "-1",
        "--read-1",
        metavar="<FILE>",
        help=wrap_help("Input read 1 file (gzipped fastq format; good for when you have one sample; incompatible with --reads-dir)"),
        type=Path,
    )
    required.add_argument(
        "-2",
        "--read-2",
        metavar="<FILE>",
        help=wrap_help("Input read 2 file (gzipped fastq format; good for when you have one sample; incompatible with --reads-dir)"),
        type=Path,
    )
    required.add_argument(
        "-r",
        "--reads-dir",
        metavar="<DIR>",
        help=wrap_help("Directory containing gzipped read files (helpful when you have multiple samples; incompatible with -1 and -2)"),
        type=Path,
    )
    general.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help=wrap_help("Directory for the output files (default: 'assembly/')"),
        type=Path,
        default=Path("assembly"),
    )
    general.add_argument(
        "-t",
        "--threads",
        metavar="<INT>",
        help=wrap_help("Number of threads/cpus to pass to QC and assembly commands (default: 1; may be multiplied by number of snakemake jobs)"),
        default=1,
        type=int,
    )
    general.add_argument(
        "--run-fastp",
        action="store_true",
        help=wrap_help("Run fastp quality trimming/filtering"),
    )
    general.add_argument(
        "--run-bbnorm",
        action="store_true",
        help=wrap_help("Run bbnorm digital normalization"),
    )

    add_help(general)

    assembly.add_argument(
        "-a",
        "--assembler",
        choices=["megahit", "spades"],
        help=wrap_help("Assembler to use (default: megahit)"),
        default="megahit",
    )
    assembly.add_argument(
        "--min-contig-len",
        metavar="<INT>",
        help=wrap_help("Minimum contig length to report (default: 500)"),
        default=500,
        type=int,
    )
    assembly.add_argument(
        "--isolate",
        action="store_true",
        help=wrap_help("Run assembly in '--isolate' mode rather than '--meta' (spades assembler only)"),
    )
    assembly.add_argument(
        "-m",
        "--memory",
        metavar="<MEM>",
        help=wrap_multiline_help(
            "Memory setting (may be multiplied by number of snakemake jobs)\n"
            "    For megahit: max memory to use (default: '20e9'); may be a fraction (e.g., '0.5') or bytes (e.g., '20e9' would be 20 GB)\n"
            "    For spades:  max memory in Gb (terminates if exceeded; default: '250')"
        ),
    )

    add_common_snakemake_arguments(snakemake)

    return parser


def main():

    parser = build_parser()
    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    check_required_inputs(args)

    args = check_other_settings(args)

    full_cmd_executed = reconstruct_invocation(parser, args)

    from bit.modules.assemble import run_assembly

    run_assembly(args, full_cmd_executed)


def check_required_inputs(args):
    if args.reads_dir and (args.read_1 or args.read_2):
        report_message("We cannot accept a --reads-dir AND -1 and/or -2 inputs. "
                       "Please use one method or the other to specify inputs. -1 "
                       "and -2 are good when you have a single sample to do. The "
                       "--reads-dir is good for when you have many in a directory.",
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


def check_other_settings(args):

        if args.assembler != "spades" and args.isolate:
            report_message("The --isolate option is only compatible with the spades assembler. Please remove the --isolate option or switch to the spades assembler.",
                            initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()

        args = normalize_memory(args)

        return args

def normalize_memory(args):
    if args.assembler == "megahit":
        args.memory = args.memory or "20e9"
    elif args.assembler == "spades":
        val = args.memory or 250
        try:
            args.memory = int(val)
        except ValueError:
            report_message(f"Memory value '{val}' is not a valid integer for SPAdes. Please provide an integer value in Gb (e.g., '250').",
                           initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()
    return args
