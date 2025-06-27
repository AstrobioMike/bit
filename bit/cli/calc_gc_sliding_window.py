import sys
import argparse
from bit.modules.seqs import calc_gc_sliding_window
from bit.cli.common import CustomRichHelpFormatter

def main():

    desc = """
        This script is for nucleotide multifastas and will return
        a tab-delimited file with 4 columns: header, sequence length,
        gc of whole sequence, and gc of each window of the specified
        window size (-w) for each step of the specified step size (-s).
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-calc-gc-per-seq -i input.fasta -o output.tsv",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar="<FILE>",
        help='Input fasta file',
        required=True,
    )
    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv file (default: "GC-windows-out.tsv")',
        default="GC-windows-out.tsv",
    )
    optional.add_argument(
        "-w",
        "--window-size",
        metavar="<INT>",
        help="Desired size of sliding window (default: 100)",
        default=100,
        type=int,
    )
    optional.add_argument(
        "-s",
        "--step-size",
        metavar="<INT>",
        help="Desired size of steps between each window (default: 10)",
        default=10,
        type=int,
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

    gc_window_dict = calc_gc_sliding_window(input_fasta=args.input_fasta,
                                            window = args.window_size,
                                            step = args.step_size)
    with open(args.output_tsv, "w") as out_tsv:
        out_tsv.write(f"header\tlength\tgc\tgc_per_window_size_{args.window_size}_with_step_size_{args.step_size}\n")
        for item in gc_window_dict:
            windows_string = ", ".join(f"{x:.2f}" for x in item['gc_of_windows'])
            out_tsv.write(f"{item['header']}\t{item['length']}\t{item['gc']}\t{windows_string}\n")
