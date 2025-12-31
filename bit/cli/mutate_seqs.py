import sys
import argparse
import datetime
import random
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import check_files_are_found
from Bio import SeqIO
from bit.modules.seqs import mutate_seq


def build_parser():

    desc = """
        This script will mutate all sequences of a nucleotide or amino-acid multifasta with
        the specified mutation rate. It does not take into consideration transition/transversion
        rates. By default it only swaps bases, but it can optionally introduce indels also.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-mutate-seqs -i input.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar = "<FILE>",
        help = "Starting fasta file",
        action = "store",
        required = True
    )

    optional.add_argument(
        "-m",
        "--mutation-rate",
        metavar = "<FLOAT>",
        help = 'Wanted mutation rate, must be a float between 0 and 1 (default: "0.01")',
        action = "store",
        type = float,
        default = 0.01
    )

    optional.add_argument(
        "-t",
        "--molecule-type",
        help = "'NT' for nucleotides or 'AA' for amino acids (default = 'NT')",
        action = "store",
        choices = ["NT", "AA"],
        default = 'NT'
    )

    optional.add_argument(
        "-o",
        "--output-fasta",
        metavar = "<FILE>",
        help = 'Output mutated fasta file (default: "mutated.fasta").',
        action = "store",
        default = "mutated.fasta"
    )

    optional.add_argument(
        "-l",
        "--output-log",
        metavar = "<FILE>",
        help = 'Output file reporting the seed used and how many characters were changed in each input sequence (default: "mutated.log")',
        action = "store",
        default = "mutated.log"
    )

    optional.add_argument(
        "-I",
        "--indel-rate",
        metavar = "<FLOAT>",
        help = """
            Specify the frequency of indels wanted, must be a float between 0 and 1. This
            is a subset of the overall mutation rate provided via the '-m' argument.
            E.g., providing '-f 0.2' means 1 out of every 5 generated mutations will be an indel.
            By default, no indels are incorporated (default = 0.0). The likelihood
            of adding an insertion/deletion is 50/50, with all characters equally likely
            to be inserted.
            """,
        action = "store",
        type = float,
        default = 0.0
    )

    optional.add_argument(
        "-s",
        "--seed",
        metavar = "<INT>",
        help = "Optionally set the randomization seed",
        action = "store",
        type = int
    )

    add_help(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    available_substitutions, seed = preflight_checks_and_setup(args)

    # doing mutations per sequence
    with (open (args.input_fasta, "r") as in_fasta,
          open(args.output_fasta, "w") as out_fasta,
          open(args.output_log, "w") as log):

        # starting log file
        log.write(f"# seed used: {seed}\n")
        log.write(f"# indel frequency utilized: {args.indel_rate}\n")
        log.write("# the rest below is a tab-delimited table summarizing the changes introduced to each input sequence\n")
        log.write("seq_id\tseq_length\tnum_total_changes_made\tnum_substitutions\tnum_total_indels\tnum_insertions\tnum_deletions\n")

        for seq_record in SeqIO.parse(in_fasta, "fasta"):

            (seq, total_num_mutations, num_substitutions, num_indels,
            num_insertions, num_deletions) = mutate_seq(seq_record.seq, available_substitutions,
                                                                              args.mutation_rate, args.indel_rate)

            out_fasta.write(f">{seq_record.id}\n")
            out_fasta.write(f"{seq}\n")

            # writing out log file
            log.write(f"{seq_record.id}\t{len(seq_record)}\t{total_num_mutations}\t{num_substitutions}\t{num_indels}\t{num_insertions}\t{num_deletions}\n")


def preflight_checks_and_setup(args):

    check_files_are_found([args.input_fasta])

    # checking mutation rate is within 0 and 1
    validate_mutation_rate(args.mutation_rate)

    # checking indel frequency rate is within 0 and 1
    validate_indel_rate(args.indel_rate)

    # setting substitution type
    available_substitutions = get_available_substitutions(args.molecule_type)

    # setting seed if provided, and recording it for the output log either way
    seed = set_seed(args.seed)

    return available_substitutions, seed


def validate_mutation_rate(mutation_rate):

    if not 0 <= mutation_rate <= 1:
        print("\n    The '--mutation-rate' argument needs to be between 0 and and 1.\n")
        print("  Exiting for now.\n")
        sys.exit(1)


def validate_indel_rate(indel_rate):

    if not 0 <= indel_rate <= 1:
        print("\n    The '--indel-rate' argument needs to be between 0 and and 1.\n")
        print("  Exiting for now.\n")
        sys.exit(1)


def get_available_substitutions(molecule_type):

    if molecule_type == 'NT':
        return ['A', 'T', 'C', 'G']

    else:
        return ['A', 'R', 'N', 'D',
                'C', 'E', 'Q', 'G',
                'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S',
                'T', 'W', 'Y', 'V']


def set_seed(input_seed):

    if input_seed:

        seed = input_seed

    else:
        # setting seed manually (pseudo-randomly) BECAUSE it seems if you don't set it, you can't
        # pull it out otherwise (and i want to log it whether it was specified or not by the user)
        time = datetime.datetime.now()
        seed = time.hour * 10000 + time.minute * 100 + time.second

    random.seed(seed)

    return seed
