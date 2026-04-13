import sys
import argparse
import datetime
import random
from bit.cli.common import CustomRichHelpFormatter, add_help, reconstruct_invocation
from bit.modules.general import check_files_are_found, report_message


def build_parser():

    desc = """
        This script will mutate all sequences of a nucleotide or amino-acid multifasta with
        the specified mutation rate. By default it only swaps bases, but it can optionally introduce indels also.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-mutate-seqs -i input.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-i",
        "--input-fasta",
        metavar = "<FILE>",
        help = "Starting fasta file",
        action = "store",
        required = True
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
        "-m",
        "--mutation-rate",
        metavar = "<FLOAT>",
        help = 'Wanted mutation rate, must be a float between 0 and 1 (default: "0.01")',
        action = "store",
        type = float,
        default = 0.01
    )

    optional.add_argument(
        "-T",
        "--ti-tv-ratio",
        metavar = "<FLOAT>",
        help = "Specify a transition/transversion ratio for nucleotide sequences (default = 1.0).",
        action = "store",
        type = float,
        default = 1.0
    )

    optional.add_argument(
        "-I",
        "--indel-rate",
        metavar = "<FLOAT>",
        help = """
            Specify the frequency of indels wanted, must be a float between 0 and 1 (default = 0.0).
            This is a subset of the overall mutation rate provided via the '-m' argument.
            E.g., providing '-I 0.2' means 1 out of every 5 generated mutations will be an indel.
            By default, no indels are incorporated. The likelihood of adding an insertion/deletion
            is 50/50, with all characters equally likely to be inserted.
            """,
        action = "store",
        type = float,
        default = 0.0
    )

    optional.add_argument(
        "-o",
        "--output-prefix",
        metavar = "<STR>",
        help = 'Output prefix (default: "mutated").',
        action = "store",
        default = "mutated"
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

    full_command = reconstruct_invocation(parser, args)

    available_substitutions, seed = preflight_checks_and_setup(args)

    from Bio import SeqIO # type: ignore
    from bit.modules.seqs import mutate_seq

    # doing mutations per sequence
    with (open (args.input_fasta, "r") as in_fasta,
          open(f"{args.output_prefix}.fasta", "w") as out_fasta,
          open(f"{args.output_prefix}.tsv", "w") as tsv,
          open(f"{args.output_prefix}.log", "w") as log):

        # writing primary log file
        log.write(f"Full command executed: {full_command}\n")
        log.write(f"Mutation rate utilized: {args.mutation_rate}\n")
        if args.molecule_type == "NT":
            log.write(f"Transition/transversion ratio utilized: {args.ti_tv_ratio}\n")
        log.write(f"Indel rate utilized: {args.indel_rate}\n")
        log.write(f"Seed used: {seed}\n")

        tsv.write("seq_id\tseq_length\tnum_total_changes_made\tnum_substitutions\tnum_transitions\tnum_transversions\tnum_total_indels\tnum_insertions\tnum_deletions\n")

        for seq_record in SeqIO.parse(in_fasta, "fasta"):

            (seq, total_num_mutations,
             num_substitutions, num_indels,
             num_insertions, num_deletions,
             num_transitions, num_transversions) = mutate_seq(seq_record.seq,
                                                              args.molecule_type,
                                                              available_substitutions,
                                                              args.mutation_rate,
                                                              args.ti_tv_ratio,
                                                              args.indel_rate)

            out_fasta.write(f">{seq_record.id}\n")
            out_fasta.write(f"{seq}\n")

            # writing out log file
            tsv.write(f"{seq_record.id}\t{len(seq_record)}\t{total_num_mutations}\t{num_substitutions}\t{num_transitions}\t{num_transversions}\t{num_indels}\t{num_insertions}\t{num_deletions}\n")

    report_outputs(out_fasta, tsv, log)


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


def report_outputs(out_fasta, tsv, log):

    report_message("Seqs mutated!", color = "green", initial_indent="    ")
    report_message(f"Output fasta:        {out_fasta.name}", color = "none", initial_indent="        ", subsequent_indent="        ", join=False)
    report_message(f"Output summary tsv:  {tsv.name}", color = "none", initial_indent="        ", subsequent_indent="        ", join=False, leading_newline=False)
    report_message(f"Output log:          {log.name}", color = "none", initial_indent="        ", subsequent_indent="        ", join=False, leading_newline=False)
    print()
