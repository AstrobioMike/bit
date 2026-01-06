from Bio import SeqIO
import random
import os
import sys
import subprocess
from tqdm import tqdm

def generate_reads(args):

    preflight_checks(args)

    simulate_paired_end_reads(args)

    compress_with_pigz(args.output_prefix)


def preflight_checks(args):

    missing_files = [f for f in args.input_fastas if not os.path.exists(f)]

    if args.proportions_file and not os.path.exists(args.proportions_file):
        missing_files.append(args.proportions_file)

    if missing_files:
        print(f"\n    Specified input files not found: {', '.join(missing_files)}")
        print("\n    Exiting for now :(\n")
        sys.exit(1)


def parse_proportions_file(proportions_file, input_fastas):

    if proportions_file:
        proportions = {}
        total_proportion = 0
        with open(proportions_file, 'r') as f:
            for line in f:
                fasta_file, proportion = line.strip().split('\t')
                proportion = float(proportion)
                proportions[fasta_file] = proportion
                total_proportion += proportion

        # Normalize proportions to ensure they sum to 1
        for key in proportions:
            proportions[key] /= total_proportion

        return proportions

    else:
        # Assign equal proportions if no proportions file is provided
        num_fastas = len(input_fastas)
        return {fasta_file: 1 / num_fastas for fasta_file in input_fastas}


def simulate_paired_end_reads(args):

    proportions = parse_proportions_file(args.proportions_file, args.input_fastas)

    if args.seed is not None:
        random.seed(args.seed)

    # ensure all input FASTA files are accounted for in the proportions file
    for fasta_file in args.input_fastas:
        if fasta_file not in proportions:
            raise ValueError(f"{fasta_file} is not specified in the proportions file.")

    forward_reads_file = f"{args.output_prefix}_R1.fastq"
    reverse_reads_file = f"{args.output_prefix}_R2.fastq"

    with open(forward_reads_file, 'w') as fw, open(reverse_reads_file, 'w') as rw:

        print("")

        for fasta_file in tqdm(args.input_fastas, desc = "    Generating reads from each input FASTA file"):

            total_length = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))

            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id
                sequence = str(record.seq)
                seq_length = len(sequence)

                # proportional reads for this sequence entry
                entry_reads = int((seq_length / total_length) * args.num_read_pairs * proportions[fasta_file])

                # generate paired-end reads
                for _ in range(entry_reads):
                    start = random.randint(0, max(0, seq_length - args.fragment_size))
                    fragment = sequence[start:start + args.fragment_size]

                    if len(fragment) < args.fragment_size:
                        continue

                    forward_read = fragment[:args.read_length]
                    reverse_read = fragment[-args.read_length:][::-1].translate(str.maketrans("ACGT", "TGCA"))

                    # generating FASTQ quality scores
                    quality_scores = "I" * args.read_length

                    # writing forward read
                    fw.write(f"@{seq_id}_{start}/1\n")
                    fw.write(f"{forward_read}\n")
                    fw.write(f"+\n")
                    fw.write(f"{quality_scores}\n")

                    # writing reverse read
                    rw.write(f"@{seq_id}_{start}/2\n")
                    rw.write(f"{reverse_read}\n")
                    rw.write(f"+\n")
                    rw.write(f"{quality_scores}\n")


def compress_with_pigz(output_prefix):

    forward_reads_file = f"{output_prefix}_R1.fastq"
    reverse_reads_file = f"{output_prefix}_R2.fastq"

    print("\n    Compressing output FASTQ files...")

    try:
        subprocess.run(["pigz", "-f", forward_reads_file], check = True)
        subprocess.run(["pigz", "-f", reverse_reads_file], check = True)
        print(f"\n    Compressed files written: {forward_reads_file}.gz, {reverse_reads_file}.gz\n")
    except FileNotFoundError:
        print("pigz not found. Please install pigz or compress files manually.")
    except subprocess.CalledProcessError as e:
        print(f"Error during compression: {e}")
