from Bio import SeqIO # type: ignore
import random
import os
import sys
import subprocess
from tqdm import tqdm # type: ignore
from bit.modules.general import color_text

def generate_reads(args):

    preflight_checks(args)

    proportions = get_proportions(args)

    gen_reads(args, proportions)

    compress_with_pigz(args.output_prefix, read_type=args.type)


def preflight_checks(args):

    missing_files = [f for f in args.input_fastas if not os.path.exists(f)]

    if args.coverage:
        try:
            float(args.coverage)
        except ValueError:
            if not os.path.exists(args.coverage):
                missing_files.append(args.coverage)

    if args.proportions_file and not os.path.exists(args.proportions_file):
        missing_files.append(args.proportions_file)

    if missing_files:
        print(f"\n    Specified input files not found: {', '.join(missing_files)}")
        print("\n    Exiting for now :(\n")
        sys.exit(1)

    if args.proportions_file and args.coverage:
        print("\n    Cannot use both `--proportions-file` and `--coverage` options together.")
        print("\n    Exiting for now :(\n")
        sys.exit(1)


def get_proportions(args):

    if args.coverage:
        return compute_reads_from_coverage(args)

    proportions = parse_proportions_file(args.proportions_file, args.input_fastas)

    for fasta_file in args.input_fastas:
        if fasta_file not in proportions:
            print(f"\n   {fasta_file} is not specified in the proportions file.")
            print("\n    Exiting for now :(\n")
            sys.exit(1)

    return proportions


def compute_reads_from_coverage(args):

    coverages = parse_coverages_input(args.coverage, args.input_fastas)

    for fasta_file in args.input_fastas:
        if fasta_file not in coverages:
            print(f"\n    {fasta_file} is not specified in the coverage file.")
            print("\n    Exiting for now :(\n")
            sys.exit(1)

    # compute reads needed per genome: coverage * genome_size / read_length
    reads_per_file = {}

    for fasta_file in args.input_fastas:
        genome_size = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))

        if args.type == "paired-end":
            bases_per_fragment = min(2 * args.read_length, args.fragment_size)
            fragments_needed = round(coverages[fasta_file] * genome_size / bases_per_fragment)
            reads_per_file[fasta_file] = 2 * fragments_needed
        else:
            reads_per_file[fasta_file] = round(coverages[fasta_file] * genome_size / args.read_length)

    total_reads = sum(reads_per_file.values())

    if total_reads == 0:
        print("\n    Processing the coverage file resulted in 0 total reads. Check coverage values and genome sizes.")
        print("\n    Exiting for now :(\n")
        sys.exit(1)

    # override num_reads with computed total
    args.num_reads = total_reads

    print(f"\n    {color_text('Running in coverage-specified mode:', 'yellow')} making {total_reads:,} total reads from {len(args.input_fastas)} input file(s)")

    # return proportions derived from the per-file read counts
    proportions = {fasta_file: reads / total_reads for fasta_file, reads in reads_per_file.items()}

    return proportions


def parse_coverages_input(coverage_input, input_fastas):

    try:
        coverage_value = float(coverage_input)
        return {fasta_file: coverage_value for fasta_file in input_fastas}

    except ValueError:
        pass

    coverages = {}

    with open(coverage_input, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fasta_file, coverage = line.split('\t')
            coverages[fasta_file] = float(coverage)

    return coverages


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

        # normalizing proportions to ensure they sum to 1
        for key in proportions:
            proportions[key] = proportions[key] / total_proportion
        return proportions

    else:

        # assigning equal proportions if no proportions file is provided
        num_fastas = len(input_fastas)
        return {fasta_file: 1 / num_fastas for fasta_file in input_fastas}


def extract_subsequence(sequence, seq_length, subseq_len, circularize, include_Ns=False, max_attempts=100):

    for _ in range(max_attempts):

        if circularize:
            start = random.randint(0, seq_length - 1)
            end = start + subseq_len
            if end <= seq_length:
                subseq = sequence[start:end]
            else:
                subseq = sequence[start:] + sequence[: end - seq_length]
        elif subseq_len >= seq_length:
            # when the requested length meets or exceeds the contig length,
            # pick any start position and truncate at the contig end
            start = random.randint(0, seq_length - 1)
            subseq = sequence[start:]
        else:
            start = random.randint(0, seq_length - subseq_len)
            subseq = sequence[start:start + subseq_len]

        if include_Ns or 'N' not in subseq.upper():
            return subseq, start

    # fallback: return last attempt if all had Ns
    return subseq, start


def gen_reads(args, proportions):

    if args.seed is not None:
        random.seed(args.seed)

    if args.type in ("single-end", "long"):
        gen_single_reads(args, proportions)
    else:
        gen_paired_reads(args, proportions)


def gen_paired_reads(args, proportions):

    if args.num_reads % 2 != 0:
        args.num_reads += 1
        print(f"\n    Note: Paired-end mode requires an even number of reads. Rounding up to {args.num_reads}.")

    forward_reads_file = f"{args.output_prefix}_R1.fastq"
    reverse_reads_file = f"{args.output_prefix}_R2.fastq"

    pct = args.fragment_size_range / 100
    min_frag = max(1, int(args.fragment_size * (1 - pct)))
    max_frag = int(args.fragment_size * (1 + pct))

    with open(forward_reads_file, 'w') as fw, open(reverse_reads_file, 'w') as rw:

        print("")

        num_fragments = args.num_reads // 2
        reads_remaining = num_fragments
        remainder = 0.0
        read_count = 0

        pbar = tqdm(total=len(args.input_fastas), desc = "    Generating reads from fasta file(s)", ncols=90)
        for fasta_file in args.input_fastas:

            total_length = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))

            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id
                sequence = str(record.seq).upper()
                seq_length = len(sequence)

                exact = (seq_length / total_length) * num_fragments * proportions[fasta_file] + remainder
                entry_reads = min(round(exact), reads_remaining)
                remainder = exact - round(exact)

                # generate paired-end reads
                for _ in range(entry_reads):
                    frag_len = min(random.randint(min_frag, max_frag), seq_length)

                    fragment, start = extract_subsequence(sequence, seq_length, frag_len, args.circularize, args.include_Ns)

                    forward_read = fragment[:args.read_length]
                    reverse_read = fragment[-args.read_length:][::-1].translate(str.maketrans("ACGT", "TGCA"))

                    fwd_quality_scores = "I" * len(forward_read)
                    rev_quality_scores = "I" * len(reverse_read)

                    read_count += 1
                    fw.write(f"@{seq_id}_{read_count}_{start}/1\n")
                    fw.write(f"{forward_read}\n")
                    fw.write(f"+\n")
                    fw.write(f"{fwd_quality_scores}\n")

                    rw.write(f"@{seq_id}_{read_count}_{start}/2\n")
                    rw.write(f"{reverse_read}\n")
                    rw.write(f"+\n")
                    rw.write(f"{rev_quality_scores}\n")

                reads_remaining = reads_remaining - entry_reads
                if reads_remaining <= 0:
                    break

            pbar.update(1)
            if reads_remaining <= 0:
                break
        pbar.update(pbar.total - pbar.n)
        pbar.close()

        # if any reads remain due to rounding, adding them from the last contig
        if reads_remaining > 0 and sequence:

            for _ in range(reads_remaining):
                frag_len = min(random.randint(min_frag, max_frag), seq_length)
                fragment, start = extract_subsequence(sequence, seq_length, frag_len, args.circularize, args.include_Ns)

                forward_read = fragment[:args.read_length]
                reverse_read = fragment[-args.read_length:][::-1].translate(str.maketrans("ACGT", "TGCA"))
                quality_scores = "I" * args.read_length

                read_count += 1
                fw.write(f"@{seq_id}_{read_count}_{start}/1\n")
                fw.write(f"{forward_read}\n")
                fw.write(f"+\n")
                fw.write(f"{quality_scores}\n")

                rw.write(f"@{seq_id}_{read_count}_{start}/2\n")
                rw.write(f"{reverse_read}\n")
                rw.write(f"+\n")
                rw.write(f"{quality_scores}\n")


def gen_single_reads(args, proportions):

    reads_file = f"{args.output_prefix}.fastq"

    with open(reads_file, 'w') as fw:

        print("")

        if args.type == "long":
            pct = args.long_read_length_range / 100
            min_len = max(1, int(args.read_length * (1 - pct)))
            max_len = int(args.read_length * (1 + pct))

        reads_remaining = args.num_reads
        remainder = 0.0
        read_count = 0

        pbar = tqdm(total=len(args.input_fastas), desc = "    Generating reads from fasta file(s)", ncols=90)
        for fasta_file in args.input_fastas:

            total_length = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))

            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id
                sequence = str(record.seq)
                seq_length = len(sequence)

                exact = (seq_length / total_length) * args.num_reads * proportions[fasta_file] + remainder
                entry_reads = min(round(exact), reads_remaining)
                remainder = exact - round(exact)

                for _ in range(entry_reads):

                    if args.type == "long":
                        read_len = min(random.randint(min_len, max_len), seq_length)
                    else:
                        read_len = min(args.read_length, seq_length)

                    read, start = extract_subsequence(sequence, seq_length, read_len, args.circularize, args.include_Ns)

                    quality_scores = "I" * len(read)

                    read_count += 1
                    fw.write(f"@{seq_id}_{read_count}_{start}\n")
                    fw.write(f"{read}\n")
                    fw.write(f"+\n")
                    fw.write(f"{quality_scores}\n")

                reads_remaining = reads_remaining - entry_reads
                if reads_remaining <= 0:
                    break

            pbar.update(1)
            if reads_remaining <= 0:
                break
        pbar.update(pbar.total - pbar.n)
        pbar.close()

        # if any reads remain due to rounding, adding them from the last contig
        if reads_remaining > 0 and sequence:

            for _ in range(reads_remaining):
                if args.type == "long":
                    read_len = min(random.randint(min_len, max_len), seq_length)
                else:
                    read_len = min(args.read_length, seq_length)

                read, start = extract_subsequence(sequence, seq_length, read_len, args.circularize, args.include_Ns)

                quality_scores = "I" * len(read)
                read_count += 1
                fw.write(f"@{seq_id}_{read_count}_{start}\n")
                fw.write(f"{read}\n")
                fw.write(f"+\n")
                fw.write(f"{quality_scores}\n")


def compress_with_pigz(output_prefix, read_type="paired-end"):

    if read_type in ("single-end", "long"):
        reads_file = f"{output_prefix}.fastq"
        print("\n    Compressing output FASTQ file...")
        try:
            subprocess.run(["pigz", "-f", reads_file], check = True)
            print(f"\n    Compressed file written: {reads_file}.gz\n")
        except FileNotFoundError:
            print("pigz not found. You're on your own for compression!")
        except subprocess.CalledProcessError as e:
            print(f"Error during compression: {e}")
    else:
        forward_reads_file = f"{output_prefix}_R1.fastq"
        reverse_reads_file = f"{output_prefix}_R2.fastq"

        print("\n    Compressing output fastq files...")

        try:
            subprocess.run(["pigz", "-f", forward_reads_file], check = True)
            subprocess.run(["pigz", "-f", reverse_reads_file], check = True)
            print(f"\n    Compressed files written: {forward_reads_file}.gz, {reverse_reads_file}.gz\n")
        except FileNotFoundError:
            print("pigz not found. You're on your own for compression!")
        except subprocess.CalledProcessError as e:
            print(f"Error during compression: {e}")
