from Bio import SeqIO # type: ignore
import gzip
import random
import os
import sys
import subprocess
from tqdm import tqdm # type: ignore
from bit.modules.general import color_text


def parse_fasta(fasta_file):
    """ parse a fasta that may be plain text or gzip-compressed (.gz). """
    if str(fasta_file).endswith(".gz"):
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file, "r")
    try:
        for record in SeqIO.parse(handle, "fasta"):
            yield record
    finally:
        handle.close()

def _status(args, msg):
    """ print a status line unless the caller (e.g., gen-metagenome) set quiet """
    if not getattr(args, "quiet", False):
        print(msg)


def generate_reads(args):

    preflight_checks(args)

    proportions = get_proportions(args)

    gen_reads(args, proportions)

    compress_with_pigz(args.output_prefix, read_type=args.type, quiet=getattr(args, 'quiet', False))


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
        genome_size = sum(len(record.seq) for record in parse_fasta(fasta_file))

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

    _status(args, f"\n    {color_text('Running in coverage-specified mode:', 'yellow')} making {total_reads:,} total reads from {len(args.input_fastas)} input file(s)")

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


def compute_mate_coords(frag_start, fragment, fwd_len, rev_len, seq_length):
    """
    given the fragment's reference start, the realized fragment string, and the
    realized forward/reverse read lengths, return each mate's own 0-based,
    end-exclusive reference span plus a wrap flag.

    these spans are computed in the reference's forward orientation: r1 covers the
    low-coordinate end of the fragment footprint, r2 the high-coordinate end. when the
    sampled fragment comes from the minus strand, the caller is responsible for swapping
    the pairs so each reported span matches the mate that actually covers it (see gen_paired_reads)
    """
    frag_len_actual = len(fragment)

    r1_start = frag_start
    r1_end = frag_start + fwd_len

    r2_end = frag_start + frag_len_actual
    r2_start = r2_end - rev_len

    wrapped = (frag_start + frag_len_actual) > seq_length

    return r1_start, r1_end, r2_start, r2_end, wrapped


def format_header(read_id, source_tsv, contig=None, start=None, end=None,
                  strand=None, wrapped=False):
    """
    when a source TSV is being written, the header is just the unique read id so
    read names stay clean for downstream tools. otherwise we're writing the full provenance 
    in the comment field (after a space) as key=value pairs
    """
    if source_tsv:
        return read_id

    comment = f"contig={contig};start={start};end={end};strand={strand}"
    if wrapped:
        comment += ";wrapped=true"

    return f"{read_id} {comment}"


def write_source_row(source_tsv, read_id, source_fasta, contig, start, end,
                     strand, wrapped):
    """ writes one per-read provenance row to the source TSV """

    source_tsv.write(f"{read_id}\t{source_fasta}\t{contig}\t{start}\t{end}\t"
                     f"{strand}\t{str(wrapped).lower()}\n")


def reverse_complement(seq):

    return seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


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
        _status(args, f"\n    Note: Paired-end mode requires an even number of reads. Rounding up to {args.num_reads}.")

    forward_reads_file = f"{args.output_prefix}_R1.fastq"
    reverse_reads_file = f"{args.output_prefix}_R2.fastq"

    pct = args.fragment_size_range / 100
    min_frag = max(1, int(args.fragment_size * (1 - pct)))
    max_frag = int(args.fragment_size * (1 + pct))

    source_tsv = open(f"{args.output_prefix}-read-sources.tsv", "w") if args.source_tsv else None
    if source_tsv:
        source_tsv.write("read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n")

    with open(forward_reads_file, 'w') as fw, open(reverse_reads_file, 'w') as rw:

        print("")

        num_fragments = args.num_reads // 2
        reads_remaining = num_fragments
        remainder = 0.0
        read_count = 0

        if getattr(args, "quiet", False):
            desc = "    Progress"
            ncols = 78
        else:
            desc = "    Generating reads from fasta file(s)"
            ncols = 90

        pbar = tqdm(total=len(args.input_fastas), desc = desc, unit="file", ncols=ncols)

        for fasta_file in args.input_fastas:

            total_length = sum(len(record.seq) for record in parse_fasta(fasta_file))

            for record in parse_fasta(fasta_file):
                seq_id = record.id
                sequence = str(record.seq).upper()
                seq_length = len(sequence)

                exact = (seq_length / total_length) * num_fragments * proportions[fasta_file] + remainder
                entry_reads = min(round(exact), reads_remaining)
                remainder = exact - round(exact)

                # generate paired-end reads
                for _ in range(entry_reads):

                    frag_len = min(max(min_frag, round(random.gauss(args.fragment_size, args.fragment_size * pct / 3))), max_frag)
                    frag_len = min(frag_len, seq_length)

                    fragment, start = extract_subsequence(sequence, seq_length, frag_len, args.circularize, args.include_Ns)

                    # real libraries draw the fragment from either strand ~50/50;
                    # R2 is always the opposite orientation of its R1 mate. strands are
                    # reported relative to the input reference, so a minus-strand fragment
                    # also swaps which reference span each mate covers
                    if random.random() < 0.5:
                        fragment_strand = "+"
                    else:
                        fragment_strand = "-"
                        fragment = reverse_complement(fragment)

                    forward_read = fragment[:args.read_length]
                    reverse_read = fragment[-args.read_length:][::-1].translate(str.maketrans("ACGT", "TGCA"))

                    fwd_quality_scores = "I" * len(forward_read)
                    rev_quality_scores = "I" * len(reverse_read)

                    r1_start, r1_end, r2_start, r2_end, wrapped = compute_mate_coords(
                        start, fragment, len(forward_read), len(reverse_read), seq_length)

                    if fragment_strand == "-":
                        r1_start, r1_end, r2_start, r2_end = r2_start, r2_end, r1_start, r1_end

                    r1_strand = fragment_strand
                    r2_strand = "-" if fragment_strand == "+" else "+"

                    read_count += 1
                    read_id = f"r{read_count}"

                    fw.write(f"@{format_header(read_id + '/1', source_tsv, seq_id, r1_start, r1_end, r1_strand, wrapped)}\n")
                    fw.write(f"{forward_read}\n")
                    fw.write(f"+\n")
                    fw.write(f"{fwd_quality_scores}\n")

                    rw.write(f"@{format_header(read_id + '/2', source_tsv, seq_id, r2_start, r2_end, r2_strand, wrapped)}\n")
                    rw.write(f"{reverse_read}\n")
                    rw.write(f"+\n")
                    rw.write(f"{rev_quality_scores}\n")

                    if source_tsv:
                        write_source_row(source_tsv, read_id + "/1", fasta_file, seq_id, r1_start, r1_end, r1_strand, wrapped)
                        write_source_row(source_tsv, read_id + "/2", fasta_file, seq_id, r2_start, r2_end, r2_strand, wrapped)

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

                frag_len = min(max(min_frag, round(random.gauss(args.fragment_size, args.fragment_size * pct / 3))), max_frag)
                frag_len = min(frag_len, seq_length)

                fragment, start = extract_subsequence(sequence, seq_length, frag_len, args.circularize, args.include_Ns)

                # real libraries draw the fragment from either strand ~50/50;
                # R2 is always the opposite orientation of its R1 mate. strands are
                # reported relative to the input reference, so a minus-strand fragment
                # also swaps which reference span each mate covers
                if random.random() < 0.5:
                    fragment_strand = "+"
                else:
                    fragment_strand = "-"
                    fragment = reverse_complement(fragment)

                forward_read = fragment[:args.read_length]
                reverse_read = fragment[-args.read_length:][::-1].translate(str.maketrans("ACGT", "TGCA"))
                fwd_quality_scores = "I" * len(forward_read)
                rev_quality_scores = "I" * len(reverse_read)

                r1_start, r1_end, r2_start, r2_end, wrapped = compute_mate_coords(
                    start, fragment, len(forward_read), len(reverse_read), seq_length)

                if fragment_strand == "-":
                    r1_start, r1_end, r2_start, r2_end = r2_start, r2_end, r1_start, r1_end

                r1_strand = fragment_strand
                r2_strand = "-" if fragment_strand == "+" else "+"

                read_count += 1
                read_id = f"r{read_count}"

                fw.write(f"@{format_header(read_id + '/1', source_tsv, seq_id, r1_start, r1_end, r1_strand, wrapped)}\n")
                fw.write(f"{forward_read}\n")
                fw.write(f"+\n")
                fw.write(f"{fwd_quality_scores}\n")

                rw.write(f"@{format_header(read_id + '/2', source_tsv, seq_id, r2_start, r2_end, r2_strand, wrapped)}\n")
                rw.write(f"{reverse_read}\n")
                rw.write(f"+\n")
                rw.write(f"{rev_quality_scores}\n")

                if source_tsv:
                    write_source_row(source_tsv, read_id + "/1", fasta_file, seq_id, r1_start, r1_end, r1_strand, wrapped)
                    write_source_row(source_tsv, read_id + "/2", fasta_file, seq_id, r2_start, r2_end, r2_strand, wrapped)

        if source_tsv:
            source_tsv.close()


def gen_single_reads(args, proportions):

    reads_file = f"{args.output_prefix}.fastq"

    source_tsv = open(f"{args.output_prefix}-read-sources.tsv", "w") if args.source_tsv else None
    if source_tsv:
        source_tsv.write("read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n")

    with open(reads_file, 'w') as fw:

        print("")

        if args.type == "long":
            pct = args.long_read_length_range / 100
            min_len = max(1, int(args.read_length * (1 - pct)))
            max_len = int(args.read_length * (1 + pct))

        reads_remaining = args.num_reads
        remainder = 0.0
        read_count = 0

        if getattr(args, "quiet", False):
            desc = "    Progress"
            ncols = 78
        else:
            desc = "    Generating reads from fasta file(s)"
            ncols = 90

        pbar = tqdm(total=len(args.input_fastas), desc = desc, unit="file", ncols=ncols)

        for fasta_file in args.input_fastas:

            total_length = sum(len(record.seq) for record in parse_fasta(fasta_file))

            for record in parse_fasta(fasta_file):
                seq_id = record.id
                sequence = str(record.seq)
                seq_length = len(sequence)

                exact = (seq_length / total_length) * args.num_reads * proportions[fasta_file] + remainder
                entry_reads = min(round(exact), reads_remaining)
                remainder = exact - round(exact)

                for _ in range(entry_reads):

                    if args.type == "long":
                        read_len = min(max(min_len, round(random.gauss(args.read_length, args.read_length * pct / 3))), max_len)
                        read_len = min(read_len, seq_length)
                    else:
                        read_len = min(args.read_length, seq_length)

                    read, start = extract_subsequence(sequence, seq_length, read_len, args.circularize, args.include_Ns)

                    read_end = start + len(read)
                    wrapped = read_end > seq_length

                    # real single-end / long-read libraries sequence both strands ~50/50
                    if random.random() < 0.5:
                        strand = "+"
                    else:
                        strand = "-"
                        read = reverse_complement(read)

                    quality_scores = "I" * len(read)

                    read_count += 1
                    read_id = f"r{read_count}"
                    fw.write(f"@{format_header(read_id, source_tsv, seq_id, start, read_end, strand, wrapped)}\n")
                    fw.write(f"{read}\n")
                    fw.write(f"+\n")
                    fw.write(f"{quality_scores}\n")

                    if source_tsv:
                        write_source_row(source_tsv, read_id, fasta_file, seq_id, start, read_end, strand, wrapped)

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
                    read_len = min(max(min_len, round(random.gauss(args.read_length, args.read_length * pct / 3))), max_len)
                    read_len = min(read_len, seq_length)
                else:
                    read_len = min(args.read_length, seq_length)

                read, start = extract_subsequence(sequence, seq_length, read_len, args.circularize, args.include_Ns)

                read_end = start + len(read)
                wrapped = read_end > seq_length

                if random.random() < 0.5:
                    strand = "+"
                else:
                    strand = "-"
                    read = reverse_complement(read)

                quality_scores = "I" * len(read)

                read_count += 1
                read_id = f"r{read_count}"
                fw.write(f"@{format_header(read_id, source_tsv, seq_id, start, read_end, strand, wrapped)}\n")
                fw.write(f"{read}\n")
                fw.write(f"+\n")
                fw.write(f"{quality_scores}\n")

                if source_tsv:
                    write_source_row(source_tsv, read_id, fasta_file, seq_id, start, read_end, strand, wrapped)

        if source_tsv:
            source_tsv.close()


def compress_with_pigz(output_prefix, read_type="paired-end", quiet=False):

    if read_type in ("single-end", "long"):
        reads_file = f"{output_prefix}.fastq"
        if not quiet:
            print("\n    Compressing output FASTQ file...")
        try:
            subprocess.run(["pigz", "-f", reads_file], check = True)
            if not quiet:
                print(f"\n    Compressed file written: {reads_file}.gz\n")
        except FileNotFoundError:
            print("pigz not found. You're on your own for compression!")
        except subprocess.CalledProcessError as e:
            print(f"Error during compression: {e}")
    else:
        forward_reads_file = f"{output_prefix}_R1.fastq"
        reverse_reads_file = f"{output_prefix}_R2.fastq"

        if not quiet:
            print("\n    Compressing output fastq files...")

        try:
            subprocess.run(["pigz", "-f", forward_reads_file], check = True)
            subprocess.run(["pigz", "-f", reverse_reads_file], check = True)
            if not quiet:
                print(f"\n    Compressed files written: {forward_reads_file}.gz, {reverse_reads_file}.gz\n")
        except FileNotFoundError:
            print("pigz not found. You're on your own for compression!")
        except subprocess.CalledProcessError as e:
            print(f"Error during compression: {e}")
