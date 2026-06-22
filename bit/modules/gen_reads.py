from Bio import SeqIO # type: ignore
import gzip
import random
import os
import sys
import shutil
import tempfile
import multiprocessing
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
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

    # gen-metagenome already measures every genome; if it passed those sizes
    # through (keyed by the fasta path), we reuse them
    known_sizes = getattr(args, "genome_sizes", None) or {}

    for fasta_file in args.input_fastas:
        genome_size = known_sizes.get(fasta_file)
        if genome_size is None:
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


def extract_subsequence(sequence, seq_length, subseq_len, circularize, include_Ns=False, max_attempts=100, rng=None):

    # rng is a random.Random instance for per-file reproducibility; falls back to
    # the module-global random for any legacy callers
    if rng is None:
        rng = random

    for _ in range(max_attempts):

        if circularize:
            start = rng.randint(0, seq_length - 1)
            end = start + subseq_len
            if end <= seq_length:
                subseq = sequence[start:end]
            else:
                subseq = sequence[start:] + sequence[: end - seq_length]
        elif subseq_len >= seq_length:
            # when the requested length meets or exceeds the contig length,
            # pick any start position and truncate at the contig end
            start = rng.randint(0, seq_length - 1)
            subseq = sequence[start:]
        else:
            start = rng.randint(0, seq_length - subseq_len)
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


def apportion_units(input_fastas, total_units, proportions):
    """
    Apportion an integer number of `units` across files by their proportions using
    the largest-remainder (Hamilton) method. `units` are fragments for paired-end
    and reads for single-end/long. The per-file counts sum exactly to total_units
    and are independent of contig structure and execution order, which is what makes
    per-file generation parallelizable. Returns {fasta: units} in input order.

    Note: per-file counts can differ by at most 1 from the previous serial
    per-contig remainder logic for multi-file inputs; single-file inputs are
    unchanged. The grand total always equals total_units.
    """
    # proportions may not sum to exactly 1 (callers/tests can pass raw weights);
    # normalize so the largest-remainder distribution is well-defined. if they sum
    # to <= 0, fall back to equal weighting
    weight_total = sum(proportions[f] for f in input_fastas)
    if weight_total <= 0:
        norm = {f: 1.0 / len(input_fastas) for f in input_fastas}
    else:
        norm = {f: proportions[f] / weight_total for f in input_fastas}

    raw = {f: total_units * norm[f] for f in input_fastas}
    floor_counts = {f: int(raw[f]) for f in input_fastas}
    leftover = total_units - sum(floor_counts.values())

    out = dict(floor_counts)
    if leftover > 0:
        # rank by fractional remainder (desc); ties broken by input order for determinism
        order = sorted(
            range(len(input_fastas)),
            key=lambda i: (raw[input_fastas[i]] - floor_counts[input_fastas[i]], -i),
            reverse=True,
        )
        for k in range(leftover):
            out[input_fastas[order[k % len(order)]]] += 1
    return out


def compute_id_offsets(input_fastas, per_file_units):
    """
    Return {fasta: starting_read_index}, where the index is the number of base
    read ids already consumed by earlier files. The base read-id counter advances
    once per unit in every mode: in paired-end a single fragment (unit) produces
    r{n}/1 and r{n}/2 sharing one base number, and in single/long each read is its
    own unit. Workers emit r{offset+1}, r{offset+2}, ... so the concatenated output
    keeps the global r{N} numbering scheme.
    """
    offsets = {}
    running = 0
    for f in input_fastas:
        offsets[f] = running
        running += per_file_units[f]
    return offsets


def _local_seed(base_seed, file_index):
    """
    Per-file seed so each worker's RNG is independent of execution order yet fully
    reproducible. When no base seed is given, return None so workers use system
    entropy (matching the previous unseeded, non-reproducible behavior).
    """
    if base_seed is None:
        return None
    return base_seed + file_index


def distribute_units_across_contigs(contig_lengths, file_budget):
    """
    Split a file's unit budget across its own contigs proportional to length, using
    a remainder carry that is local to this file. Any leftover from rounding is added
    to the last contig. Returns a list of per-contig unit counts (same order as input).
    """
    total_length = sum(contig_lengths)
    if total_length == 0 or file_budget <= 0:
        return [0] * len(contig_lengths)

    counts = []
    remaining = file_budget
    remainder = 0.0
    for seq_length in contig_lengths:
        exact = (seq_length / total_length) * file_budget + remainder
        n = min(round(exact), remaining)
        remainder = exact - round(exact)
        counts.append(n)
        remaining -= n
        if remaining <= 0:
            # pad the rest with zeros so the list length matches contig count
            counts.extend([0] * (len(contig_lengths) - len(counts)))
            break

    # leftover from rounding goes on the last contig that received any reads,
    # falling back to the final contig
    if remaining > 0:
        for i in range(len(counts) - 1, -1, -1):
            if counts[i] > 0:
                counts[i] += remaining
                break
        else:
            counts[-1] += remaining

    return counts


def gen_reads(args, proportions):

    jobs = max(1, int(getattr(args, "jobs", 1) or 1))

    # paired-end works in fragments (each yields 2 reads); single/long in reads
    if args.type in ("single-end", "long"):
        total_units = args.num_reads
    else:
        if args.num_reads % 2 != 0:
            args.num_reads += 1
            _status(args, f"\n    Note: Paired-end mode requires an even number of reads. Rounding up to {args.num_reads}.")
        total_units = args.num_reads // 2

    per_file_units = apportion_units(args.input_fastas, total_units, proportions)
    id_offsets = compute_id_offsets(args.input_fastas, per_file_units)

    if jobs == 1 or len(args.input_fastas) == 1:
        _gen_reads_serial(args, per_file_units, id_offsets)
    else:
        _gen_reads_parallel(args, per_file_units, id_offsets, jobs)


def _open_pair_writers(prefix):
    """ open R1/R2 fastq writers for a given output prefix """
    fw = open(f"{prefix}_R1.fastq", "w")
    rw = open(f"{prefix}_R2.fastq", "w")
    return fw, rw


def _open_source_tsv(prefix, write_header):
    """ open a per-read provenance TSV, optionally writing the header row """
    fh = open(f"{prefix}-read-sources.tsv", "w")
    if write_header:
        fh.write("read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n")
    return fh


def gen_one_file_paired(args, fasta_file, file_budget, id_offset, local_seed,
                        fw, rw, source_tsv):
    """
    Generate `file_budget` fragments (2 reads each) from a single fasta, writing to
    the already-open R1/R2 handles (and optional source TSV). Read ids start at
    id_offset+1 so concatenated output keeps global r{N} numbering. Uses a private
    RNG seeded with local_seed for order-independent reproducibility.

    Returns the number of reads written (== file_budget * 2).
    """
    rng = random.Random(local_seed)

    pct = args.fragment_size_range / 100
    min_frag = max(1, int(args.fragment_size * (1 - pct)))
    max_frag = int(args.fragment_size * (1 + pct))

    records = list(parse_fasta(fasta_file))
    contig_lengths = [len(r.seq) for r in records]
    per_contig = distribute_units_across_contigs(contig_lengths, file_budget)

    read_count = id_offset
    reads_written = 0

    for record, entry_reads in zip(records, per_contig):
        if entry_reads <= 0:
            continue
        seq_id = record.id
        sequence = str(record.seq).upper()
        seq_length = len(sequence)

        for _ in range(entry_reads):

            frag_len = min(max(min_frag, round(rng.gauss(args.fragment_size, args.fragment_size * pct / 3))), max_frag)
            frag_len = min(frag_len, seq_length)

            fragment, start = extract_subsequence(sequence, seq_length, frag_len, args.circularize, args.include_Ns, rng=rng)

            # real libraries draw the fragment from either strand ~50/50;
            # R2 is always the opposite orientation of its R1 mate. strands are
            # reported relative to the input reference, so a minus-strand fragment
            # also swaps which reference span each mate covers
            if rng.random() < 0.5:
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

            reads_written += 2

    return reads_written


def gen_one_file_single(args, fasta_file, file_budget, id_offset, local_seed,
                        fw, source_tsv):
    """
    Generate `file_budget` reads (single-end or long) from a single fasta, writing to
    the already-open handle (and optional source TSV). Read ids start at id_offset+1.
    Uses a private RNG seeded with local_seed.

    Returns the number of reads written (== file_budget).
    """
    rng = random.Random(local_seed)

    if args.type == "long":
        pct = args.long_read_length_range / 100
        min_len = max(1, int(args.read_length * (1 - pct)))
        max_len = int(args.read_length * (1 + pct))

    records = list(parse_fasta(fasta_file))
    contig_lengths = [len(r.seq) for r in records]
    per_contig = distribute_units_across_contigs(contig_lengths, file_budget)

    read_count = id_offset
    reads_written = 0

    for record, entry_reads in zip(records, per_contig):
        if entry_reads <= 0:
            continue
        seq_id = record.id
        sequence = str(record.seq)
        seq_length = len(sequence)

        for _ in range(entry_reads):

            if args.type == "long":
                read_len = min(max(min_len, round(rng.gauss(args.read_length, args.read_length * pct / 3))), max_len)
                read_len = min(read_len, seq_length)
            else:
                read_len = min(args.read_length, seq_length)

            read, start = extract_subsequence(sequence, seq_length, read_len, args.circularize, args.include_Ns, rng=rng)

            read_end = start + len(read)
            wrapped = read_end > seq_length

            # real single-end / long-read libraries sequence both strands ~50/50
            if rng.random() < 0.5:
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

            reads_written += 1

    return reads_written


def _progress_bar(args, total):
    if getattr(args, "quiet", False):
        desc, ncols = "    Progress", 78
    else:
        desc, ncols = "    Generating reads from fasta file(s)", 90
    return tqdm(total=total, desc=desc, unit=" file", ncols=ncols)


def _gen_reads_serial(args, per_file_units, id_offsets):
    """ single-process path: generate each file in turn, writing to the final files """

    paired = args.type not in ("single-end", "long")

    if paired:
        fw, rw = _open_pair_writers(args.output_prefix)
    else:
        fw = open(f"{args.output_prefix}.fastq", "w")
        rw = None

    source_tsv = _open_source_tsv(args.output_prefix, write_header=True) if args.source_tsv else None

    print("")
    pbar = _progress_bar(args, len(args.input_fastas))

    try:
        for idx, fasta_file in enumerate(args.input_fastas):
            local_seed = _local_seed(args.seed, idx)
            if paired:
                gen_one_file_paired(args, fasta_file, per_file_units[fasta_file],
                                    id_offsets[fasta_file], local_seed, fw, rw, source_tsv)
            else:
                gen_one_file_single(args, fasta_file, per_file_units[fasta_file],
                                    id_offsets[fasta_file], local_seed, fw, source_tsv)
            pbar.update(1)
    finally:
        pbar.update(pbar.total - pbar.n)
        pbar.close()
        fw.close()
        if rw is not None:
            rw.close()
        if source_tsv:
            source_tsv.close()


def _worker_generate_file(spec):
    """
    Picklable worker run in a separate process. Generates reads for one fasta into
    temp files under spec['tmp_dir'] and returns the temp paths so the parent can
    concatenate them in input order. No shared state, no shared handles.
    """
    args = spec["args"]
    fasta_file = spec["fasta_file"]
    idx = spec["index"]
    file_budget = spec["budget"]
    id_offset = spec["id_offset"]
    tmp_dir = spec["tmp_dir"]

    paired = args.type not in ("single-end", "long")
    local_seed = _local_seed(args.seed, idx)

    tmp_prefix = os.path.join(tmp_dir, f"part_{idx:06d}")

    source_tsv = _open_source_tsv(tmp_prefix, write_header=False) if args.source_tsv else None

    if paired:
        fw, rw = _open_pair_writers(tmp_prefix)
        try:
            gen_one_file_paired(args, fasta_file, file_budget, id_offset, local_seed,
                                fw, rw, source_tsv)
        finally:
            fw.close()
            rw.close()
            if source_tsv:
                source_tsv.close()
        return {
            "index": idx,
            "r1": f"{tmp_prefix}_R1.fastq",
            "r2": f"{tmp_prefix}_R2.fastq",
            "source": f"{tmp_prefix}-read-sources.tsv" if args.source_tsv else None,
        }
    else:
        fw = open(f"{tmp_prefix}.fastq", "w")
        try:
            gen_one_file_single(args, fasta_file, file_budget, id_offset, local_seed,
                                fw, source_tsv)
        finally:
            fw.close()
            if source_tsv:
                source_tsv.close()
        return {
            "index": idx,
            "reads": f"{tmp_prefix}.fastq",
            "source": f"{tmp_prefix}-read-sources.tsv" if args.source_tsv else None,
        }


def _concat_files(parts, dest):
    """ append each part file to dest (binary copy, order as given) """
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as src:
                shutil.copyfileobj(src, out)


def _gen_reads_parallel(args, per_file_units, id_offsets, jobs):
    """
    multi-process path: each file is generated independently into temp files, then
    concatenated in input order. Deterministic given a seed (per-file seeds), and
    independent of which worker finishes first.
    """

    paired = args.type not in ("single-end", "long")
    jobs = min(jobs, len(args.input_fastas))

    print("")
    pbar = _progress_bar(args, len(args.input_fastas))

    tmp_dir = tempfile.mkdtemp(prefix="bit-gen-reads-")
    results = {}

    try:
        specs = []
        for idx, fasta_file in enumerate(args.input_fastas):
            specs.append({
                "args": args,
                "fasta_file": fasta_file,
                "index": idx,
                "budget": per_file_units[fasta_file],
                "id_offset": id_offsets[fasta_file],
                "tmp_dir": tmp_dir,
            })

        # use 'spawn' rather than the default 'fork': gen-reads can be called from
        # multi-threaded contexts (e.g. gen-metagenome), where forking risks deadlocks
        # and emits a DeprecationWarning on Python 3.12+. spawn re-imports this module
        # in each worker, which is fine since the worker and its helpers are top-level.
        mp_ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=jobs, mp_context=mp_ctx) as ex:
            futures = {ex.submit(_worker_generate_file, s): s["index"] for s in specs}
            for fut in as_completed(futures):
                res = fut.result()
                results[res["index"]] = res
                pbar.update(1)

        pbar.update(pbar.total - pbar.n)
        pbar.close()

        ordered = [results[i] for i in range(len(args.input_fastas))]

        if paired:
            _concat_files([r["r1"] for r in ordered], f"{args.output_prefix}_R1.fastq")
            _concat_files([r["r2"] for r in ordered], f"{args.output_prefix}_R2.fastq")
        else:
            _concat_files([r["reads"] for r in ordered], f"{args.output_prefix}.fastq")

        if args.source_tsv:
            dest = f"{args.output_prefix}-read-sources.tsv"
            with open(dest, "w") as out:
                out.write("read_id\tsource_fasta\tcontig\tstart\tend\tstrand\twrapped\n")
                for r in ordered:
                    with open(r["source"], "r") as src:
                        shutil.copyfileobj(src, out)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


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
