import random
import datetime
from Bio import SeqIO
from collections import OrderedDict
from bit.modules.general import report_message, notify_premature_exit



def add_insertion(args):

    seed = set_seed(args.seed)

    contigs = parse_fasta(args.input_fasta)
    insertion_seq = parse_insertion_fasta(args.insertion_fasta)

    if args.position is not None:
        validate_position(args.position, contigs)

    if args.back_to_back and args.num_insertions > 1 and args.position is None:
        contig_name, pos = pick_random_position(contigs)
        combined = insertion_seq * args.num_insertions
        contigs[contig_name] = contigs[contig_name][:pos] + combined + contigs[contig_name][pos:]
        positions_used = [(contig_name, pos)] * args.num_insertions

    elif args.back_to_back and args.num_insertions > 1 and args.position is not None:
        contig_name, pos = parse_position_string(args.position)
        combined = insertion_seq * args.num_insertions
        contigs[contig_name] = contigs[contig_name][:pos] + combined + contigs[contig_name][pos:]
        positions_used = [(contig_name, pos)] * args.num_insertions

    elif args.position is not None:
        contig_name, pos = parse_position_string(args.position)
        # first insertion at the specified position
        contigs[contig_name] = contigs[contig_name][:pos] + insertion_seq + contigs[contig_name][pos:]
        positions_used = [(contig_name, pos)]
        # any additional insertions go to random positions
        for _ in range(args.num_insertions - 1):
            contig_name, pos = pick_random_position(contigs)
            contigs[contig_name] = contigs[contig_name][:pos] + insertion_seq + contigs[contig_name][pos:]
            positions_used.append((contig_name, pos))

    elif args.num_insertions == 1:
        contig_name, pos = pick_random_position(contigs)
        contigs[contig_name] = contigs[contig_name][:pos] + insertion_seq + contigs[contig_name][pos:]
        positions_used = [(contig_name, pos)]

    else:
        # multiple random insertions, spaced out
        positions_used = []
        for _ in range(args.num_insertions):
            contig_name, pos = pick_random_position(contigs)
            contigs[contig_name] = contigs[contig_name][:pos] + insertion_seq + contigs[contig_name][pos:]
            positions_used.append((contig_name, pos))

    write_fasta(contigs, args.output_fasta)

    return {
        "seed": seed,
        "positions": positions_used,
        "insertion_length": len(insertion_seq),
        "num_insertions": args.num_insertions,
    }


def set_seed(input_seed):
    if input_seed is not None:
        seed = input_seed
    else:
        time = datetime.datetime.now()
        seed = time.hour * 10000 + time.minute * 100 + time.second
    random.seed(seed)
    return seed


def parse_fasta(fasta_path):

    contigs = OrderedDict()
    for record in SeqIO.parse(fasta_path, "fasta"):
        contigs[record.id] = str(record.seq)
    if not contigs:
        report_message("No sequences were found in the input fasta.")
        notify_premature_exit()
    return contigs


def parse_insertion_fasta(fasta_path):

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) == 0:
        report_message("No sequences were found in the insertion fasta.")
        notify_premature_exit()
    if len(records) > 1:
        report_message("Multiple sequences were found in the insertion fasta, but currently bit can only work with one. \
                       Cut it down to one or reach out demanding support for multiple!")
        notify_premature_exit()
    return str(records[0].seq)


def parse_position_string(position_str):

    parts = position_str.rsplit(":", 1)
    if len(parts) != 2:
        report_message(f"Invalid position format: '{position_str}'")
        report_message("Expected 'contig_name:position'")
        notify_premature_exit()
    contig_name = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        report_message(
            f"Position value '{parts[1]}' is not a valid integer."
        )
        notify_premature_exit()
    if pos < 0:
        report_message(f"Position must be >= 0, got {pos}.")
        notify_premature_exit()
    return contig_name, pos


def validate_position(position_str, contigs):

    contig_name, pos = parse_position_string(position_str)
    if contig_name not in contigs:
        report_message(f"'{contig_name}' not found in input fasta.")
        notify_premature_exit()
    if pos > len(contigs[contig_name]):
        report_message(f"Position {pos} is beyond the length of '{contig_name}' "
                       f"(length: {len(contigs[contig_name])}).")
        notify_premature_exit()


def pick_random_position(contigs):
    """pick a random contig (weighted by length) and a random position within it"""
    contig_names = list(contigs.keys())
    lengths = [len(contigs[name]) for name in contig_names]
    contig_name = random.choices(contig_names, weights=lengths, k=1)[0]
    pos = random.randint(0, len(contigs[contig_name]))
    return contig_name, pos


def write_fasta(contigs, output_path, line_width=60):

    with open(output_path, "w") as f:
        for header, seq in contigs.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")
