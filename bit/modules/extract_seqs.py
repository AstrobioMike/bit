import os
import edlib # type: ignore
from pybedtools import BedTool # type: ignore
from bit.modules.general import report_message
from bit.modules.seqs import revcomp, read_fasta

def extract_seqs_by_coords(args):
    coordinates_file = BedTool(args.bed_file)
    fasta = BedTool(args.input_fasta)
    seq = coordinates_file.sequence(fi = fasta)

    content = open(seq.seqfn).read()

    if not content.strip():
        report_message("No sequences were extracted based on the provided coordinates.", trailing_newline=True)
        return

    with open(args.output_fasta, "w") as out_fasta:
        out_fasta.write(content)
    report_message(f"Extracted sequences based on the provided coordinates were written to:")
    report_message(args.output_fasta, color="green", initial_indent="    ", leading_newline=False, trailing_newline=True)


def extract_seqs_by_primers(args):
    in_fasta = args.input_fasta
    fwd = args.forward_primer.upper().strip()
    rev = args.reverse_primer.upper().strip()

    hit_count = 0

    with open(args.output_fasta, "w") as out_fasta:

        for header, seq in read_fasta(in_fasta):
            seq = seq.upper()
            amplicons = find_amplicons(seq, fwd, rev, args.max_mismatches)

            for _, (left_label, right_label, left_start, right_end, amplicon, length) in enumerate(amplicons):

                out_header = f"{header}|{left_label}-to-{right_label}|{left_start}-{right_end}|{length}"
                out_seq = amplicon

                out_fasta.write(f">{out_header}\n")
                out_fasta.write(f"{out_seq}\n")
                hit_count += 1

    if hit_count == 0:
        os.remove(args.output_fasta)
        report_message("No sequences were found based on the provided primers.", trailing_newline=True)
    else:
        report_message(f"Extracted {hit_count} sequence(s) based on the provided primers, written to:")
        report_message(args.output_fasta, color="green", initial_indent="    ", leading_newline=False, trailing_newline=True)


def find_all_primer_hits(seq, fwd, rev, max_mismatches = 0):

    primers = [
        ("fwd", fwd),
        ("rev", rev),
        ("fwd_rc", revcomp(fwd)),
        ("rev_rc", revcomp(rev)),
    ]

    hits = []

    for label, primer in primers:

        result = edlib.align(
            primer,
            seq,
            mode = "HW",            # search within sequence
            task = "locations",     # return match positions
            k = max_mismatches,
        )

        if result["editDistance"] == -1:
            continue

        for start, end_inclusive in result["locations"]:
            end = end_inclusive + 1

            hits.append(
                (
                    label,
                    start,
                    end,
                    seq[start:end]
                )
            )

    hits.sort(key = lambda x: x[1])

    # remove potential duplicates
    seen = set()
    unique_hits = []

    for hit in hits:
        key = (hit[1], hit[2])  # using start and end positions as the key
        if key not in seen:
            seen.add(key)
            unique_hits.append(hit)

    return unique_hits


def is_forward_label(label):
    return label in {"fwd", "fwd_rc"}


def is_reverse_label(label):
    return label in {"rev", "rev_rc"}


def find_amplicons(seq, fwd, rev, max_mismatches = 0):
    hits = find_all_primer_hits(seq, fwd, rev, max_mismatches)
    results = []

    for i, left in enumerate(hits):
        left_label, left_start, left_end, left_primer = left

        for right in hits[i + 1:]:
            right_label, right_start, right_end, right_primer = right

            if right_start < left_end:
                continue

            left_is_forward = is_forward_label(left_label)
            right_is_forward = is_forward_label(right_label)

            if left_is_forward == right_is_forward:
                continue

            amplicon = seq[left_start:right_end]
            length = len(amplicon)

            results.append(
                (
                    left_label,
                    right_label,
                    left_start,
                    right_end,
                    amplicon,
                    length
                )
            )

    return results
