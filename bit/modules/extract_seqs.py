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
    report_message(f"Extracted sequences based on the provided coordinates were written to:", color="none")
    report_message(args.output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)


def extract_seqs_by_headers(args):

    if args.headers:
        headers_of_interest = set(line.strip() for line in args.headers)
    else:
        headers_of_interest = set(line.strip() for line in open(args.file_with_headers, "r"))

    seqs_pulled = 0
    num_targets = len(headers_of_interest)

    with open(args.output_fasta, "w") as out_fasta:

        if not args.inverse:

            for header, seq in read_fasta(args.input_fasta):
                if header in headers_of_interest:
                    seqs_pulled += 1
                    out_fasta.write(f">{header}\n")
                    out_fasta.write(f"{seq}\n")

            report_by_header_results(seqs_pulled, num_targets, args.output_fasta)

        else:

            for header, seq in read_fasta(args.input_fasta):
                if header not in headers_of_interest:
                    seqs_pulled += 1
                    out_fasta.write(f">{header}\n")
                    out_fasta.write(f"{seq}\n")

            report_by_header_inverse_results(seqs_pulled, args.output_fasta)


def report_by_header_results(seqs_pulled, num_targets, output_fasta):

    if seqs_pulled == 0:
        report_message("No sequences were found based on the provided headers.", initial_indent="    ",
                        subsequent_indent="    ", trailing_newline=True)
        os.remove(output_fasta)

    if seqs_pulled == num_targets:
        report_message(f"Extracted all {seqs_pulled} target sequence(s), written to:\n", color="none")
        report_message(output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)

    if seqs_pulled < num_targets:
        report_message(f"Extracted {seqs_pulled} out of {num_targets} target sequence(s) based on the provided headers, written to:\n", color="none")
        report_message(output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)


def report_by_header_inverse_results(seqs_pulled, output_fasta):

    if seqs_pulled == 0:
        report_message("No sequences were extracted based on the provided headers with the --inverse flag.",
                        initial_indent="    ", subsequent_indent="    ", trailing_newline=True)
        os.remove(output_fasta)

    if seqs_pulled > 0:
        report_message(f"Extracted {seqs_pulled} sequence(s) based on the provided headers with the --inverse flag, written to:\n", color="none")
        report_message(output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)


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
        report_message(f"Extracted {hit_count} sequence(s) based on the provided primers, written to:", color="none")
        report_message(args.output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)


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
