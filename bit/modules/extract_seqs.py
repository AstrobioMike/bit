import os
import re
import edlib # type: ignore
from pybedtools import BedTool # type: ignore
from bit.modules.general import report_message, notify_premature_exit
from bit.modules.seqs import revcomp, read_fasta


def _normalize_header(header):
    return header.strip().lstrip(">")


def extract_seqs_by_coords(args):
    coordinates = _build_coords_bedtool(args.bed)
    fasta = BedTool(args.input_fasta)
    seq = coordinates.sequence(fi = fasta)

    content = open(seq.seqfn).read()

    if not content.strip():
        report_message("No sequences were extracted based on the provided coordinates.", trailing_newline=True)
        return

    with open(args.output_fasta, "w") as out_fasta:
        out_fasta.write(content)
    report_message(f"Extracted sequences based on the provided coordinates were written to:", color="none")
    report_message(args.output_fasta, color="yellow", initial_indent="    ", leading_newline=False, trailing_newline=True)


def _build_coords_bedtool(bed):
    """
    accepts either a single bed-file path, or three strings (contig, start, end)
    specifying one region inline. returns a BedTool
    """
    if len(bed) == 1:
        path = bed[0]
        if not os.path.isfile(path):
            report_message(f"Provided coordinate input '{path}' is not an existing file. "
                           f"Provide a bed file, or three values: contig start end.",
                           color="yellow", initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()
        return BedTool(path)

    if len(bed) == 3:
        contig, start, end = bed
        try:
            start, end = int(start), int(end)
        except ValueError:
            report_message(f"Inline coordinates must be integers, but got start='{bed[1]}' and end='{bed[2]}'.",
                           color="yellow", initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()
        if start >= end:
            report_message(f"Inline start ({start}) must be less than end ({end}).",
                           color="yellow", initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()
        return BedTool(f"{contig}\t{start}\t{end}\n", from_string=True)

    report_message(f"Provide either a single bed file, or exactly three values: contig start end ",
                   color="yellow", initial_indent="    ", subsequent_indent="    ")
    notify_premature_exit()


def extract_seqs_by_headers(args):

    if len(args.headers) == 1 and os.path.isfile(args.headers[0]):
            with open(args.headers[0], "r") as in_headers:
                headers_of_interest = set(_normalize_header(line) for line in in_headers)
    else:
        if len(args.headers) > 1:
            files = [entry for entry in args.headers if os.path.isfile(entry)]
            if files:
                report_message(f"Note: We're treating all of -H/--headers as literal header strings, but "
                            f"{', '.join(files)} "
                            f"{'is an existing file' if len(files) == 1 else 'are existing files!'} "
                            f"File contents are only read when a single file path "
                            f"is the sole argument to -H. You may want to re-evaluate and re-run.",
                            initial_indent="    ", subsequent_indent="    ")
        headers_of_interest = set(_normalize_header(line) for line in args.headers)

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


IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# primer ambiguity code -> each concrete base it matches, in edlib's additionalEqualities format
IUPAC_EQUALITIES = [(code, base) for code, bases in IUPAC_CODES.items()
                    if len(bases) > 1 for base in bases]


def find_invalid_primer_chars(primer):
    return sorted(set(primer.upper()) - set(IUPAC_CODES))


def extract_seqs_by_primers(args):
    in_fasta = args.input_fasta
    fwd = args.forward_primer.upper().strip()
    rev = args.reverse_primer.upper().strip()
    max_amplicon_length = getattr(args, "max_amplicon_length", None)

    hit_count = 0
    alt_orientation_counts = {alt: 0 for alt in ORIENTATION_ALTERNATIVES}

    with open(args.output_fasta, "w") as out_fasta:

        for header, seq, hits in _primer_hits_per_record(in_fasta, fwd, rev, args.max_mismatches):
            amplicons = pair_primer_hits(seq, hits, max_amplicon_length)

            # tallies what a primer-orientation mix-up would give; it reuses the hits already
            # found, so it's cheap to do for every sequence
            for alt, count in count_alt_orientation_amplicons(seq, hits, max_amplicon_length).items():
                alt_orientation_counts[alt] += count

            for _, (left_label, right_label, left_start, right_end, amplicon, length) in enumerate(amplicons):

                out_header = f"{header}|{left_label}-to-{right_label}|{left_start}-{right_end}|{length}"
                out_seq = amplicon

                out_fasta.write(f">{out_header}\n")
                out_fasta.write(f"{out_seq}\n")
                hit_count += 1

    report_alt = should_report_alt_orientation(alt_orientation_counts, hit_count)

    if hit_count == 0:
        os.remove(args.output_fasta)
        report_message("No sequences were found based on the provided primers.",
                       trailing_newline=not report_alt)
    else:
        report_message(f"Extracted {hit_count} sequence(s) based on the provided primers, written to:", color="none")
        report_message(args.output_fasta, color="yellow", initial_indent="    ", leading_newline=False,
                       trailing_newline=not report_alt)

    if report_alt:
        report_alt_orientation_notice(alt_orientation_counts, hit_count, fwd, rev)


def _primers(fwd, rev):
    return [
        ("fwd", fwd),
        ("rev", rev),
        ("fwd_rc", revcomp(fwd)),
        ("rev_rc", revcomp(rev)),
    ]


def _primer_hits_per_record(in_fasta, fwd, rev, max_mismatches):
    """
    (header, uppercased seq, primer hits) for each record of `in_fasta` (as read_fasta() reads
    it). The primer search runs in the compiled core when it can (the same masking, IUPAC-aware
    edlib search, and repeated rounds as _find_primer_sites(), across threads); records it
    reports without two or more sites can't hold an amplicon in any orientation, so it skips
    them. Otherwise (e.g., non-ASCII input), it's done here.
    """
    primers = _primers(fwd, rev)
    _check_primer_lengths(primers, max_mismatches)
    try:
        from bit._native import fasta_primer_sites
        threads = min(8, os.cpu_count() or 1)
        native = fasta_primer_sites(str(in_fasta), [p for _, p in primers], max_mismatches,
                                    IUPAC_EQUALITIES, threads)
    except (ImportError, ValueError, OverflowError, TypeError):
        native = None

    if native is not None:
        for header, seq, sites in native:
            yield header, seq, _hits_from_sites(seq, primers, sites)
        return

    for header, seq in read_fasta(in_fasta):
        seq = seq.upper()
        yield header, seq, find_all_primer_hits(seq, fwd, rev, max_mismatches)


def _check_primer_lengths(primers, max_mismatches):
    # a primer no longer than the allowed differences would match anywhere (and the
    # search-mask-search loop below couldn't finish)
    for label, primer in primers:
        if len(primer) <= max_mismatches:
            raise ValueError(f"primers need to be longer than the allowed differences ({max_mismatches})")


def _mask_target(seq):
    # IUPAC equalities in edlib are symmetric, so any ambiguity code in the target (e.g. N-gaps
    # in scaffolds) would match primer bases; masking non-ACGT positions with a character that
    # matches nothing keeps primer-side degeneracy from bleeding into target-side ambiguity
    return re.sub(r"[^ACGT]", "X", seq.upper())


def _find_primer_sites(primer, masked_seq, max_mismatches):
    # edlib's "locations" task only returns locations at the *best* edit distance it finds, so
    # an exact site would hide every 1- or 2-difference site elsewhere in the target. To get all
    # sites within max_mismatches, sites found in each round are masked and the search is
    # repeated until none are left within the limit
    locations = []
    target = masked_seq

    if not target:
        return locations

    while True:
        result = edlib.align(
            primer,
            target,
            mode = "HW",            # search within sequence
            task = "locations",     # return match positions
            k = max_mismatches,
            additionalEqualities = IUPAC_EQUALITIES,
        )

        if result["editDistance"] == -1:
            break

        round_locations = [(start, end_inclusive + 1) for start, end_inclusive in result["locations"]]
        locations.extend(round_locations)

        pieces = []
        prev_end = 0
        for start, end in sorted(round_locations):
            start = max(start, prev_end)
            pieces.append(target[prev_end:start])
            pieces.append("X" * max(0, end - start))
            prev_end = max(prev_end, end)
        pieces.append(target[prev_end:])
        target = "".join(pieces)

    return locations


def find_all_primer_hits(seq, fwd, rev, max_mismatches = 0):

    primers = _primers(fwd, rev)
    _check_primer_lengths(primers, max_mismatches)
    masked_seq = _mask_target(seq)
    sites = [_find_primer_sites(primer, masked_seq, max_mismatches) for _, primer in primers]

    return _hits_from_sites(seq, primers, sites)


def _hits_from_sites(seq, primers, sites):
    """The sorted, de-duplicated hits from each primer's (start, end) sites."""

    hits = []

    for (label, primer), primer_sites in zip(primers, sites):

        for start, end in primer_sites:

            hits.append(
                (
                    label,
                    start,
                    end,
                    seq[start:end]
                )
            )

    # sorting by start, then end, then label so ties resolve deterministically
    hits.sort(key = lambda x: (x[1], x[2], x[0]))

    # remove potential duplicates
    # label is part of the key so a forward and reverse primer hitting identical
    # coordinates (e.g., fully overlapping primers) are both retained
    seen = set()
    unique_hits = []

    for hit in hits:
        key = (hit[0], hit[1], hit[2])
        if key not in seen:
            seen.add(key)
            unique_hits.append(hit)

    return unique_hits


# the only primer-site arrangements that can amplify: on the plus strand the forward primer
# is found as given with the reverse primer's revcomp downstream; on the minus strand it's
# the reverse primer as given with the forward primer's revcomp downstream
VALID_PAIRINGS = (("fwd", "rev_rc"), ("rev", "fwd_rc"))


def find_amplicons(seq, fwd, rev, max_mismatches = 0, max_amplicon_length = None):
    hits = find_all_primer_hits(seq, fwd, rev, max_mismatches)
    return pair_primer_hits(seq, hits, max_amplicon_length)


def pair_primer_hits(seq, hits, max_amplicon_length = None):
    results = []
    seen_regions = set()

    for left_label_wanted, right_label_wanted in VALID_PAIRINGS:

        lefts = [h for h in hits if h[0] == left_label_wanted]
        rights = [h for h in hits if h[0] == right_label_wanted]

        # each left site pairs only with its nearest downstream partner. primers may overlap
        # (partially or fully), but the right primer can't start before or end before the left one
        nearest_right = {}
        for left in lefts:
            _, left_start, left_end, _ = left
            candidates = [r for r in rights if r[1] >= left_start and r[2] >= left_end]
            if candidates:
                nearest_right[left] = min(candidates, key = lambda r: (r[2], r[1]))

        # and when several left sites share the same nearest partner (e.g., overlapping or
        # nearby repeated primer sites), only the innermost one is kept, giving the shortest product
        innermost_left = {}
        for left, right in nearest_right.items():
            if right not in innermost_left or left[1] > innermost_left[right][1]:
                innermost_left[right] = left

        for right, left in innermost_left.items():
            left_label, left_start, left_end, left_primer = left
            right_label, right_start, right_end, right_primer = right

            # same region reached via different label combos (e.g., a palindromic primer
            # hitting as both fwd and fwd_rc) is reported once
            if (left_start, right_end) in seen_regions:
                continue
            seen_regions.add((left_start, right_end))

            amplicon = seq[left_start:right_end]
            length = len(amplicon)

            if max_amplicon_length is not None and length > max_amplicon_length:
                continue

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

    results.sort(key = lambda x: (x[2], x[3], x[0]))

    return results


# if nothing amplifies with the primers as given, these are the orientation mix-ups worth checking
# for. Each is just a relabeling of the hits already found (e.g., reverse complementing the reverse
# primer swaps which sites are "rev" and which are "rev_rc"), so no additional searching is needed
ORIENTATION_ALTERNATIVES = {
    "reverse": {"rev": "rev_rc", "rev_rc": "rev"},
    "forward": {"fwd": "fwd_rc", "fwd_rc": "fwd"},
    "both": {"fwd": "fwd_rc", "fwd_rc": "fwd", "rev": "rev_rc", "rev_rc": "rev"},
}


def count_alt_orientation_amplicons(seq, hits, max_amplicon_length = None):
    counts = {}
    for alt, relabel in ORIENTATION_ALTERNATIVES.items():
        relabeled = [(relabel.get(label, label), start, end, hit_seq) for label, start, end, hit_seq in hits]
        counts[alt] = len(pair_primer_hits(seq, relabeled, max_amplicon_length))
    return counts


def should_report_alt_orientation(alt_orientation_counts, hit_count):
    # correctly oriented primers can only produce flipped-orientation products by spanning between
    # copies of a target that sit on opposite strands, which gives fewer products than the real ones.
    # Misoriented primers can still produce some (those same cross-copy spans), so rather than only
    # checking for zero products, this flags whenever a flipped orientation would give more
    return max(alt_orientation_counts.values(), default = 0) > hit_count


def report_alt_orientation_notice(alt_orientation_counts, hit_count, fwd, rev):

    # most products wins; ties go to the most common mix-up (reverse primer, then forward, then both)
    alt = max(ORIENTATION_ALTERNATIVES, key = lambda a: alt_orientation_counts[a])
    count = alt_orientation_counts[alt]

    if alt == "reverse":
        which, verb = "the reverse primer", "was"
        suggestion = f"-r {revcomp(rev)}"
    elif alt == "forward":
        which, verb = "the forward primer", "was"
        suggestion = f"-f {revcomp(fwd)}"
    else:
        which, verb = "both primers", "were"
        suggestion = f"-f {revcomp(fwd)} -r {revcomp(rev)}"

    if hit_count == 0:
        as_given = "No sequences were extracted with the primers as given"
    else:
        as_given = (f"{hit_count} sequence(s) were extracted with the primers as given (misoriented "
                    f"primers can still span between target copies on opposite strands)")

    report_message(f"NOTICE: {as_given}, but {count} would be if {which} were reverse complemented. "
                   f"Both primers should be provided 5'->3' as they would be ordered. If {which} "
                   f"{verb} entered in the wrong orientation, try:",
                   leading_newline=True)
    report_message(suggestion, color="none", initial_indent="    ", leading_newline=False,
                   trailing_newline=True)
