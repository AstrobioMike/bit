from Bio import SeqIO  # type: ignore
from Bio.Align import PairwiseAligner, substitution_matrices  # type: ignore


def run_aa_diff(args):
    """Align a query AA sequence to a reference AA sequence and report differences."""
    ref_record, query_record = _load_sequences(args)

    ref_seq = str(ref_record.seq).upper().rstrip("*")
    query_seq = str(query_record.seq).upper().rstrip("*")

    alignment = _align(ref_seq, query_seq)
    ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref_seq, query_seq)
    positions, insertions = _parse_alignment(ref_gapped, qry_gapped)
    mutations = _collect_mutations(positions, insertions)

    prefix = args.output_prefix
    _write_tsv(positions, insertions, f"{prefix}-all-positions.tsv")
    _write_mutations(mutations, f"{prefix}-mutations.txt")
    _write_alignment(ref_gapped, qry_gapped, ref_record.id, query_record.id, f"{prefix}-alignment.txt")
    _report_summary(positions, insertions, mutations, prefix)


def _load_sequences(args):
    ref_record = next(SeqIO.parse(args.ref_faa, "fasta"))
    query_record = next(SeqIO.parse(args.input_query_fa, "fasta"))
    return ref_record, query_record


def _align(ref_seq, query_seq):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    alignments = aligner.align(ref_seq, query_seq)
    return next(iter(alignments))


def _get_gapped_seqs(alignment, ref_seq, query_seq):
    """Reconstruct gap-inserted strings from PairwiseAligner alignment coordinates."""
    ref_blocks = alignment.aligned[0]
    qry_blocks = alignment.aligned[1]

    ref_parts = []
    qry_parts = []
    ref_pos = 0
    qry_pos = 0

    for (r_start, r_end), (q_start, q_end) in zip(ref_blocks, qry_blocks):
        if r_start > ref_pos:
            # deletion in query: ref has residues, query gets dashes
            ref_parts.append(ref_seq[ref_pos:r_start])
            qry_parts.append("-" * (r_start - ref_pos))
        if q_start > qry_pos:
            # insertion in query: query has residues, ref gets dashes
            ref_parts.append("-" * (q_start - qry_pos))
            qry_parts.append(query_seq[qry_pos:q_start])

        ref_parts.append(ref_seq[r_start:r_end])
        qry_parts.append(query_seq[q_start:q_end])
        ref_pos = r_end
        qry_pos = q_end

    return "".join(ref_parts), "".join(qry_parts)


def _parse_alignment(ref_gapped, qry_gapped):
    """
    Walk gapped alignment strings and return:
      positions - one dict per ref residue: ref_pos, ref_aa, query_aa, change_type
      insertions - one dict per inserted run: after_ref_pos, inserted_seq
    """
    positions = []
    insertions = []
    ref_pos = 0
    qry_pos = 0
    pending_ins = []

    for ref_char, qry_char in zip(ref_gapped, qry_gapped):
        if ref_char == "-":
            # query insertion relative to ref
            qry_pos += 1
            pending_ins.append(qry_char)
        else:
            if pending_ins:
                insertions.append({
                    "after_ref_pos": ref_pos,
                    "inserted_seq": "".join(pending_ins),
                })
                pending_ins = []

            ref_pos += 1

            if qry_char == "-":
                change = "deletion"
                this_qry_pos = None
            else:
                qry_pos += 1
                this_qry_pos = qry_pos
                change = "match" if ref_char == qry_char else "substitution"

            positions.append({
                "ref_pos": ref_pos,
                "ref_aa": ref_char,
                "query_pos": this_qry_pos,
                "query_aa": qry_char,
                "change_type": change,
            })

    if pending_ins:
        insertions.append({
            "after_ref_pos": ref_pos,
            "inserted_seq": "".join(pending_ins),
        })

    return positions, insertions


def _collect_mutations(positions, insertions):
    """Build mutation strings: substitutions as A210R, deletions as A210del, insertions as ins210:KL."""
    mutations = []

    for p in positions:
        if p["change_type"] == "substitution":
            mutations.append(f"{p['ref_aa']}{p['ref_pos']}{p['query_aa']}")
        elif p["change_type"] == "deletion":
            mutations.append(f"{p['ref_aa']}{p['ref_pos']}del")

    for ins in insertions:
        mutations.append(f"ins{ins['after_ref_pos']}:{ins['inserted_seq']}")

    return mutations


def _write_tsv(positions, insertions, tsv_path):
    ins_lookup = {ins["after_ref_pos"]: ins["inserted_seq"] for ins in insertions}

    with open(tsv_path, "w") as f:
        f.write("ref_pos\tref_aa\tquery_aa\tquery_pos\tchange_type\tinserted_before\n")
        for p in positions:
            qpos_str = str(p["query_pos"]) if p["query_pos"] is not None else "-"
            ins_before = ins_lookup.get(p["ref_pos"] - 1, "-")
            f.write(f"{p['ref_pos']}\t{p['ref_aa']}\t{p['query_aa']}\t{qpos_str}\t{p['change_type']}\t{ins_before}\n")


def _write_mutations(mutations, mut_path):
    with open(mut_path, "w") as f:
        if mutations:
            for m in mutations:
                f.write(m + "\n")
        else:
            f.write("# No mutations detected\n")


def _write_alignment(ref_gapped, qry_gapped, ref_id, query_id, aln_path, width=60):
    """Write a BLAST-style pairwise alignment with per-block position markers."""
    blosum62 = substitution_matrices.load("BLOSUM62")

    match_chars = []
    for r, q in zip(ref_gapped, qry_gapped):
        if r == "-" or q == "-":
            match_chars.append(" ")
        elif r == q:
            match_chars.append("|")
        else:
            try:
                score = blosum62[r, q]
            except KeyError:
                score = 0
            match_chars.append("+" if score > 0 else " ")
    match_line = "".join(match_chars)

    ref_label = f"Ref ({ref_id[:20]})"
    qry_label = f"Query ({query_id[:20]})"
    label_w = max(len(ref_label), len(qry_label))

    ref_pos = 0
    qry_pos = 0
    total_len = len(ref_gapped)

    with open(aln_path, "w") as f:
        for i in range(0, total_len, width):
            chunk_r = ref_gapped[i:i + width]
            chunk_q = qry_gapped[i:i + width]
            chunk_m = match_line[i:i + width]

            ref_start = ref_pos + 1
            qry_start = qry_pos + 1
            ref_end = ref_pos + sum(1 for c in chunk_r if c != "-")
            qry_end = qry_pos + sum(1 for c in chunk_q if c != "-")

            ruler_chars = []
            _rpos = ref_start - 1
            for _c in chunk_r:
                if _c == "-":
                    ruler_chars.append(" ")
                else:
                    _rpos += 1
                    ruler_chars.append("|" if _rpos % 10 == 0 else " ")
            ruler = "".join(ruler_chars)

            f.write(f"{'':>{label_w}}  {'':>6}  {ruler}\n")
            f.write(f"{ref_label:<{label_w}}  {ref_start:>6}  {chunk_r}  {ref_end}\n")
            f.write(f"{'':>{label_w}}  {'':>6}  {chunk_m}\n")
            f.write(f"{qry_label:<{label_w}}  {qry_start:>6}  {chunk_q}  {qry_end}\n")

            qruler_chars = []
            _qpos = qry_start - 1
            for _c in chunk_q:
                if _c == "-":
                    qruler_chars.append(" ")
                else:
                    _qpos += 1
                    qruler_chars.append("|" if _qpos % 10 == 0 else " ")
            f.write(f"{'':>{label_w}}  {'':>6}  {''.join(qruler_chars)}\n")

            f.write("\n")

            ref_pos = ref_end
            qry_pos = qry_end


def _report_summary(positions, insertions, mutations, prefix):
    from bit.modules.general import color_text

    n_sub = sum(1 for p in positions if p["change_type"] == "substitution")
    n_ins = len(insertions)
    ref_len = len(positions)

    del_positions = [p["ref_pos"] for p in positions if p["change_type"] == "deletion"]
    n_del_aa = len(del_positions)
    n_del_runs = 0
    if del_positions:
        n_del_runs = 1
        for i in range(1, len(del_positions)):
            if del_positions[i] != del_positions[i - 1] + 1:
                n_del_runs += 1

    n_aligned_to_ref = sum(1 for p in positions if p["change_type"] != "deletion")

    n_inserted_total = sum(len(ins["inserted_seq"]) for ins in insertions)
    query_len = n_aligned_to_ref + n_inserted_total

    print("")
    print(f"  Reference length:    {ref_len} AAs\n")
    print(f"  Query length:        {query_len} AAs")
    print(f"    Aligned to ref:    {n_aligned_to_ref} AAs\n")
    print(f"  Total mutations:     {n_sub + n_del_runs + n_ins}")
    print(f"    Substitutions:     {n_sub}")
    print(f"    Deletions:         {n_del_runs} ({n_del_aa} AAs total)")
    print(f"    Insertions:        {n_ins} ({n_inserted_total} AAs total)")
    print("")
    print(color_text(f"  All-positions table:  {prefix}-all-positions.tsv", "yellow"))
    print(color_text(f"  Mutations file:       {prefix}-mutations.txt", "yellow"))
    print(color_text(f"  Alignment file:       {prefix}-alignment.txt", "yellow"))
    print("")
