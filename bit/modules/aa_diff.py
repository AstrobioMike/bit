from Bio import SeqIO  # type: ignore
from Bio.Align import PairwiseAligner, substitution_matrices  # type: ignore


def run_aa_diff(input_query_fa, ref_faa, seq_type, output_prefix):

    ref_record, query_record = _load_sequences(ref_faa, input_query_fa)

    ref_seq = str(ref_record.seq).upper().rstrip("*")
    frameshifts = []
    translated_path = None
    nt_query_len = None
    cds_path = None

    if seq_type == "nt":
        query_seq, frameshifts, translated_path, nt_query_len, cds_path = _run_miniprot_and_translate(
            input_query_fa, ref_faa, ref_seq, output_prefix
        )
    else:
        query_seq = str(query_record.seq).upper().rstrip("*")

    alignment = _align(ref_seq, query_seq)
    ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref_seq, query_seq)
    positions, insertions = _parse_alignment(ref_gapped, qry_gapped)
    mutations = _collect_mutations(positions, insertions, frameshifts)

    prefix = output_prefix
    _write_tsv(positions, insertions, f"{prefix}-all-positions.tsv")
    _write_mutations(mutations, f"{prefix}-mutations.txt")
    summary_text = _report_summary(positions, insertions, mutations, frameshifts, prefix, translated_path,
                                     nt_query_len=nt_query_len, cds_path=cds_path)
    _write_summary(summary_text, f"{prefix}-summary.txt")
    _write_alignment(ref_gapped, qry_gapped, ref_record.id, query_record.id, f"{prefix}-alignment.txt")


def _load_sequences(ref_faa, input_query_fa):
    ref_record = next(SeqIO.parse(ref_faa, "fasta"))
    query_record = next(SeqIO.parse(input_query_fa, "fasta"))
    return ref_record, query_record


def _run_miniprot_and_translate(query_nt_path, ref_faa_path, ref_seq, prefix):
    """
    Run miniprot to align a nucleotide query against the reference protein.
    Parses the cs tag to reconstruct the translated query AA sequence and
    capture frameshift events.

    Returns (translated_aa_seq, frameshifts, translated_path) where frameshifts
    is a list of dicts with keys ref_pos (1-based) and type ('+1' or '+2').
    """
    import re
    import subprocess
    from bit.modules.general import report_message, notify_premature_exit

    result = subprocess.run(
        ["miniprot", "--cs", "-N", "1", query_nt_path, ref_faa_path],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        report_message(f"miniprot exited with an error:\n{result.stderr.strip()}")
        notify_premature_exit()

    # Take the first non-comment PAF line (best hit with -N 1)
    hit = None
    for line in result.stdout.splitlines():
        if not line.startswith("#") and line.strip():
            hit = line.strip().split("\t")
            break

    if hit is None:
        report_message(
            f"miniprot found no alignment for '{query_nt_path}' against '{ref_faa_path}'."
        )
        notify_premature_exit()

    # In miniprot PAF: protein is the query (cols 0-4), nucleotide is the target (cols 5-8).
    # hit[2] is the protein start position (0-based); hit[7] is the nucleotide start.
    ref_start = int(hit[2])   # 0-based start in reference protein

    # Pull cs tag from optional fields
    cs_string = None
    for field in hit[12:]:
        if field.startswith("cs:Z:"):
            cs_string = field[5:]
            break

    if cs_string is None:
        report_message(
            "miniprot output did not contain a cs tag for some reason."
        )
        notify_premature_exit()

    translated_aa, frameshifts = _parse_cs_tag(cs_string, ref_seq, ref_start)

    nt_query_len = int(hit[6])   # nucleotide sequence length (PAF target length)
    query_id = hit[0]

    translated_path = f"{prefix}-inferred-protein.faa"
    with open(translated_path, "w") as f:
        f.write(f">{query_id} translated-by-miniprot\n")
        f.write(translated_aa + "\n")

    # Extract the coding nucleotide sequence using PAF target coordinates
    target_name = hit[5]
    target_start = int(hit[7])
    target_end = int(hit[8])
    strand = hit[4]
    cds_path = None
    for record in SeqIO.parse(query_nt_path, "fasta"):
        if record.id == target_name:
            from Bio.Seq import Seq as _Seq  # type: ignore
            cds_raw = str(record.seq)[target_start:target_end]
            cds_seq = str(_Seq(cds_raw).reverse_complement()) if strand == "-" else cds_raw
            cds_path = f"{prefix}-inferred-CDS.fasta"
            with open(cds_path, "w") as f:
                f.write(f">{query_id} coding-seq-from-miniprot\n")
                f.write(cds_seq + "\n")
            break

    return translated_aa, frameshifts, translated_path, nt_query_len, cds_path


def _parse_cs_tag(cs_string, ref_seq, ref_start):
    """
    Parse a miniprot cs tag and reconstruct the translated query protein sequence.

    miniprot cs tag operations (uppercase = AA level, lowercase = nucleotide level):
      :N      N consecutive identical AA matches
      *XY     AA substitution: ref has X, query has Y  (uppercase)
      +SEQ    AA insertion in query                     (uppercase)
      -SEQ    AA deletion from query                    (uppercase)
      ~s_l_c  intron: strand s, nt length l, splice c  (no AA contribution)
      fN      +1 frameshift: 1 extra nt N skipped       (lowercase nucleotide)
      bNN     +2 frameshift: 2 extra nt NN skipped      (lowercase nucleotides)

    Returns (translated_seq, frameshifts) where frameshifts is a list of dicts
    with keys ref_pos (1-based) and type ('+1' or '+2').
    """
    import re

    CS_OP = re.compile(
        r":(\d+)"                    # :N  matches
        r"|\*([acgt]{3})([A-Z])"     # *codon+qry_aa: 3-nt ref codon + 1-letter query AA
        r"|\+([A-Z]+)"               # +SEQ AA insertion
        r"|-([A-Z]+)"                # -SEQ AA deletion (uppercase = amino acids)
        r"|-([acgt]+)"               # -nt(s) frameshift: lowercase nt(s) in nucleotide with no protein residue
        r"|~[+\-]\d+[a-z]+"          # ~intron  (no capture needed)
        r"|(f[acgt])"                # fN  +1 frameshift (alternative encoding)
        r"|(b[acgt]{2})"             # bNN +2 frameshift (alternative encoding)
    )

    translated = []
    frameshifts = []
    ref_pos = ref_start  # 0-based index into ref_seq

    for m in CS_OP.finditer(cs_string):

        match_n, subst_codon, subst_aa, ins, deletion, fs_del, fs_f, fs_b = m.groups()

        if match_n:
            n = int(match_n)
            translated.append(ref_seq[ref_pos:ref_pos + n].upper())
            ref_pos += n
        elif subst_codon:
            # subst_codon is the 3-nt codon from the nucleotide (lowercase) — minimap2 cs convention:
            # lowercase = target (nucleotide), uppercase = query (reference protein).
            # Translate the codon to get what the nucleotide query actually encodes.
            from Bio.Seq import Seq  # type: ignore
            qry_aa = str(Seq(subst_codon.upper()).translate())
            translated.append(qry_aa)
            ref_pos += 1
        elif ins:
            translated.append(ins.upper())         # no ref advance
        elif deletion:
            ref_pos += len(deletion)               # no query AA output
        elif fs_del:
            # 1 or 2 lowercase nucleotides: extra nt(s) in nucleotide = frameshift
            # ref_pos does NOT advance (no protein residue consumed)
            fs_type = "+1" if len(fs_del) % 3 == 1 else "+2"
            frameshifts.append({"ref_pos": ref_pos + 1, "type": fs_type})
        elif fs_f:
            frameshifts.append({"ref_pos": ref_pos + 1, "type": "+1"})
        elif fs_b:
            frameshifts.append({"ref_pos": ref_pos + 1, "type": "+2"})
        # intron: matched by regex but no capture group fires → nothing to do

    return "".join(translated), frameshifts


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


def _collect_mutations(positions, insertions, frameshifts=None):
    """Build mutation strings: substitutions as A210R, deletions as A210del, insertions as ins210:KL, frameshifts as fs210+1."""
    mutations = []

    for p in positions:
        if p["change_type"] == "substitution":
            mutations.append(f"{p['ref_aa']}{p['ref_pos']}{p['query_aa']}")
        elif p["change_type"] == "deletion":
            mutations.append(f"{p['ref_aa']}{p['ref_pos']}del")

    for ins in insertions:
        mutations.append(f"ins{ins['after_ref_pos']}:{ins['inserted_seq']}")

    if frameshifts:
        for fs in frameshifts:
            mutations.append(f"fs{fs['ref_pos']}{fs['type']}")

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


def _write_summary(stats_text, summary_path):
    with open(summary_path, "w") as f:
        f.write(stats_text + "\n")


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


def _report_summary(positions, insertions, mutations, frameshifts, prefix, translated_path=None, nt_query_len=None, cds_path=None):
    from bit.modules.general import color_text

    n_sub = sum(1 for p in positions if p["change_type"] == "substitution")
    n_ins = len(insertions)
    n_fs  = len(frameshifts)
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

    mut_lines = [
        f"  Total AA mutations:  {n_sub + n_del_runs + n_ins:,}",
        f"    Substitutions:     {n_sub:,}",
        f"    Deletions:         {n_del_runs:,} ({n_del_aa:,} AAs total)",
        f"    Insertions:        {n_ins:,} ({n_inserted_total:,} AAs total)",
    ]


    if nt_query_len is not None:
        stats_lines = [
            f"  Reference length:           {ref_len:,} AAs",
            "",
            f"  Query nt length:            {nt_query_len:,} NTs",
            f"    Inferred protein length:  {query_len:,} AAs",
            f"      Aligned to ref:         {n_aligned_to_ref:,} AAs",
            "",
        ] + mut_lines + [
            "",
            f"  Frameshifts (nt):    {n_fs:,}" + (" (ignored in inferred and aligned protein)" if n_fs > 0 else ""),
        ]
    else:
        stats_lines = [
            f"  Reference length:    {ref_len:,} AAs",
            "",
            f"  Query length:        {query_len:,} AAs",
            f"    Aligned to ref:    {n_aligned_to_ref:,} AAs",
            "",
        ] + mut_lines

    stats_text = "\n".join(stats_lines)

    print()
    print(color_text("                      SUMMARY                      ", "yellow"))
    print(color_text("---------------------------------------------------", "yellow"))
    print(stats_text)
    print()
    print(color_text("                      OUTPUTS                      ", "yellow"))
    print(color_text("---------------------------------------------------", "yellow"))
    if cds_path:
        print(f"  Inferred CDS:         {cds_path}")
    if translated_path:
        print(f"  Inferred protein:     {translated_path}")
        print()
    elif cds_path:
        print()
    print(f"  Summary file:         {prefix}-summary.txt")
    print(f"  All-positions table:  {prefix}-all-positions.tsv")
    print(f"  Mutations file:       {prefix}-mutations.txt")
    print(f"  Alignment file:       {prefix}-alignment.txt")
    print()

    return stats_text
