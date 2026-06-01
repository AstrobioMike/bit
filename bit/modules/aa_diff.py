from Bio import SeqIO  # type: ignore
from Bio.Align import PairwiseAligner, substitution_matrices  # type: ignore


def run_aa_diff(input_query_fa, ref_faa, seq_type, output_prefix, min_perc_id=30, min_perc_ref_cov=25):

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
    _check_alignment_thresholds(positions, min_perc_id, min_perc_ref_cov)
    mutations = _collect_mutations(positions, insertions, frameshifts)

    prefix = output_prefix
    _write_tsv(positions, insertions, f"{prefix}-all-positions.tsv")
    _write_mutations(mutations, f"{prefix}-mutations.txt")
    n_stops = query_seq.count("*") if seq_type == "nt" else 0
    summary_text = _report_summary(positions, insertions, mutations, frameshifts, prefix, translated_path,
                                     nt_query_len=nt_query_len, cds_path=cds_path, n_stops=n_stops)
    _write_summary(summary_text, f"{prefix}-summary.txt")
    _write_alignment(ref_gapped, qry_gapped, ref_record.id, query_record.id, f"{prefix}-alignment.txt")


def _load_sequences(ref_faa, input_query_fa):
    ref_record = next(SeqIO.parse(ref_faa, "fasta"))
    query_record = next(SeqIO.parse(input_query_fa, "fasta"))
    return ref_record, query_record


def _run_miniprot_and_translate(query_nt_path, ref_faa_path, ref_seq, prefix):
    """
    this runs miniprot to align a nucleotide query against the reference protein

    it then parses the cs tag to reconstruct the translated query-AA sequence and
    capture frameshift events

    it returns (translated_aa_seq, frameshifts, translated_path) where 'frameshifts'
    is a list of dicts with keys ref_pos (1-based) and type ('+1' or '+2')
    """
    import subprocess
    from bit.modules.general import report_message, notify_premature_exit

    result = subprocess.run(
        ["miniprot", "-L", "30", "--outc", "0", query_nt_path, ref_faa_path],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        report_message(f"miniprot exited with an error:\n{result.stderr.strip()}")
        notify_premature_exit()

    # collect all PAF hits
    all_hits = []
    for line in result.stdout.splitlines():
        if not line.startswith("#") and line.strip():
            all_hits.append(line.strip().split("\t"))

    if not all_hits:
        report_message(
            f"miniprot found no alignment between '{query_nt_path}' and '{ref_faa_path}'."
        )
        notify_premature_exit()

    # parse each hit and compute stats
    parsed_hits = []
    for h in all_hits:
        parsed = _parse_hit_stats(h, ref_seq)
        if parsed is not None:
            parsed_hits.append(parsed)

    if not parsed_hits:
        report_message(
            "miniprot output did not contain a cs tag for any hit."
        )
        notify_premature_exit()

    # write alignment summary table only when multiple hits are found
    if len(parsed_hits) > 1:
        # enrich each hit with its CDS nucleotide sequence before writing
        from Bio.Seq import Seq as _Seq  # type: ignore
        _nt_seqs = {record.id: str(record.seq) for record in SeqIO.parse(query_nt_path, "fasta")}
        for ph in parsed_hits:
            cds_raw = _nt_seqs.get(ph["hit"][5], "")[ph["nt_start"]:ph["nt_end"]]
            ph["cds_seq"] = str(_Seq(cds_raw).reverse_complement()) if ph["strand"] == "-" else cds_raw
        hits_table_path = f"{prefix}-alignment-summaries.tsv"
        _write_miniprot_hits_table(parsed_hits, hits_table_path)

    # use the top-ranked hit for the analysis
    top = parsed_hits[0]
    hit = top["hit"]
    translated_aa = top["translated_aa"]
    frameshifts = top["frameshifts"]

    # query nt seq, used for clip-extension and CDS extraction below
    target_name = hit[5]
    query_nt_seq = next(
        (str(r.seq) for r in SeqIO.parse(query_nt_path, "fasta") if r.id == target_name),
        "",
    )

    # recover residues miniprot soft-clipped past a terminal indel
    translated_aa, cds_seq = _extend_over_clips(translated_aa, hit, len(ref_seq), query_nt_seq)

    nt_query_len = int(hit[6])
    query_id = hit[0]

    translated_path = f"{prefix}-inferred-protein.faa"
    with open(translated_path, "w") as f:
        f.write(f">{query_id} translated-by-miniprot\n")
        f.write(translated_aa + "\n")

    cds_path = f"{prefix}-inferred-cds.fasta"
    with open(cds_path, "w") as f:
        f.write(f">{query_id} coding-seq-from-miniprot\n")
        f.write(cds_seq + "\n")

    return translated_aa, frameshifts, translated_path, nt_query_len, cds_path


def _parse_hit_stats(hit, ref_seq):
    """
    this extracts alignment statistics from a single miniprot PAF hit line

    in miniprot PAF: protein is the "query" (our ref; cols 0-4), nucleotide is the "target" (our query; cols 5-8):
      hit[0]  protein name (our ref)
      hit[2]  protein alignment start (0-based, AA units)
      hit[3]  protein alignment end (AA units)
      hit[4]  strand on the nucleotide
      hit[5]  nucleotide sequence name (our query)
      hit[6]  nucleotide length
      hit[7]  nucleotide alignment start
      hit[8]  nucleotide alignment end
      hit[9]  number of matching nucleotide bases (NOT AA matches — 3x AA count)
      hit[11] mapping quality
      AS:i:   alignment score

    returns None if the hit has no cs tag
    """
    cs_string = None
    score = None
    for field in hit[12:]:
        if field.startswith("cs:Z:"):
            cs_string = field[5:]
        elif field.startswith("AS:i:"):
            score = int(field[5:])

    if cs_string is None:
        return None

    ref_start = int(hit[2])
    translated_aa, frameshifts = _parse_cs_tag(cs_string, ref_seq, ref_start)

    ref_aa_span = int(hit[3]) - int(hit[2])

    # run the same alignment pipeline the main flow uses so perc_id/perc_ref_cov
    # are computed identically to what _check_alignment_thresholds sees
    alignment = _align(ref_seq, translated_aa)
    ref_gapped, qry_gapped = _get_gapped_seqs(alignment, ref_seq, translated_aa)
    positions, aln_insertions = _parse_alignment(ref_gapped, qry_gapped)
    stats = _calc_stats_from_positions(positions)

    n_sub = sum(1 for p in positions if p["change_type"] == "substitution")
    del_pos_list = [p["ref_pos"] for p in positions if p["change_type"] == "deletion"]
    n_del_aa = len(del_pos_list)
    n_del_runs = 0
    if del_pos_list:
        n_del_runs = 1
        for i in range(1, len(del_pos_list)):
            if del_pos_list[i] != del_pos_list[i - 1] + 1:
                n_del_runs += 1
    n_ins = len(aln_insertions)
    n_ins_aa = sum(len(ins["inserted_seq"]) for ins in aln_insertions)

    return {
        "hit": hit,
        "query_seq": hit[5],
        "nt_start": int(hit[7]),
        "nt_end": int(hit[8]),
        "strand": hit[4],
        "score": score if score is not None else stats["n_match"],
        "aln_aa_len": ref_aa_span,
        "perc_id": stats["perc_id"],
        "perc_ref_cov": stats["perc_ref_cov"],
        "n_frameshifts": len(frameshifts),
        "n_stops": translated_aa.count("*"),
        "n_total_mut": n_sub + n_del_runs + n_ins,
        "n_sub": n_sub,
        "n_del_runs": n_del_runs,
        "n_del_aa": n_del_aa,
        "n_ins": n_ins,
        "n_ins_aa": n_ins_aa,
        "translated_aa": translated_aa,
        "frameshifts": frameshifts,
    }


def _write_miniprot_hits_table(parsed_hits, tsv_path):
    """
    this writes a ranked summary table of all miniprot alignments to a TSV file
    and prints a brief notice to stdout
    """
    from bit.modules.general import color_text

    headers = ["score_rank", "query_seq", "nt_start", "nt_end", "strand", "score",
               "aln_aa_len", "perc_id", "perc_ref_cov", "total_mutations", "substitutions",
               "deletions", "insertions", "frameshifts", "stops",
               "inferred_cds", "inferred_protein"]

    with open(tsv_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for i, h in enumerate(parsed_hits, 1):
            row = [
                str(i),
                h["query_seq"],
                str(h["nt_start"]),
                str(h["nt_end"]),
                h["strand"],
                str(h["score"]),
                str(h["aln_aa_len"]),
                f"{h['perc_id']:.1f}",
                f"{h['perc_ref_cov']:.1f}",
                str(h["n_total_mut"]),
                str(h["n_sub"]),
                f"{h['n_del_runs']} ({h['n_del_aa']} AAs total)",
                f"{h['n_ins']} ({h['n_ins_aa']} AAs total)",
                str(h["n_frameshifts"]),
                str(h["n_stops"]),
                h.get("cds_seq", ""),
                h["translated_aa"],
            ]
            f.write("\t".join(row) + "\n")

    n = len(parsed_hits)
    print()
    print(color_text("----------------------------------- NOTICE -----------------------------------", "yellow"))
    print(color_text(f"  {n} miniprot alignments found", "yellow"))
    print(color_text("  The primary outputs and summary are for the highest-scoring alignment only", "yellow"))
    print(color_text(f"  All alignment summaries written to: {tsv_path}", "yellow"))
    print(color_text("------------------------------------------------------------------------------", "yellow"))
    print()


def _extend_over_clips(translated_aa, hit, ref_len, query_nt_seq):
    """
    miniprot soft-clips the alignment when extending past a terminal indel costs
    more than it gains, so terminal reference residues never appear in the cs tag
    and the clipped query residues are missing from the reconstruction

    here we translate the query nucleotides flanking the aligned region (in frame,
    up to the next stop codon or the end of the sequence) and prepend/append them
    to both the protein and the inferred CDS, so the two stay in sync and the
    downstream re-alignment resolves terminal indels instead of inventing them.

    returns (extended_aa, extended_cds), both in the protein's coding orientation.
    """

    from Bio.Seq import Seq  # type: ignore

    prot_start, prot_end = int(hit[2]), int(hit[3])
    strand = hit[4]
    nt_start, nt_end = int(hit[7]), int(hit[8])

    if strand == "-":
        coding = str(Seq(query_nt_seq).reverse_complement())
        L = len(query_nt_seq)
        a, b = L - nt_end, L - nt_start
    else:
        coding = query_nt_seq
        a, b = nt_start, nt_end

    n_ext = c_ext = n_ext_nt = c_ext_nt = ""

    if prot_end < ref_len:                     # C-terminal clip
        tail = coding[b:]
        tail = tail[: len(tail) // 3 * 3]       # whole codons only
        c_ext = str(Seq(tail).translate()).split("*")[0]
        c_ext_nt = tail[: 3 * len(c_ext)]

    if prot_start > 0:                         # N-terminal clip
        head = coding[a % 3 : a]                # keep frame so codons end at `a`
        head = head[: len(head) // 3 * 3]
        n_ext = str(Seq(head).translate()).split("*")[-1]
        n_ext_nt = head[len(head) - 3 * len(n_ext):]

    extended_aa  = n_ext + translated_aa + c_ext
    extended_cds = n_ext_nt + coding[a:b] + c_ext_nt

    # include the terminating stop codon if the query has one immediately after
    # the CDS (the protein deliberately stops before it). this also adds the stop
    # in the no-clip case, where the CDS otherwise ends at the last aligned codon.
    cds_end = b + len(c_ext_nt)
    next_codon = coding[cds_end : cds_end + 3]
    if len(next_codon) == 3 and str(Seq(next_codon).translate()) == "*":
        extended_cds += next_codon

    return extended_aa, extended_cds


def _parse_cs_tag(cs_string, ref_seq, ref_start):
    """
    this parses miniprot's cs tag and reconstructs the translated query-protein sequence

    miniprot cs tag operations (uppercase = AA level, lowercase = nucleotide level):
      :N        N consecutive identical AA matches
      *cccX     AA substitution: ccc is the 3-nt codon from the input nt query (lowercase),
                X is the ref protein AA (uppercase)
      +SEQ      ref residues absent from the query = DELETION (uppercase amino acids)
      -SEQ      uppercase form: not emitted by miniprot in practice
      -nts      query nucleotides absent from the ref = INSERTION (lowercase);
                a 1-2 nt remainder with no protein residue is a frameshift
      ~s_l_c    intron: strand s, nt length l, splice c     (no AA contribution)
      fN        +1 frameshift: 1 extra nt N skipped
      bNN       +2 frameshift: 2 extra nts NN skipped

    it returns (translated_seq, frameshifts) where 'frameshifts' is a list of dicts
    with keys ref_pos (1-based) and type ('+1' or '+2')
    """

    import re
    from Bio.Seq import Seq  # type: ignore

    CS_OP = re.compile(
        r":(\d+)"                    # :N  matches
        r"|\*([acgt]{3})([A-Z])"     # *cccX: 3-nt nt-query codon (lowercase) + ref protein AA (uppercase)
        r"|\+([A-Z]+)"               # +SEQ: uppercase ref AA(s) absent from query = DELETION
        r"|-([A-Z]+)"                # -SEQ: uppercase form, not emitted by miniprot in practice
        r"|-([acgt]+)"               # -nt(s): lowercase query nt absent from ref = INSERTION (+ frameshift remainder)
        r"|~[+\-]\d+[a-z]+"          # ~intron  (no capture needed)
        r"|(f[acgt])"                # fN  +1 frameshift (alternative encoding)
        r"|(b[acgt]{2})"             # bNN +2 frameshift (alternative encoding)
    )

    translated = []
    frameshifts = []
    ref_pos = ref_start

    for m in CS_OP.finditer(cs_string):

        match_n, subst_codon, subst_aa, plus_del, minus_seq, fs_del, fs_f, fs_b = m.groups()

        if match_n:
            n = int(match_n)
            translated.append(ref_seq[ref_pos:ref_pos + n].upper())
            ref_pos += n
        elif subst_codon:
            # subst_codon is the 3-nt codon from the nucleotide (lowercase); in our use here:
            # lowercase = query nucleotide, uppercase = ref protein AA
            # we translate the codon to get what the nucleotide query encodes
            qry_aa = str(Seq(subst_codon.upper()).translate())
            translated.append(qry_aa)
            ref_pos += 1
        elif plus_del:
            # '+' + uppercase ref AA(s) = residues missing from the query = deletion
            ref_pos += len(plus_del)        # skip the deleted ref residues, emit nothing
        elif minus_seq:
            # '-' + uppercase: not emitted by miniprot in practice; advance to be safe
            ref_pos += len(minus_seq)

        elif fs_del:  # '-' + lowercase nt = query sequence absent from ref = insertion
            nt = fs_del
            n_codons = len(nt) // 3
            if n_codons:
                translated.append(str(Seq(nt[:n_codons * 3].upper()).translate()))
            rem = len(nt) % 3
            if rem:
                frameshifts.append({"ref_pos": ref_pos + 1, "type": f"+{rem}"})
            # no ref advance — an insertion doesn't consume the reference
        elif fs_f:
            frameshifts.append({"ref_pos": ref_pos + 1, "type": "+1"})
        elif fs_b:
            frameshifts.append({"ref_pos": ref_pos + 1, "type": "+2"})
        # nothing to do for introns here

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
    """
    this reconstructs gap-inserted strings from the PairwiseAligner alignment coordinates
    """

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

    # append any trailing portions of either sequence not covered by the last block
    if ref_pos < len(ref_seq):
        ref_parts.append(ref_seq[ref_pos:])
        qry_parts.append("-" * (len(ref_seq) - ref_pos))
    if qry_pos < len(query_seq):
        ref_parts.append("-" * (len(query_seq) - qry_pos))
        qry_parts.append(query_seq[qry_pos:])

    return "".join(ref_parts), "".join(qry_parts)


def _parse_alignment(ref_gapped, qry_gapped):
    """
    this walks gapped-alignment strings and returns
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


def _calc_stats_from_positions(positions):
    """
    this computes percent identity and percent reference coverage from parsed alignment positions
    """
    ref_len = len(positions)
    n_match = sum(1 for p in positions if p["change_type"] == "match")
    n_aligned_to_ref = sum(1 for p in positions if p["change_type"] != "deletion")
    return {
        "ref_len": ref_len,
        "n_match": n_match,
        "n_aligned_to_ref": n_aligned_to_ref,
        "perc_id": n_match / ref_len * 100 if ref_len > 0 else 0.0,
        "perc_ref_cov": n_aligned_to_ref / ref_len * 100 if ref_len > 0 else 0.0,
    }


def _check_alignment_thresholds(positions, min_perc_id, min_perc_ref_cov):
    """
    this checks that the alignment meets minimum percent identity and percent reference
    coverage thresholds, and exits with a helpful message if not
    """
    from bit.modules.general import report_message, notify_premature_exit

    if not positions:
        report_message("No alignment positions were found.")
        notify_premature_exit()

    stats = _calc_stats_from_positions(positions)
    pct_id = stats["perc_id"]
    pct_ref_cov = stats["perc_ref_cov"]

    if pct_id < min_perc_id or pct_ref_cov < min_perc_ref_cov:
        report_message("No found alignment passed the minimum thresholds:")
        report_message(f"Percent identity:       {pct_id:.1f}%  (min: {min_perc_id}%)", initial_indent="    ", color="none")
        report_message(f"Percent ref covered:    {pct_ref_cov:.1f}%  (min: {min_perc_ref_cov}%)", initial_indent="    ", leading_newline=False, color="none")
        report_message("The sequences may be unrelated or highly divergent. If wanted, thresholds can be adjusted with the '--min-perc-id' and '--min-perc-ref-cov' parameters.", color="none")
        notify_premature_exit()


def _collect_mutations(positions, insertions, frameshifts=None):
    """
    this builds the mutation strings:
        substitutions as A210R
        deletions as A210del
        insertions as ins210:KL
        frameshifts as fs210+1
    """

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
    """
    this writes a BLAST-style pairwise alignment with Mike's handy position markers
    """

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


def _report_summary(positions, insertions, mutations, frameshifts, prefix, translated_path=None, nt_query_len=None, cds_path=None, n_stops=0):
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

    stats = _calc_stats_from_positions(positions)
    perc_id = stats["perc_id"]
    perc_ref_cov = stats["perc_ref_cov"]

    mut_lines = [
        f"  Total AA mutations:  {n_sub + n_del_runs + n_ins:,}",
        f"    Substitutions:     {n_sub:,}",
        f"    Deletions:         {n_del_runs:,} ({n_del_aa:,} AAs total)",
        f"    Insertions:        {n_ins:,} ({n_inserted_total:,} AAs total)",
    ]
    if n_stops > 0:
        mut_lines.append(f"    Stop codons:       {n_stops:,}")


    fs_line = f"  Frameshifts (nt):    {n_fs:,}" + (" (ignored in inferred and aligned protein)" if n_fs > 0 else "")

    if nt_query_len is not None:
        stats_lines = [
            f"  Reference length:           {ref_len:,} AAs",
            "",
            f"  Query nt length:            {nt_query_len:,} NTs",
            f"    Inferred protein length:  {query_len:,} AAs",
            f"      Aligned to ref:         {n_aligned_to_ref:,} AAs",
            f"        Percent identity:     {perc_id:.1f}%",
            f"        Percent ref covered:  {perc_ref_cov:.1f}%",
            "",
        ] + mut_lines + [
            "",
            fs_line,
        ]
        print_lines = stats_lines[:-1] + [color_text(fs_line, "orange") if n_fs > 0 else fs_line]
    else:
        stats_lines = [
            f"  Reference length:         {ref_len:,} AAs",
            "",
            f"  Query length:             {query_len:,} AAs",
            f"    Aligned to ref:         {n_aligned_to_ref:,} AAs",
            f"      Percent identity:     {perc_id:.1f}%",
            f"      Percent ref covered:  {perc_ref_cov:.1f}%",
            "",
        ] + mut_lines
        print_lines = stats_lines

    stats_text = "\n".join(stats_lines)
    print_text = "\n".join(print_lines)

    print()
    print(color_text("                      SUMMARY                      ", "yellow"))
    print(color_text("---------------------------------------------------", "yellow"))
    print(print_text)
    print("\n")
    print(color_text("                  PRIMARY OUTPUTS                  ", "yellow"))
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
