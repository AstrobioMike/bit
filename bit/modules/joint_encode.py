import sys
from Bio import SeqIO  # type: ignore
from bit.modules.seqs import calc_variation_in_msa, get_gap_fracs_per_col
from bit.modules.general import spinner
from colorama import Fore, init  # type: ignore

init(autoreset=True)


def joint_encode(args):

    print(f"\n  {Fore.YELLOW}Generating variability info for alignments...", file=sys.stderr)

    with spinner("", ""):
        variability_3di, variability_aa, gap_fracs = get_variability_and_gaps(args.input_3di_alignment,
                                                                          args.input_aa_alignment)

    print(f"\n  {Fore.YELLOW}Trimmin' and swappin'...", file=sys.stderr)

    with spinner("", ""):

        (actions, num_cols, num_kept, num_swapped,
         num_trimmed, threedi_cols, aa_cols) = determine_action_per_column(variability_3di,
                                                                        variability_aa,
                                                                        gap_fracs, args)

        output_alignment_file = args.output_prefix + "-alignment.fasta"
        output_partitions_file = args.output_prefix + "-partitions.nex"

        write_out_alignment(args, actions, output_alignment_file)

        write_out_partitions_file(threedi_cols, aa_cols, output_partitions_file)

    report_results(num_cols, num_trimmed, num_kept,
                   num_swapped, output_alignment_file,
                   output_partitions_file)


def get_variability_and_gaps(input_3di_alignment, input_aa_alignment):

    variability_3di = calc_variation_in_msa(input_3di_alignment, type="3Di")
    variability_aa = calc_variation_in_msa(input_aa_alignment, type="Protein")
    gap_fracs = get_gap_fracs_per_col(input_3di_alignment)

    return variability_3di, variability_aa, gap_fracs



def determine_action_per_column(variability_3di, variability_aa, gap_fracs, args):

    """
    Trimming/swapping logic
    1. if gap_fracs >= gap_threshold -> trim
    2. if var_3di >= 3di-upper-bound:
           if var_aa < AA-upper-bound -> swap to AA
           else -> trim
    3. if var_3di <= 3di-lower-bound:
           if var_aa > var_3di AND var_aa < AA-upper-bound -> swap to AA
           else -> keep 3Di
    4. Otherwise (3di-lower-bound < var_3di < 3di-upper-bound) -> keep 3Di
    """

    var_3di = variability_3di["variation"].to_numpy()
    var_aa = variability_aa["variation"].to_numpy()
    num_cols = len(var_3di)

    # determining the action to be taken for every column
    actions = []
    for i in range(num_cols):
        v3 = var_3di[i]
        va = var_aa[i]
        gf = gap_fracs[i]

        if gf >= args.gap_threshold:
            actions.append("trim")
        elif v3 >= args.threedi_upper_bound:
            if va < args.AA_upper_bound:
                actions.append("swap")
            else:
                actions.append("trim")
        elif v3 <= args.threedi_lower_bound:
            if va > v3 and va < args.AA_upper_bound:
                actions.append("swap")
            else:
                actions.append("keep")
        else:
            actions.append("keep")

    num_kept = sum(1 for a in actions if a != "trim")
    num_swapped = actions.count("swap")
    num_trimmed = actions.count("trim")

    threedi_cols, aa_cols = build_output_column_index_lists(actions)

    return actions, num_cols, num_kept, num_swapped, num_trimmed, threedi_cols, aa_cols


def build_output_column_index_lists(actions):

    # building output-column index lists (1-based, after trimming) for each character type

    out_col = 0
    threedi_cols = []
    aa_cols = []
    for action in actions:
        if action == "trim":
            continue
        out_col += 1
        if action == "swap":
            aa_cols.append(out_col)
        else:
            threedi_cols.append(out_col)

    return threedi_cols, aa_cols


def write_out_alignment(args, actions, output_alignment_file):

    with open(output_alignment_file, "w") as out_fh:
        for rec_3di, rec_aa in zip(
            SeqIO.parse(args.input_3di_alignment, "fasta"),
            SeqIO.parse(args.input_aa_alignment, "fasta"),
        ):
            seq_3di = str(rec_3di.seq)
            seq_aa = str(rec_aa.seq)

            new_seq = []
            for col_idx, action in enumerate(actions):
                if action == "trim":
                    continue
                elif action == "swap":
                    new_seq.append(seq_aa[col_idx])
                else:
                    new_seq.append(seq_3di[col_idx])

            out_fh.write(f">{rec_3di.id}\n{''.join(new_seq)}\n")


def write_out_partitions_file(threedi_cols, aa_cols, output_partitions_file):

    with open(output_partitions_file, "w") as part_fh:
        threedi_model = "GTR20+FO+R6"
        AA_model = "LG+FO+R6"
        part_fh.write("#nexus\nbegin sets;\n")
        if threedi_cols:
            part_fh.write(f"  charset 3di [3DI] = {_cols_to_ranges(threedi_cols)};\n")
        if aa_cols:
            part_fh.write(f"  charset aa [AA] = {_cols_to_ranges(aa_cols)};\n\n")
        if threedi_cols and aa_cols:
            part_fh.write(f"  charpartition mike = {threedi_model}:3di, {AA_model}:aa;\n")
        elif threedi_cols and not aa_cols:
            part_fh.write(f"  charpartition mike = {threedi_model}:3di;\n")
        elif not threedi_cols and aa_cols:
            part_fh.write(f"  charpartition mike = {AA_model}:aa;\n")
        part_fh.write("end;\n")


def _cols_to_ranges(cols):

    if not cols:
        return ""
    ranges = []
    start = cols[0]
    end = cols[0]
    for c in cols[1:]:
        if c == end + 1:
            end = c
        else:
            ranges.append(str(start) if start == end else f"{start}-{end}")
            start = end = c
    ranges.append(str(start) if start == end else f"{start}-{end}")
    return " ".join(ranges)


def report_results(num_cols, num_trimmed, num_kept, num_swapped, output_alignment_file, output_partitions_file):

    print(f"\n  {Fore.GREEN}Done!\n", file=sys.stderr)
    print(f"    {Fore.YELLOW}Summary:")
    print(f"        {'Total initial columns:':<24} {num_cols:,}")
    print(f"        {'Trimmed:':<24} {num_trimmed:,} ({num_trimmed/num_cols:.1%})")
    print(f"        {'Retained:':<24} {num_kept:,} ({num_kept/num_cols:.1%})")
    print(f"        {'  3di-to-AA swapped:':<24} {num_swapped:,} ({num_swapped/num_kept:.1%} of retained)\n")
    print(f"    {Fore.YELLOW}Output files:")
    print(f"        Alignment:   {output_alignment_file}")
    print(f"        Partitions:  {output_partitions_file}\n")
    print(f"  Note: you may want to modify the 'charpartition' block in the nexus file prior to use.\n")
