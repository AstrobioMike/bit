import numpy as np
from Bio import SeqIO  # type: ignore
from bit.modules.seqs import calc_variation_in_msa, get_gap_fracs_per_col
from bit.modules.general import spinner
from colorama import Fore, init  # type: ignore

init(autoreset=True)


def joint_encode(args):

    """
    Trimming/swapping logic
    1. gap_fracs >= gap_threshold -> trim
    2. var_3di > percentile cutoff AND var_aa > cutoff -> trim (both too variable)
    3. var_3di > percentile cutoff (but AA is not) -> swap to AA
    4. var_3di <= swap_threshold AND var_aa >= percentile cutoff -> trim (low 3Di var but AA too variable)
    5. var_3di <= swap_threshold AND var_aa >= var_3di -> swap to AA
    6. Otherwise -> keep 3Di
    """

    print(f"\n  {Fore.YELLOW}Generating variability info for alignments...")
    with spinner("", ""):
        variability_3di = calc_variation_in_msa(args.input_3di_alignment, type="3Di")
        variability_aa = calc_variation_in_msa(args.input_aa_alignment, type="Protein")
        gap_fracs = get_gap_fracs_per_col(args.input_3di_alignment)

    print(f"\n  {Fore.YELLOW}Trimmin' and swappin'...")
    with spinner("", ""):

        var_3di = variability_3di["variation"].to_numpy()
        var_aa = variability_aa["variation"].to_numpy()
        num_cols = len(var_3di)

        trim_cutoff_3di = np.percentile(var_3di, args.trim_variability_upper_percentile)
        trim_cutoff_aa = np.percentile(var_aa, args.trim_variability_upper_percentile)

        # determine the action for every column once, up front
        # priority: gap-trim > variability-trim > swap > keep
        actions = []
        for i in range(num_cols):
            v3 = var_3di[i]
            va = var_aa[i]
            gf = gap_fracs[i]

            if gf >= args.gap_threshold:
                actions.append("trim")
            elif v3 >= trim_cutoff_3di and va >= trim_cutoff_aa:
                actions.append("trim")
            elif v3 >= trim_cutoff_3di:
                actions.append("swap")
            elif v3 <= args.swap_3di_variability_threshold and va >= trim_cutoff_aa:
                actions.append("trim")
            elif v3 <= args.swap_3di_variability_threshold and va > v3:
                actions.append("swap")
            else:
                actions.append("keep")

        num_kept = sum(1 for a in actions if a != "trim")
        num_swapped = actions.count("swap")
        num_trimmed = actions.count("trim")

        output_alignment_file = args.output_prefix + "-alignment.fasta"
        output_partitions_file = args.output_prefix + "-partitions.nex"

        # build output-column index lists (1-based, after trimming) for each character type
        out_col = 0
        threeddi_cols = []
        aa_cols = []
        for action in actions:
            if action == "trim":
                continue
            out_col += 1
            if action == "swap":
                aa_cols.append(out_col)
            else:
                threeddi_cols.append(out_col)

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

        with open(output_partitions_file, "w") as part_fh:
            part_fh.write("#nexus\nbegin sets;\n")
            if threeddi_cols:
                part_fh.write(f"  charset 3di [3DI] = {_cols_to_ranges(threeddi_cols)};\n")
            if aa_cols:
                part_fh.write(f"  charset aa [AA] = {_cols_to_ranges(aa_cols)};\n")
            part_fh.write("end;\n")


    print(f"\n  {Fore.GREEN}Done!\n")
    print(f"    Summary:")
    print(f"        - {num_kept:,} / {num_cols:,} columns retained ({num_trimmed:,} trimmed)")
    print(f"        - {num_swapped:,} columns swapped from 3di to AA\n")
    print(f"    Output alignment written to:       {output_alignment_file}")
    print(f"    Output partitions file written to: {output_partitions_file}")
    print("        (Note: you may want to modify the 'charpartition' block of the partitions file prior to use)\n")


def _cols_to_ranges(cols):
    """Convert a sorted list of 1-based column indices into a space-separated range string, e.g. '1-5 8 10-12'."""
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
