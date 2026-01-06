import os
import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help)
import pandas as pd
import pyfastx


def build_parser():

    desc = """
        This script outputs general summary stats for assemblies provided in fasta
        format. If an output file is specified, it writes the results there as a tsv.
        Otherwise it prints the results to the screen. "Ambiguous characters" reports
        total counts of any letter that is not "A", "T", "C", or "G". For version
        info, run `bit-version`.
    """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-summarize-assembly assembly.fasta",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )
    required = parser.add_argument_group("REQUIRED PARAMETERS")
    optional = parser.add_argument_group("OPTIONAL PARAMETERS")

    required.add_argument(
        "input_assemblies",
        metavar="<FILE(s)>",
        nargs="+",
        help="Input assembly file(s)"
    )

    optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv file (if none provided, prints to screen)',
        default=False
    )
    optional.add_argument(
        "-t",
        "--transpose-output-tsv",
        help='Set this flag if we want to have the output table have genomes as rows rather than columns.',
        action="store_true"
    )

    add_help(optional)

    return parser


def main():
    parser = build_parser()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    summarize_assemblies(
        input_assemblies = args.input_assemblies,
        output_tsv = args.output_tsv,
        transpose_output_tsv = args.transpose_output_tsv
    )


def summarize_assemblies(input_assemblies, output_tsv=False, transpose_output_tsv=False):
    df, use_paths_instead_of_basenames = setup_master_df(input_assemblies)
    df = summarize(df, input_assemblies, use_paths_instead_of_basenames)

    if output_tsv:
        if transpose_output_tsv:
            df = df.T
            df.to_csv(output_tsv, sep="\t", index=False, na_rep = "NA")
        else:
            df.to_csv(output_tsv, sep="\t", header=False, na_rep = "NA")

    else:
        print("")
        print(df.to_string(header=False))
        print("")


def setup_master_df(input_assemblies):
    df_colnames = []

    for assembly in input_assemblies:
        assembly_base = os.path.basename(assembly)
        df_colnames.append(assembly_base.rsplit(".", 1)[0])

    # checking for a situation where inputs may have the same basename, due to being from different directories
    # if so, setting a flag and reporting them as the full input paths instead of just basenames
    use_paths_instead_of_basenames = False
    for assembly in df_colnames:
        num_occurences = 0
        for assembly_2 in df_colnames:
            if assembly_2 == assembly:
                num_occurences += 1
        if num_occurences > 1:
            use_paths_instead_of_basenames = True
    if use_paths_instead_of_basenames:
        df_colnames = []
        for assembly in input_assemblies:
            df_colnames.append(assembly.rsplit(".", 1)[0])

    df_index = ["Assembly", "Total contigs", "Total length", "Ambiguous characters",
                "GC content", "Maximum contig length", "Minimum contig length", "N50",
                "N75", "N90", "L50", "L75", "L90", "Num. contigs >= 100",
                "Num. contigs >= 500", "Num. contigs >= 1000", "Num. contigs >= 5000",
                "Num. contigs >= 10000", "Num. contigs >= 50000", "Num. contigs >= 100000"]

    return pd.DataFrame(columns=df_colnames, index=df_index), use_paths_instead_of_basenames


def summarize(df, input_assemblies, use_paths_instead_of_basenames):
    for assembly in input_assemblies:
        if use_paths_instead_of_basenames:
            assembly_name = assembly.rsplit(".", 1)[0]
        else:
            assembly_base = os.path.basename(assembly)
            assembly_name = assembly_base.rsplit(".", 1)[0]

        df.at["Assembly", str(assembly_name)] = assembly_name
        # putting in a catch if file is empty (which can happen if an assembly produced no contigs)
        # this will leave it in the table, but with NAs (written out as "NA")
        if os.stat(assembly).st_size == 0:
            continue

        fasta = pyfastx.Fasta(assembly)

        df.at["Total contigs", str(assembly_name)] = len(fasta)
        df.at["Total length", str(assembly_name)] = fasta.size

        num_ambiguous_chars = 0
        for key in fasta.composition:
            if key not in ["A","T","G","C"]:
                num_ambiguous_chars += fasta.composition[key]

        df.at["Ambiguous characters", str(assembly_name)] = num_ambiguous_chars
        df.at["GC content", str(assembly_name)] = round(fasta.gc_content, 2)
        df.at["Maximum contig length", str(assembly_name)] = len(fasta.longest)
        df.at["Minimum contig length", str(assembly_name)] = len(fasta.shortest)

        info_at_50 = fasta.nl(50)
        info_at_75 = fasta.nl(75)
        info_at_90 = fasta.nl(90)
        df.at["N50", str(assembly_name)] = info_at_50[0]
        df.at["N75", str(assembly_name)] = info_at_75[0]
        df.at["N90", str(assembly_name)] = info_at_90[0]
        df.at["L50", str(assembly_name)] = info_at_50[1]
        df.at["L75", str(assembly_name)] = info_at_75[1]
        df.at["L90", str(assembly_name)] = info_at_90[1]

        df.at["Num. contigs >= 100", str(assembly_name)] = fasta.count(100)
        df.at["Num. contigs >= 500", str(assembly_name)] = fasta.count(500)
        df.at["Num. contigs >= 1000", str(assembly_name)] = fasta.count(1000)
        df.at["Num. contigs >= 5000", str(assembly_name)] = fasta.count(5000)
        df.at["Num. contigs >= 10000", str(assembly_name)] = fasta.count(10000)
        df.at["Num. contigs >= 50000", str(assembly_name)] = fasta.count(50000)
        df.at["Num. contigs >= 100000", str(assembly_name)] = fasta.count(100000)

        # removing intermediate index file
        os.remove(assembly + ".fxi")

    return df
