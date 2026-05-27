import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg

def build_parser(parent_subparsers=None):

    desc = """
        This program performs various operations on fasta files. See subcommand-specific
        help menus for more info.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "fasta",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    add_help(parser)

    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers

    ### shared args ###
    def add_common_required_arguments(group):
        group.add_argument(
            "input_fasta",
            metavar = "<FILE>",
            nargs = '?',
            default = sys.stdin,
            help = "Input fasta file or stdin if none provided"
        )

    #################################################
    ### subcommand cli for calculating gc content ###
    #################################################
    calc_gc_desc = """
        This subcommand takes a nucleotide fasta or multifasta and returns GC content either per
        sequence or per sliding window.
        """

    calc_gc_parser = subparsers.add_parser(
        "calc-gc",
        help="Calculate GC content (per seq or per sliding window)",
        description=calc_gc_desc,
        epilog="Ex. usage: `bit-fasta calc-gc input.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    calc_gc_required = calc_gc_parser.add_argument_group("Required Parameters")
    calc_gc_optional = calc_gc_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(calc_gc_required)

    calc_gc_optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tsv file (default: "gc-out.tsv"; when there is only one '
             'input sequence and no window size set, this is only written if specified)',
        default=None,
    )
    calc_gc_optional.add_argument(
        "-w",
        "--window-size",
        metavar="<INT>",
        help="Sliding window size; if set, GC is calculated per window across each sequence",
        default=None,
        type=int,
    )
    calc_gc_optional.add_argument(
        "-s",
        "--step-size",
        metavar="<INT>",
        help="Step size between windows; only used when -w is set (default: 10)",
        default=10,
        type=int,
    )

    add_help(calc_gc_optional)

    add_version_arg(calc_gc_optional)

    calc_gc_parser.set_defaults(func="calc_gc")

    #######################################################
    ### subcommand cli for calculating variation in msa ###
    #######################################################
    calc_var_desc = """
        This subcommand takes an alignment in fasta format as input and returns the Shannon uncertainty values for each column
        (using: https://scikit.bio/docs/dev/generated/skbio.alignment.TabularMSA.conservation.html). In output, a "variation" value of 0 would
        mean the same character in all sequences for that position (highest conservation); 1 would mean equal probability of any character
        (greatest variability). "Conservation" column is inverse. As written, any ambiguous bases or residues are converted to gap characters.
        """

    calc_var_parser = subparsers.add_parser(
        "calc-var-in-msa",
        help="Calculate variation in a multiple-sequence alignment",
        description=calc_var_desc,
        epilog="Ex. usage: `bit-fasta calc-var-in-msa alignment.fasta -o variation.tsv`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    calc_var_required = calc_var_parser.add_argument_group("Required Parameters")
    calc_var_optional = calc_var_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(calc_var_required)

    calc_var_optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tab-separated file (default: "variation.tsv")',
        action="store",
        default="variation.tsv"
    )

    calc_var_optional.add_argument(
        "-t",
        "--type",
        help='Molecule type (default: "Protein")',
        choices=["Protein", "DNA", "3Di"],
        action="store",
        default="Protein"
    )

    calc_var_optional.add_argument(
        "-g",
        "--gap-treatment",
        help='How to treat gaps (default: "include")',
        choices=["include", "nan", "ignore", "error"],
        action="store",
        default="include"
    )

    add_help(calc_var_optional)

    add_version_arg(calc_var_optional)

    calc_var_parser.set_defaults(func="calc_var_in_msa")

    ##################################################
    ### subcommand cli for counting bases and seqs ###
    ##################################################
    count_desc = """
        This subcommand takes a fasta as input and returns the total number of characters if the input
        holds a single sequence, or some summary stats if it is a multifasta. If you specify an output
        file, it also produces a tab-delimited file with two columns (header and number of characters for
        each sequence).
        """

    count_parser = subparsers.add_parser(
        "count",
        help="Count number of seqs and characters",
        description=count_desc,
        epilog="Ex. usage: `bit-fasta count input.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    count_required = count_parser.add_argument_group("Required Parameters")
    count_optional = count_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(count_required)

    count_optional.add_argument(
        "-o",
        "--output-tsv",
        metavar="<FILE>",
        help='Name of output tab-delimited file if you want all sequence lengths written out'
    )

    add_help(count_optional)

    add_version_arg(count_optional)

    count_parser.set_defaults(func="count")

    ##############################################################
    ### subcommand cli for extracting sequences by coordinates ###
    ##############################################################
    extract_by_coords_desc = """
        This subcommand takes a fasta file and tab-delimited (bed) file specifying which contigs
        and coordinates are wanted, and it returns a fasta of the chopped out sequences.
        """

    extract_by_coords_parser = subparsers.add_parser(
        "extract-by-coords",
        help="Extract sequences based on coordinates provided in a bed file",
        description=extract_by_coords_desc,
        epilog="Ex. usage: `bit-fasta extract-by-coords -i input.fasta -b targets.bed`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    extract_by_coords_required = extract_by_coords_parser.add_argument_group("Required Parameters")
    extract_by_coords_optional = extract_by_coords_parser.add_argument_group("Optional Parameters")

    extract_by_coords_required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )
    extract_by_coords_required.add_argument(
        "-b",
        "--bed-file",
        help="Tab-delimited bed file of desired contigs and coordinates (3 columns - contig, start, end - no header, 0-based counting)",
        metavar="<FILE>"
    )

    extract_by_coords_optional.add_argument(
        "-o",
        "--output-fasta",
        help='Output fasta file (default: "extracted-seqs.fasta")',
        metavar="<FILE>",
        default="extracted-seqs.fasta"
    )

    add_help(extract_by_coords_optional)

    add_version_arg(extract_by_coords_optional)

    extract_by_coords_parser.set_defaults(func="extract_by_coords")


    ##############################################################
    ### subcommand cli for extracting sequences by header(s)  ###
    ##############################################################
    extract_by_headers_desc = """
        This subcommand takes a fasta file and specified headers and extracts the sequences
        with those headers (or does the inverse if wanted).
        """

    extract_by_headers_parser = subparsers.add_parser(
        "extract-by-headers",
        help="Extract sequences based on specified headers",
        description=extract_by_headers_desc,
        epilog="Ex. usage: `bit-fasta extract-by-headers -i input.fasta -H contig-1 contig-2`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    extract_by_headers_required = extract_by_headers_parser.add_argument_group("Required Parameters (choose one of `--headers` or `--file-with-headers`)")
    extract_by_headers_optional = extract_by_headers_parser.add_argument_group("Optional Parameters")

    extract_by_headers_required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )
    extract_by_headers_required.add_argument(
        "-H",
        "--headers",
        help="Headers of sequences (space-delimited if more than one)",
        metavar="<STR>",
        nargs="+"
    )
    extract_by_headers_required.add_argument(
        "-f",
        "--file-with-headers",
        help="File with headers of sequences (one header per line)",
        metavar="<FILE>"
    )

    extract_by_headers_optional.add_argument(
        "-o",
        "--output-fasta",
        help='Output fasta file (default: "extracted-seqs.fasta")',
        metavar="<FILE>",
        default="extracted-seqs.fasta"
    )
    extract_by_headers_optional.add_argument(
        "--inverse",
        help="If specified, we will extract all sequences [bold]other[/bold] than the provided headers (default: False)",
        action="store_true"
    )

    add_help(extract_by_headers_optional)

    add_version_arg(extract_by_headers_optional)

    extract_by_headers_parser.set_defaults(func="extract_by_headers")


    ###########################################################
    ### subcommand cli for extracting sequences by primers  ###
    ###########################################################
    extract_by_primers_desc = """
        This subcommand takes a fasta file and forward and reverse primer sequences, and it
        returns a multifasta of the sequences including the specified primers. It currently doesn't
        allow for degenerate bases in the primers, but it does allow for up to 2 mismatches.
        """

    extract_by_primers_parser = subparsers.add_parser(
        "extract-by-primers",
        help="Extract sequences based on forward and reverse primer sequences",
        description=extract_by_primers_desc,
        epilog="Ex. usage: `bit-fasta extract-by-primers -i input.fasta -f ForwardPrimerSeq -r ReversePrimerSeq`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    extract_by_primers_required = extract_by_primers_parser.add_argument_group("Required Parameters")
    extract_by_primers_optional = extract_by_primers_parser.add_argument_group("Optional Parameters")

    extract_by_primers_required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )
    extract_by_primers_required.add_argument(
        "-f",
        "--forward-primer",
        help="Forward primer sequence",
        metavar="<STR>"
    )
    extract_by_primers_required.add_argument(
        "-r",
        "--reverse-primer",
        help="Reverse primer sequence",
        metavar="<STR>"
    )

    extract_by_primers_optional.add_argument(
        "-o",
        "--output-fasta",
        help='Output fasta file (default: "extracted-seqs.fasta")',
        metavar="<FILE>",
        default="extracted-seqs.fasta"
    )
    extract_by_primers_optional.add_argument(
        "-m",
        "--max-mismatches",
        help="Maximum number of mismatches allowed between the primer and the target sequence (default: 0)",
        type=int,
        choices=[0, 1, 2],
        default=0
    )

    add_help(extract_by_primers_optional)

    add_version_arg(extract_by_primers_optional)

    extract_by_primers_parser.set_defaults(func="extract_by_primers")


    ##############################################
    ### subcommand cli for filtering by length ###
    ##############################################
    filter_by_length_desc = """
        This subcommand filters input fasta sequences based on length.
        """

    filter_by_length_parser = subparsers.add_parser(
        "filter-by-length",
        help="Filter sequences based on minimum/maximum length",
        description=filter_by_length_desc,
        epilog="Ex. usage: `bit-fasta filter-by-length -i input.fasta -m 1000 -o filtered.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    filter_by_length_required = filter_by_length_parser.add_argument_group("Required Parameters")
    filter_by_length_optional = filter_by_length_parser.add_argument_group("Optional Parameters")

    filter_by_length_required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )

    filter_by_length_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Name of output fasta file (default: "filtered.fasta")',
        default="filtered.fasta"
    )
    filter_by_length_optional.add_argument(
        "-m",
        "--min-length",
        metavar="<INT>",
        help="Minimum length retained",
        default=None
    )
    filter_by_length_optional.add_argument(
        "-M",
        "--max-length",
        metavar="<INT>",
        help="Maximum length retained",
        default="9223372036854775807",
    )

    add_help(filter_by_length_optional)

    add_version_arg(filter_by_length_optional)

    filter_by_length_parser.set_defaults(func="filter_by_length")


    ############################################
    ### subcommand cli for modifying headers ###
    ############################################
    modify_headers_desc = """
        This subcommand facilitates modifying or renaming headers in a fasta. You can specify a string to be the same
        for all headers (and a number will be appened to keep them unique), or give a prefix or suffix to be added to existing headers.
        """

    modify_headers_parser = subparsers.add_parser(
        "modify-headers",
        help="Modify or rename fasta headers",
        description=modify_headers_desc,
        epilog="Ex. usage: `bit-fasta modify-headers -i input.fasta -w contig`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    modify_headers_required = modify_headers_parser.add_argument_group("Required Parameters")
    modify_headers_optional = modify_headers_parser.add_argument_group("Optional Parameters")

    modify_headers_required.add_argument(
        "-i",
        "--input-fasta",
        help="Input fasta file",
        metavar="<FILE>",
        required=True
    )

    modify_headers_optional.add_argument(
        "-o",
        "--output-fasta",
        metavar="<FILE>",
        help='Output fasta file (default: "modified-headers.fasta")',
        default="modified-headers.fasta"
    )

    modify_headers_optional.add_argument(
        "-w",
        "--wanted-text",
        metavar="<STR>",
        help='Base name to give seqs when renaming to "<wanted-text>_<n>" (default: "s")'
    )

    modify_headers_optional.add_argument(
        "-p",
        "--prefix",
        metavar="<STR>",
        help="Prepend this text to the original header (include separator if you want one; cannot be combined with --wanted-text)",
    )

    modify_headers_optional.add_argument(
        "-s",
        "--suffix",
        metavar="<STR>",
        help="Append this text to the original header (include separator if you want one; cannot be combined with "
        "--wanted-text; enter as -s='-suffix' if wanting a dash at the front)",
    )

    add_help(modify_headers_optional)

    add_version_arg(modify_headers_optional)

    modify_headers_parser.set_defaults(func="modify_headers")


    #############################################
    ### subcommand cli for removing softwraps ###
    #############################################
    remove_wraps_desc = """
        This subcommand removes line wraps from a fasta file, joining wrapped sequence
        lines back into single lines per sequence.
        """

    remove_wraps_parser = subparsers.add_parser(
        "remove-wraps",
        help="Remove line wraps from a fasta file",
        description=remove_wraps_desc,
        epilog="Ex. usage: `bit-fasta remove-wraps input.fasta > unwrapped.fasta`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    remove_wraps_required = remove_wraps_parser.add_argument_group("Required Parameters")
    remove_wraps_optional = remove_wraps_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(remove_wraps_required)

    remove_wraps_optional.add_argument(
        "-o",
        "--output-fasta",
        metavar="<FILE>",
        help="Output fasta file (default: stdout)",
        default=None
    )

    add_help(remove_wraps_optional)

    add_version_arg(remove_wraps_optional)

    remove_wraps_parser.set_defaults(func="remove_wraps")


    ################################################
    ### subcommand cli for generating a bed file ###
    ################################################
    to_bed_desc = """
        This subcommand takes a fasta as input and returns a tab-delimited bed file.
        """

    to_bed_parser = subparsers.add_parser(
        "to-bed",
        help="Generate a bed file from a fasta",
        description=to_bed_desc,
        epilog="Ex. usage: `bit-fasta to-bed input.fasta -o output.bed`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    to_bed_required = to_bed_parser.add_argument_group("Required Parameters")
    to_bed_optional = to_bed_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(to_bed_required)

    to_bed_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Name of output bed file (default: "output.bed")',
        default="output.bed"
    )

    add_help(to_bed_optional)

    add_version_arg(to_bed_optional)

    to_bed_parser.set_defaults(func="to_bed")


    ####################################################
    ### subcommand cli for generating a genbank file ###
    ####################################################
    to_genbank_desc = """
        This subcommand takes a nucleotide fasta file and generates a file in minimal genbank format.
        """

    to_genbank_parser = subparsers.add_parser(
        "to-genbank",
        help="Generate a genbank file from a nucleotide fasta",
        description=to_genbank_desc,
        epilog="Ex. usage: `bit-fasta to-genbank input.fasta -o output.gb`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    to_genbank_required = to_genbank_parser.add_argument_group("Required Parameters")
    to_genbank_optional = to_genbank_parser.add_argument_group("Optional Parameters")
    add_common_required_arguments(to_genbank_required)

    to_genbank_optional.add_argument(
        "-o",
        "--output-file",
        metavar="<FILE>",
        help='Output genbank file (default: "output.gb")',
        default="output.gb"
    )

    add_help(to_genbank_optional)

    add_version_arg(to_genbank_optional)

    to_genbank_parser.set_defaults(func="to_genbank")

    return parser


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args when using a subcommand so approriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in ("-v", "--version"):
            from bit.modules.general import report_version
            report_version()
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            if sys.stdin.isatty():
                parser.subparsers.choices[cmd].print_help(sys.stderr)
                sys.exit(0)
            # else: stdin is being piped, fall through to parse_args()
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()

    from bit.modules.general import check_files_are_found

    if args.input_fasta is not sys.stdin:
        check_files_are_found([args.input_fasta])

    import io

    if args.input_fasta is sys.stdin:
        text = args.input_fasta.read()
        args.input_fasta = io.StringIO(text)

    func_map = {
        "calc_gc": run_calc_gc,
        "calc_var_in_msa": run_calc_var_in_msa,
        "count": run_count,
        "extract_by_coords": run_extract_by_coords,
        "extract_by_headers": run_extract_by_headers,
        "extract_by_primers": run_extract_by_primers,
        "filter_by_length": run_filter_by_length,
        "modify_headers": run_modify_headers,
        "remove_wraps": run_remove_wraps,
        "to_bed": run_to_bed,
        "to_genbank": run_to_genbank,
    }

    func = func_map[args.func]

    func(args)


def run_calc_gc(args):

    window_mode = args.window_size is not None

    if window_mode:
        from bit.modules.seqs import calc_gc_sliding_window

        results = calc_gc_sliding_window(
            input_fasta=args.input_fasta,
            window=args.window_size,
            step=args.step_size,
        )
        out_file = args.output_tsv or "gc-out.tsv"
        with open(out_file, "w") as out_tsv:
            out_tsv.write(f"header\tlength\tgc\tgc_per_window_size_{args.window_size}_with_step_size_{args.step_size}\n")
            for item in results:
                windows_string = ", ".join(f"{x:.2f}" for x in item["gc_of_windows"])
                out_tsv.write(f"{item['header']}\t{item['length']}\t{item['gc']}\t{windows_string}\n")
        print(f"\n  GC (window size: {args.window_size}, step size: {args.step_size}) written to: '{out_file}'\n")

    else:
        from bit.modules.seqs import calc_gc_per_seq

        results = calc_gc_per_seq(input_fasta=args.input_fasta)

        if len(results) == 1:
            item = results[0]
            print(f"\n    Header:   {item['header']}")
            print(f"    Length:   {item['length']:,}")
            print(f"    GC:       {item['gc']}\n")
            if args.output_tsv is not None:
                with open(args.output_tsv, "w") as out_tsv:
                    out_tsv.write("header\tlength\tgc\n")
                    out_tsv.write(f"{item['header']}\t{item['length']}\t{item['gc']}\n")
                print(f"  GC content written to: '{args.output_tsv}'\n")
        else:
            out_file = args.output_tsv or "gc-out.tsv"
            with open(out_file, "w") as out_tsv:
                out_tsv.write("header\tlength\tgc\n")
                for item in results:
                    out_tsv.write(f"{item['header']}\t{item['length']}\t{item['gc']}\n")
            print(f"\n  GC content written to: '{out_file}'\n")


def run_calc_var_in_msa(args):

    from bit.modules.seqs import calc_variation_in_msa

    df = calc_variation_in_msa(input_alignment=args.input_fasta,
                               type=args.type,
                               gap_treatment=args.gap_treatment)

    df.to_csv(args.output_tsv, sep="\t", index=False)


def run_count(args):

    from bit.modules.seqs import parse_fasta_lengths

    stats = parse_fasta_lengths(args.input_fasta)
    total_length = stats["stats"]["total_length"]

    # if only one sequence, just printing out full length
    if stats["stats"]["n_seqs"] == 1:
        print(f"\n    Total length:    {total_length}\n")
    else: # if multiple, printout summary stats
        count_print_summary(stats["stats"], args.output_tsv)
        print(f"\n    Total length:    {total_length}\n")

    if args.output_tsv:
        count_write_lengths_to_tsv(stats["lengths"], args.output_tsv)


def count_write_lengths_to_tsv(seq_lengths, output_tsv):
    with open(output_tsv, "w") as out_file:
        for seq_id, length in seq_lengths.items():
            out_file.write(f"{seq_id}\t{length}\n")

    print(f"  All seq lengths written to: '{output_tsv}'\n")


def count_print_summary(stats, output_tsv):
    print("\n    Number of seqs:  " + str(stats["n_seqs"]))
    print("    Min length:      " + str(stats["min"]))
    print("    Max length:      " + str(stats["max"]))
    print("    Mean length:     " + str(stats["mean"]))
    print("    Median length:   " + str(stats["median"]))


def run_extract_by_coords(args):

    from bit.modules.extract_seqs import extract_seqs_by_coords

    extract_seqs_by_coords(args)


def run_extract_by_headers(args):

    from bit.modules.general import report_message, notify_premature_exit

    if not args.headers and not args.file_with_headers:
        report_message("You must provide either -H/--headers or -f/--file-with-headers.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()
    if args.headers and args.file_with_headers:
        report_message("You have provided both -H/--headers and -f/--file-with-headers parameters, please only provide one.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    from bit.modules.extract_seqs import extract_seqs_by_headers

    extract_seqs_by_headers(args)


def run_extract_by_primers(args):

    from bit.modules.extract_seqs import extract_seqs_by_primers

    extract_seqs_by_primers(args)


def run_filter_by_length(args):

    in_fasta = args.input_fasta
    out_fasta = args.output_file
    min_length = int(args.min_length) if args.min_length else None
    max_length = int(args.max_length)

    from bit.modules.seqs import filter_fasta_by_length

    (num_initial_seqs, num_seqs_retained,
     num_initial_bases, num_bases_retained) = filter_fasta_by_length(in_fasta, out_fasta,
                                                                     min_length, max_length)

    perc_seqs_retained = round(float(num_seqs_retained) / float(num_initial_seqs) * 100, 2)
    perc_bases_retained = round(float(num_bases_retained) / float(num_initial_bases) * 100, 2)

    print("\n    Retained " + f"{num_seqs_retained:,}" + " sequence(s) of the initial " + f"{num_initial_seqs:,}" + " (" + str(perc_seqs_retained) + "%).\n")
    print("    Retained " + f"{num_bases_retained:,}" + " bases of the initial " + f"{num_initial_bases:,}" + " (" + str(perc_bases_retained) + "%).\n")

    if num_seqs_retained == 0:
        import os
        from bit.modules.general import report_message
        report_message("  NOTICE: No sequences were retained after filtering with the specified length parameters.", leading_newline=False, trailing_newline=True)
        os.remove(out_fasta)


def run_to_bed(args):

    from bit.modules.seqs import fasta_to_bed

    bed_records = fasta_to_bed(args.input_fasta)

    with open(args.output_file, "w") as out:
        for name, start, end in bed_records:
            out.write(f"{name}\t{start}\t{end}\n")


def run_to_genbank(args):

    from bit.modules.seqs import fasta_to_genbank

    sequences = fasta_to_genbank(args.input_fasta)

    from Bio import SeqIO # type: ignore

    with open(args.output_file, "w") as out:
        SeqIO.write(sequences, out, "genbank")


def run_remove_wraps(args):

    from bit.modules.seqs import remove_wraps

    if args.output_fasta:
        with open(args.output_fasta, "w") as out:
            remove_wraps(args.input_fasta, out)
    else:
        remove_wraps(args.input_fasta, sys.stdout)


def run_modify_headers(args):


    if args.wanted_text and (args.prefix or args.suffix):

        from bit.modules.general import report_message, notify_premature_exit

        report_message("The --wanted-text argument is incompatible with --prefix and --suffix.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    elif not args.prefix and not args.suffix and not args.wanted_text:
        args.wanted_text = "s"

    from bit.modules.seqs import modify_fasta_headers

    modify_fasta_headers(args)
