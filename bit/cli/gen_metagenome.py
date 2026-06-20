import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_help, wrap_help,
                            add_version_arg)
from bit.modules import gen_metagenome


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        Build a mock metagenome with ground-truth tables. Genomes are selected
        from the GTDB and/or supplied directly as accessions, downloaded, 
        optionally mutated to set per-genome ANI, then reads are generated at a 
        chosen abundance distribution. Outputs include reads plus per-genome, per-rank, and 
        (optionally) per-read truth tables.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "gen-metagenome", description=desc,
            formatter_class=CustomRichHelpFormatter, add_help=False)
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit gen-metagenome -n 20 --domains bacteria`",
            formatter_class=CustomRichHelpFormatter, add_help=False)

    # toggles whether detailed (fine-tuning) args show real help or are hidden
    def h(text):
        return argparse.SUPPRESS if not show_fine else text


    required = parser.add_argument_group("Required Parameters (one or both)")
    general = parser.add_argument_group("General Parameters")
    selection = parser.add_argument_group("Selection Parameters")
    abundance = parser.add_argument_group("Abundance Parameters")
    mutation = parser.add_argument_group("Mutation Parameters")
    reads = parser.add_argument_group("Read Parameters")


    ## required ##
    required.add_argument(
        "-n", 
        "--num-genomes", 
        metavar="<INT>", 
        type=int,
        help=wrap_help("Number of genomes to select from GTDB "
                       "(can be combined with --accessions to add specific genomes)"))
    required.add_argument(
        "-a",
        "--accessions", 
        metavar="<FILE>",
        help=wrap_help("File of NCBI accessions to include (one per line), or a TSV "
                       "with an 'accession' column plus optional 'rel_abundance', "
                       "'coverage', and 'mutation_rate' columns to pin per-genome "
                       "values"))

    ## general ##
    general.add_argument(
        "-o", 
        "--output-dir", 
        metavar="<DIR>", 
        default="gen-mg-outputs",
        help=wrap_help("Output directory (default: gen-mg-outputs)"))

    general.add_argument(
        "--output-prefix", 
        metavar="<STR>", 
        default="mg",
        help=wrap_help("Prefix for the output read files (default: mg)"))
    
    general.add_argument(
        "--per-read-tsv", 
        action="store_true", 
        default=False,
        help=wrap_help("Write out a read-level truth table mapping every read to its "
                       "source genome, taxonomy, and reference coordinates"))

    general.add_argument(
        "-j", 
        "--jobs", 
        metavar="<INT>", 
        type=int, 
        default=10,
        help=h(wrap_help("Concurrent downloads (capped at 20; default: 10)")))

    general.add_argument(
        "-s",
        "--seed",
        metavar="<INT>",
        help=h("Set the random seed if wanting control over the random number generator (default: None)"),
        type=int,
        default=None,
    )

    add_help(general)

    general.add_argument(
        "-H", 
        "--detailed-help", 
        action="store_true",
        help="Show detailed help, including fine-tuning parameters")

    add_version_arg(general)


    ## selection ##
    selection.add_argument(
        "-d",
        "--domains", 
        nargs="+", 
        default=["both"],
        choices=["bacteria", "archaea", "both"],
        help=wrap_help("Domain(s) to draw genomes from "
                       "(default: both). Eukaryotes can be specified via --accessions."))
    
    selection.add_argument(
        "--derep-rank", 
        default="species",
        choices=["domain", "phylum", "class", "order", "family", "genus",
                 "species", "off"],
        help=wrap_help("Keep at most one genome per unique taxon at this rank "
                       "(default: species). 'off' selects randomly with no "
                       "taxonomic dereplication."))


    ## abundance ##
    abundance.add_argument(
        "--abundance-mode", 
        default="relative", 
        choices=["relative", "coverage"],
        help=wrap_help("Whether the distribution describes relative abundance "
                       "(needs --total-reads) or per-genome fold-coverage directly "
                       "(default: relative)."))

    abundance.add_argument(
        "--abundance-dist", 
        default="lognormal",
        choices=["lognormal", "even"],
        help=wrap_help("Shape of the abundance/coverage distribution "
                       "(default: lognormal)."))

    abundance.add_argument(
        "--total-reads", 
        metavar="<INT>", 
        type=int, 
        default=10_000_000,
        help=wrap_help("Number of total reads to generate in 'relative' abundance mode (default: 10,000,000)"))

    abundance.add_argument(
        "--median-coverage", 
        metavar="<FLOAT>", 
        type=float, 
        default=30,
        help=h(wrap_help("Median per-genome coverage in 'coverage' abundance mode; scales the "
                         "distribution so 'even' gives this coverage to every genome and "
                         "'lognormal' centers around it (default: 30). Ignored in 'relative' mode.")))

    abundance.add_argument(
        "--sigma", 
        metavar="<FLOAT>", type=float, default=1.0,
        help=h(wrap_help("Sigma (spread) for the lognormal distribution; higher means "
                         "a longer abundance tail (default: 1.0)")))


    ## mutation ##
    mutation.add_argument(
        "--mutation-mode", 
        default="off",
        choices=["off", "uniform", "distributed"],
        help=wrap_help("Mutate genomes before read generation: "
                       "'uniform' where all are done at --mutation-rate; "
                       "'distributed' with each drawn between --mutation-rate-min and --mutation-rate-max; "
                       "(default: 'off')"))

    mutation.add_argument(
        "--mutation-rate", 
        metavar="<FLOAT>", 
        type=float, 
        default=0.01,
        help=h(wrap_help("Mutation rate for 'uniform' mode (default: 0.01)")))

    mutation.add_argument(
        "--mutation-rate-min", 
        metavar="<FLOAT>", 
        type=float, default=0.001,
        help=h(wrap_help("Lower bound for 'distributed' mode (default: 0.001)")))

    mutation.add_argument(
        "--mutation-rate-max", 
        metavar="<FLOAT>", 
        type=float, 
        default=0.05,
        help=h(wrap_help("Upper bound for 'distributed' mode (default: 0.05)")))

    mutation.add_argument(
        "--ti-tv-ratio", 
        metavar="<FLOAT>", 
        type=float, 
        default=1.0,
        help=h(wrap_help("Transition/transversion ratio for substitutions (default: 1.0)")))

    mutation.add_argument(
        "--indel-rate", 
        metavar="<FLOAT>", 
        type=float, default=0.0,
        help=h(wrap_help("Fraction of mutations that are indels (default: 0.0)")))

    # ---- Reads ----
    reads.add_argument(
        "-t", 
        "--type", 
        default="paired-end",
        choices=["paired-end", "single-end", "long"],
        help=wrap_help("Type of reads to generate (default: paired-end)"))
    
    reads.add_argument(
        "-r", 
        "--read-length", 
        metavar="<INT>", 
        type=int, 
        default=None,
        help=h(wrap_help("Read length (default: 150 for short reads, 5000 for long)")))
    
    reads.add_argument(
        "--fragment-size", 
        metavar="<INT>", 
        type=int, 
        default=500,
        help=h(wrap_help("Paired-end fragment size (default: 500)")))
    
    reads.add_argument(
        "--fragment-size-range", 
        metavar="<INT>", 
        type=int, 
        default=10,
        help=h(wrap_help("Percent +/- around --fragment-size (default: 10)")))

    reads.add_argument(
        "--long-read-length-range", 
        metavar="<INT>", 
        type=int, 
        default=50,
        help=h(wrap_help("Percent +/- around --read-length for long reads (default: 50)")))

    reads.add_argument(
        "--include-Ns", action="store_true", default=False,
        help=h(wrap_help("Allow reads that include Ns (default: skip N-containing regions).")))

    return parser


def main():

    argv = sys.argv[1:]
    if "-H" in argv or "--detailed-help" in argv:
        build_parser(show_fine=True).print_help(sys.stderr)
        sys.exit(0)

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.type == "long" and not args.read_length:
        args.read_length = 5000
    elif not args.read_length:
        args.read_length = 150

    if not args.num_genomes and not args.accessions:
        parser.error("provide at least one selection source: --num-genomes for "
                     "generative GTDB selection and/or --accessions for specific "
                     "genomes (one or both).")

    gen_metagenome(args)
