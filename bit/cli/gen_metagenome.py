import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_help, wrap_help,
                            add_version_arg, add_force, reconstruct_invocation)
from bit.modules.general import (notify_premature_exit, check_if_output_dir_exists,
                                 report_message)
from bit.modules.gen_metagenome import gen_metagenome


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        Build a mock metagenome with ground-truth tables. Genomes are selected
        from GTDB and/or supplied directly as accessions, downloaded,
        optionally mutated to set per-genome ANI, and then reads are generated at a
        chosen abundance distribution. Outputs include reads and per-genome, per-rank, and
        (optionally) per-read truth tables.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "gen-metagenome", description=desc,
            formatter_class=CustomRichHelpFormatter, add_help=False)
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit gen-metagenome -n 20`",
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
                       "(can be combined with --accessions to add specific genomes)")
    )

    required.add_argument(
        "-a",
        "--accessions",
        metavar="<FILE>",
        help=wrap_help("File of NCBI accessions to include (one per line), or a TSV "
                       "with an 'accession' column plus optional 'rel_abundance', "
                       "'coverage', and 'mutation_rate' columns to pin per-genome "
                       "values")
        )

    ## general ##
    general.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        default="gen-mg-outputs",
        help=wrap_help("Output directory (default: gen-mg-outputs)")
    )

    general.add_argument(
        "--output-prefix",
        metavar="<STR>",
        default="mg",
        help=wrap_help("Prefix for the output read files (default: mg)")
    )

    general.add_argument(
        "--per-read-tsv",
        action="store_true",
        default=False,
        help=wrap_help("Add this flag to write out a read-level truth table mapping every read to its "
                       "source genome, coordinates, and taxonomy (takes a lot more spacetime)")
    )

    add_force(general)

    general.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        type=int,
        default=10,
        help=wrap_help("Number of jobs to run in parallel where possible (capped at 10 for download step; default: 10)")
    )

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
        help="Show detailed help, including fine-tuning parameters"
    )

    add_version_arg(general)


    ## selection ##
    selection.add_argument(
        "-d",
        "--domain",
        default="both",
        choices=["bacteria", "archaea", "both"],
        help=wrap_help("Domain(s) to draw genomes from "
                       "(default: both), eukaryotes can be specified via --accessions")
    )

    selection.add_argument(
        "--derep-rank",
        default="species",
        choices=["domain", "phylum", "class", "order", "family", "genus",
                 "species", "off"],
        help=wrap_help("Keep at most one genome per unique taxon at this rank "
                       "(default: species), 'off' selects randomly with no "
                       "taxonomic dereplication")
    )


    ## abundance ##
    abundance.add_argument(
        "--abundance-mode",
        default="relative",
        choices=["relative", "coverage"],
        help=wrap_help("Whether the distribution describes relative abundance "
                       "(uses --total-reads) or per-genome fold-coverage (uses "
                       "--median-coverage) (default: relative)")
    )

    abundance.add_argument(
        "--abundance-dist",
        default="lognormal",
        choices=["lognormal", "even"],
        help=wrap_help("Shape of the abundance/coverage distribution "
                       "(default: lognormal)")
    )

    abundance.add_argument(
        "-R",
        "--total-reads",
        metavar="<INT>",
        type=int,
        default=None, # 10,000,000 set later
        help=wrap_help("Number of total reads to generate in 'relative' abundance mode (default: 10,000,000)")
    )

    abundance.add_argument(
        "-c",
        "--median-coverage",
        metavar="<FLOAT>",
        type=float,
        default=None,  # 30 set later
        help=h(wrap_help("Median per-genome coverage in 'coverage' abundance mode; scales the "
                         "distribution so 'even' gives this coverage to every genome and "
                         "'lognormal' centers around it (default: 30). Ignored in 'relative' mode"))
        )

    abundance.add_argument(
        "--sigma",
        metavar="<FLOAT>",
        type=float,
        default=None,  # 1.0 set later
        help=h(wrap_help("Sigma (spread) for the lognormal distribution; higher means "
                         "a longer abundance tail (default: 1.0)"))
        )


    ## mutation ##
    mutation.add_argument(
        "--mutation-mode",
        default="off",
        choices=["off", "uniform", "distributed"],
        help=wrap_help("Mutate genomes before read generation: "
                       "'uniform' where all are done at --mutation-rate; "
                       "'distributed' with each drawn between --mutation-rate-min and --mutation-rate-max; "
                       "(default: 'off')")
        )

    mutation.add_argument(
        "--mutation-rate",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.01 set later
        help=h(wrap_help("Mutation rate for 'uniform' mode (default: 0.01)"))
    )

    mutation.add_argument(
        "--mutation-rate-min",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.001 set later
        help=h(wrap_help("Lower bound for 'distributed' mode (default: 0.001)"))
    )

    mutation.add_argument(
        "--mutation-rate-max",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.05 set later
        help=h(wrap_help("Upper bound for 'distributed' mode (default: 0.05)"))
    )

    mutation.add_argument(
        "--ti-tv-ratio",
        metavar="<FLOAT>",
        type=float,
        default=None, # 1.0 set later
        help=h(wrap_help("Transition/transversion ratio for substitutions (default: 1.0)"))
    )

    mutation.add_argument(
        "--indel-rate",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.0 set later
        help=h(wrap_help("Fraction of mutations that are indels (default: 0.0)"))
    )

    # ---- Reads ----
    reads.add_argument(
        "-t",
        "--type",
        default="paired-end",
        choices=["paired-end", "single-end", "long"],
        help=wrap_help("Type of reads to generate (default: paired-end)")
    )

    reads.add_argument(
        "-r",
        "--read-length",
        metavar="<INT>",
        type=int,
        default=None,
        help=h(wrap_help("Read length (default: 150 for short reads, 5000 for long)"))
    )

    reads.add_argument(
        "--fragment-size",
        metavar="<INT>",
        type=int,
        default=500,
        help=h(wrap_help("Paired-end fragment size (default: 500)"))
    )

    reads.add_argument(
        "--fragment-size-range",
        metavar="<INT>",
        type=int,
        default=10,
        help=h(wrap_help("Percent +/- around --fragment-size (default: 10)"))
    )

    reads.add_argument(
        "--long-read-length-range",
        metavar="<INT>",
        type=int,
        default=50,
        help=h(wrap_help("Percent +/- around --read-length for long reads (default: 50)"))
    )

    reads.add_argument(
        "--include-Ns", action="store_true", default=False,
        help=h(wrap_help("Allow reads that include Ns (default: skip N-containing regions)"))
    )

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

    preflight_checks(args)

    # these are all set to None above and ultimately set here so that any
    # incompatible user inputs can be detected and reported in preflight_checks
    args.total_reads = args.total_reads or 10_000_000
    args.median_coverage = args.median_coverage or 30
    args.sigma = args.sigma or 1.0
    args.mutation_rate = args.mutation_rate or 0.01
    args.mutation_rate_min = args.mutation_rate_min or 0.001
    args.mutation_rate_max = args.mutation_rate_max or 0.05
    args.ti_tv_ratio = args.ti_tv_ratio or 1.0
    args.indel_rate = args.indel_rate or 0.0

    if args.type == "long" and not args.read_length:
        args.read_length = 5000
    elif not args.read_length:
        args.read_length = 150

    args.full_cmd_executed = reconstruct_invocation(parser, args)

    gen_metagenome(args)


def preflight_checks(args):

    if not args.num_genomes and not args.accessions:
        report_message("You must specify `--num-genomes` and/or provide input `--accessions` to enjoy this ride.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.output_dir:
        check_if_output_dir_exists(args.output_dir, force_overwrite=args.force_overwrite)


    # ---- mode/param contradictions ----

    if args.total_reads and args.abundance_mode == "coverage":
        report_message("Parameter `--total-reads` is incompatible with `--abundance-mode coverage`.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.median_coverage and args.abundance_mode == "relative":
        report_message("Parameter `--median-coverage` is incompatible with `--abundance-mode relative`.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.sigma and args.abundance_dist == "even":
        report_message("Parameter `--sigma` is incompatible with `--abundance-dist even`.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.mutation_mode == "off" and (args.mutation_rate or args.mutation_rate_min
                                        or args.mutation_rate_max or args.ti_tv_ratio
                                        or args.indel_rate):
        report_message("You must specify a `--mutation-mode` if wanting any other mutation parameters.`",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()


    # ---- mode-specific mutation params (set the wrong knob for the chosen mode) ----

    if args.mutation_mode == "uniform" and (args.mutation_rate_min is not None
                                            or args.mutation_rate_max is not None):
        report_message("Parameters `--mutation-rate-min`/`--mutation-rate-max` apply to "
                       "`--mutation-mode distributed`, not `uniform` (which uses `--mutation-rate`).",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.mutation_mode == "distributed" and args.mutation_rate is not None:
        report_message("Parameter `--mutation-rate` applies to `--mutation-mode uniform`, not "
                       "`distributed` (which uses `--mutation-rate-min`/`--mutation-rate-max`).",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    # ---- range / validity checks ----
    # only validate values the user actually supplied (None == unset, gets a
    # valid default later in main()); a real 0 survives here to be checked.

    def _fail(msg):
        report_message(msg, initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    def _check_range(val, name, lo=None, hi=None, lo_inclusive=True, hi_inclusive=True):
        """ fatal if a supplied value falls outside [lo, hi] (bounds optional). """
        if val is None:
            return
        if lo is not None:
            if (val < lo) or (val == lo and not lo_inclusive):
                _fail(f"Parameter `{name}` must be "
                      f"{'>= ' if lo_inclusive else '> '}{lo} (got {val}).")
        if hi is not None:
            if (val > hi) or (val == hi and not hi_inclusive):
                _fail(f"Parameter `{name}` must be "
                      f"{'<= ' if hi_inclusive else '< '}{hi} (got {val}).")

    # counts must be positive
    _check_range(args.num_genomes, "--num-genomes", lo=0, lo_inclusive=False)
    _check_range(args.total_reads, "--total-reads", lo=0, lo_inclusive=False)
    _check_range(args.jobs, "--jobs", lo=0, lo_inclusive=False)

    # coverage / spread must be positive
    _check_range(args.median_coverage, "--median-coverage", lo=0, lo_inclusive=False)
    _check_range(args.sigma, "--sigma", lo=0, lo_inclusive=False)
    _check_range(args.ti_tv_ratio, "--ti-tv-ratio", lo=0, lo_inclusive=False)

    # rates are fractions in [0, 1]
    _check_range(args.mutation_rate, "--mutation-rate", lo=0, hi=1)
    _check_range(args.mutation_rate_min, "--mutation-rate-min", lo=0, hi=1)
    _check_range(args.mutation_rate_max, "--mutation-rate-max", lo=0, hi=1)
    _check_range(args.indel_rate, "--indel-rate", lo=0, hi=1)

    # distributed bounds must be ordered (min < max); only when both supplied
    if (args.mutation_rate_min is not None and args.mutation_rate_max is not None
            and args.mutation_rate_min > args.mutation_rate_max):
        _fail(f"`--mutation-rate-min` ({args.mutation_rate_min}) must be "
              f"<= `--mutation-rate-max` ({args.mutation_rate_max}).")

    # read geometry
    _check_range(args.read_length, "--read-length", lo=0, lo_inclusive=False)
    _check_range(args.fragment_size, "--fragment-size", lo=0, lo_inclusive=False)
    _check_range(args.fragment_size_range, "--fragment-size-range", lo=0, hi=100, hi_inclusive=False)
    _check_range(args.long_read_length_range, "--long-read-length-range", lo=0, hi=100, hi_inclusive=False)

    # paired-end fragments must be at least as long as a read
    if (args.type == "paired-end" and args.read_length is not None
            and args.fragment_size is not None and args.fragment_size < args.read_length):
        _fail(f"`--fragment-size` ({args.fragment_size}) must be >= `--read-length` "
              f"({args.read_length}) for paired-end reads.")
