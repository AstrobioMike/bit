import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_help, wrap_help,
                            add_version_arg, add_force, reconstruct_invocation)
from bit.modules.general import (notify_premature_exit, check_if_output_dir_exists,
                                 report_message)
from bit.modules.gen_metagenome import gen_metagenome


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        Build a mock metagenome with ground-truth tables. Genomes are: selected
        from GTDB and/or supplied directly as accessions; downloaded;
        optionally mutated to specified per-genome values; and then reads are generated at
        chosen abundance or coverage distributions. Outputs include reads, a 
        ground-truth assembly fasta, and per-genome, per-rank, and (optionally) 
        per-read truth tables with GTDB and NCBI taxonomy info.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "gen-mg", description=desc,
            formatter_class=CustomRichHelpFormatter, add_help=False)
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit gen-mg -n 20`",
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
        help=wrap_help("File of NCBI accessions to include (one per line), or a tsv "
                       "with an 'accession' column and optional 'rel_abundance', "
                       "'coverage', and/or 'mutation_rate' columns to specify per-genome "
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
        help=wrap_help("Number of jobs to run in parallel where possible (capped at 20 for the download step; default: 10)")
    )

    general.add_argument(
        "-s",
        "--seed",
        metavar="<INT>",
        help=("Set the random seed if wanting control over the random number generator (default: None)"),
        type=int,
        default=None,
    )

    add_help(general)

    add_version_arg(general)


    ## selection ##
    selection.add_argument(
        "-d",
        "--domain",
        default="both",
        choices=["bacteria", "archaea", "both"],
        help=wrap_help("Domain(s) to randomly draw genomes from "
                       "(default: both), eukaryotes can be specified via --accessions")
    )

    selection.add_argument(
        "--derep-rank",
        default="species",
        choices=["domain", "phylum", "class", "order", "family", "genus",
                 "species", "off"],
        help=wrap_help("Keep at most one genome per unique taxon at this rank "
                       "(default: species); 'off' selects randomly with no "
                       "taxonomic dereplication")
    )


    ## abundance ##
    abundance.add_argument(
        "--abundance-mode",
        default=None,  # 'relative' set later (None lets the input TSV's columns
                       # auto-select the mode, and conflicts be detected, first)
        choices=["relative", "coverage"],
        help=wrap_help("Specifies if genome abundance is governed by relative abundance "
                       "(which uses --total-reads) or coverage (which uses "
                       "--median-coverage) (default: relative)")
    )

    abundance.add_argument(
        "--abundance-dist",
        default="lognormal",
        choices=["lognormal", "even"],
        help=wrap_help("Specifies the type of abundance/coverage distribution "
                       "(default: lognormal)")
    )

    abundance.add_argument(
        "-R",
        "--total-reads",
        metavar="<INT>",
        type=int,
        default=None, # 10,000,000 set later
        help=wrap_help("Number of total reads (NOT read-pairs) to generate in 'relative' "
                       "abundance mode (default: 10,000,000)")
    )

    abundance.add_argument(
        "-c",
        "--median-coverage",
        metavar="<FLOAT>",
        type=float,
        default=None,  # 30 set later
        help=wrap_help("Median per-genome coverage in 'coverage' abundance mode; "
                       "'even' gives this coverage to every genome and "
                       "'lognormal' centers around it (default: 30)")
        )

    abundance.add_argument(
        "--spread",
        metavar="<FLOAT>",
        type=float,
        default=None,  # 1.0 set later
        dest="sigma",
        help=(wrap_help("Spread (sigma) for the lognormal distribution; higher means "
                         "a longer abundance tail (default: 1.0)"))
        )


    ## mutation ##
    mutation.add_argument(
        "--mutation-mode",
        default=None,  # 'off' set later (None lets a mutation_rate column in the
                       # input TSV auto-enable mutation, and conflicts be detected)
        choices=["off", "uniform", "varied"],
        help=wrap_help("Mutate genomes before read generation: "
                       "'uniform' means all are done at the set --mutation-rate; "
                       "'varied' means each is drawn between the set --mutation-rate-min "
                       "and --mutation-rate-max; (default: off)")
        )

    mutation.add_argument(
        "--mutation-rate",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.01 set later
        help=(wrap_help("Mutation rate for 'uniform' mode (default: 0.01)"))
    )

    mutation.add_argument(
        "--mutation-rate-min",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.001 set later
        help=(wrap_help("Lower bound for 'varied' mutation mode (default: 0.001)"))
    )

    mutation.add_argument(
        "--mutation-rate-max",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.05 set later
        help=(wrap_help("Upper bound for 'varied' mutationmode (default: 0.05)"))
    )

    mutation.add_argument(
        "--ti-tv-ratio",
        metavar="<FLOAT>",
        type=float,
        default=None, # 1.0 set later
        help=(wrap_help("Transition/transversion ratio for substitutions (default: 1.0)"))
    )

    mutation.add_argument(
        "--indel-rate",
        metavar="<FLOAT>",
        type=float,
        default=None, # 0.0 set later
        help=(wrap_help("Fraction of mutations that are indels (default: 0.0)"))
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
        help=(wrap_help("Read length (default: 150 for short reads, 5000 for long)"))
    )

    reads.add_argument(
        "--fragment-size",
        metavar="<INT>",
        type=int,
        default=500,
        help=(wrap_help("Paired-end fragment size (default: 500)"))
    )

    reads.add_argument(
        "--fragment-size-range",
        metavar="<INT>",
        type=int,
        default=10,
        help=(wrap_help("Percent +/- around --fragment-size (default: 10)"))
    )

    reads.add_argument(
        "--long-read-length-range",
        metavar="<INT>",
        type=int,
        default=50,
        help=(wrap_help("Percent +/- around --read-length for long reads (default: 50)"))
    )

    reads.add_argument(
        "--include-Ns", action="store_true", default=False,
        help=(wrap_help("Allow reads that include Ns (default: skip N-containing regions)"))
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

    # tracking whether the user explicitly set --abundance-dist before any
    # auto-resolution flips it; used so a single-genome coverage run only flips
    # to 'even' when the user left the dist at its default
    args.abundance_dist_explicit = "--abundance-dist" in argv

    preflight_checks(args)

    # these are all set to None above and ultimately set here so that any
    # incompatible user inputs can be detected and reported in preflight_checks
    args.total_reads = args.total_reads or 10_000_000
    # remember whether the user set --median-coverage before applying the default,
    # so the orchestrator can tell pinned coverages where un-pinned genomes draw
    # around the *default* (worth a notice) from an explicit choice (no notice).
    args.median_coverage_explicit = args.median_coverage is not None
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

    # resolve the seed up front (generating one pseudo-randomly when unset) so the
    # exact seed driving this run can be recorded to the runlog and reproduced.
    args.seed = _resolve_seed(args.seed)

    args.full_cmd_executed = reconstruct_invocation(parser, args)

    gen_metagenome(args)


def _resolve_seed(input_seed):
    """
    return the user's seed, or a pseudo-random one derived from the clock when
    unset — so a concrete seed always exists to log and reproduce the run (mirrors
    set_seed in mutate-seqs / add-insertion)
    """
    if input_seed is not None:
        return int(input_seed)
    import datetime
    now = datetime.datetime.now()
    return now.hour * 10000 + now.minute * 100 + now.second


def _count_accessions(path):
    """
    number of accession rows in the input file (0 if none/unreadable); handles
    both a bare one-per-line list and a tsv with an 'accession' header (which is
    not counted)
    """
    if not path:
        return 0
    try:
        lines = [ln.strip() for ln in open(path) if ln.strip()]
    except OSError:
        return 0
    if lines and "accession" in lines[0].lower():
        return len(lines) - 1   # drop header row
    return len(lines)


def _inspect_accession_columns(path):
    """
    Peek at the accessions input to see which optional pin columns it carries and,
    for mutation_rate, whether the supplied values are all equal. Returns a dict:
      {'has_rel_abundance', 'has_coverage', 'has_mutation_rate' (bools),
       'mutation_uniform' (bool|None), 'mutation_value' (float|None),
       'rel_abundance_sum' (float|None)}.
    A bare one-per-line file (no 'accession' header) has no columns.
    """
    import pandas as pd # type: ignore
    info = {"has_rel_abundance": False, "has_coverage": False,
            "has_mutation_rate": False, "mutation_uniform": None,
            "mutation_value": None, "rel_abundance_sum": None}
    try:
        first = open(path).readline()
    except OSError:
        return info
    if "accession" not in first.lower():
        return info                      # bare list, no columns
    df = pd.read_csv(path, sep="\t", dtype=str)

    def _has_values(col):
        # a column counts as present only if it exists AND has at least one
        # non-blank value; an all-empty column is treated as absent.
        if col not in df.columns:
            return False
        return pd.to_numeric(df[col], errors="coerce").notna().any()

    info["has_rel_abundance"] = _has_values("rel_abundance")
    info["has_coverage"] = _has_values("coverage")
    info["has_mutation_rate"] = _has_values("mutation_rate")
    if info["has_rel_abundance"]:
        vals = pd.to_numeric(df["rel_abundance"], errors="coerce").dropna()
        info["rel_abundance_sum"] = float(vals.sum()) if len(vals) else 0.0
    if info["has_mutation_rate"]:
        vals = pd.to_numeric(df["mutation_rate"], errors="coerce").dropna()
        uniq = set(vals.round(12).tolist())
        info["mutation_uniform"] = (len(uniq) <= 1)
        info["mutation_value"] = float(next(iter(uniq))) if len(uniq) == 1 else None
    return info


def resolve_input_driven_modes(args):
    """
    Let the input accessions TSV's columns drive the abundance and mutation modes
    when the user didn't set them explicitly, and reject explicit settings that
    contradict the file. Runs before the mode-contradiction checks so those see
    the resolved modes. Mutates args in place.
    """
    def fail(msg):
        report_message(msg, initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    cols = _inspect_accession_columns(args.accessions) if args.accessions else \
        {"has_rel_abundance": False, "has_coverage": False,
         "has_mutation_rate": False, "mutation_uniform": None,
         "mutation_value": None, "rel_abundance_sum": None}

    # ---- abundance mode ----
    if cols["has_coverage"] and not cols["has_rel_abundance"]:
        # only a coverage column -> coverage mode (unless the user said otherwise)
        if args.abundance_mode == "relative":
            fail("Input `--accessions` has a `coverage` column, which conflicts "
                 "with `--abundance-mode relative`. Use coverage mode (or omit "
                 "`--abundance-mode`) to honor the column.")
        args.abundance_mode = "coverage"
    elif cols["has_rel_abundance"] and not cols["has_coverage"]:
        if args.abundance_mode == "coverage":
            fail("Input `--accessions` has a `rel_abundance` column, which "
                 "conflicts with `--abundance-mode coverage`. Use relative mode "
                 "(or omit `--abundance-mode`) to honor the column.")
        args.abundance_mode = "relative"

    # --median-coverage is only meaningful in coverage mode; supplying it (without
    # an explicit relative mode) implies coverage mode, same as a `coverage` column.
    if getattr(args, "median_coverage", None) is not None:
        if args.abundance_mode == "relative":
            fail("`--median-coverage` conflicts with `--abundance-mode relative`. "
                "Omit `--abundance-mode` (or set it to coverage) to use it.")
        args.abundance_mode = "coverage"

    # both columns or neither -> fall through to default below; the resolved mode
    # picks which column is honored (the other is ignored).
    if args.abundance_mode is None:
        args.abundance_mode = "relative"

    # ---- pinned rel-abundances vs --num-genomes ----
    # if the user pins relative abundances that already sum to >= 1 there is no
    # room left for randomly-selected genomes, so that combination is a conflict.
    if (args.abundance_mode == "relative" and args.num_genomes
            and cols["rel_abundance_sum"] is not None
            and cols["rel_abundance_sum"] >= 1.0):
        fail(f"Input `rel_abundance` values sum to {cols['rel_abundance_sum']:.3f} "
             f"(>= 1), leaving no room for the `--num-genomes "
             f"{args.num_genomes}` randomly-selected genomes. Lower the pinned "
             f"abundances or drop `--num-genomes`.")

    # ---- mutation mode ----
    if cols["has_mutation_rate"]:
        if args.mutation_mode == "off":
            fail("Input `--accessions` has a `mutation_rate` column, which "
                 "conflicts with `--mutation-mode off`. Omit `--mutation-mode` "
                 "(or set uniform/varied) to honor the column.")
        if args.mutation_mode is None:
            # auto-enable from the column: all-equal supplied rates -> uniform at
            # that rate (random genomes get the same); differing rates -> varied,
            # with random genomes drawn from --mutation-rate-min/-max (defaults
            # unless the user set them).
            if cols["mutation_uniform"]:
                args.mutation_mode = "uniform"
                if args.mutation_rate is None and cols["mutation_value"] is not None:
                    args.mutation_rate = cols["mutation_value"]
            else:
                args.mutation_mode = "varied"
    if args.mutation_mode is None:
        args.mutation_mode = "off"

    # ---- single-genome coverage runs -> even distribution ----
    # with exactly one genome in coverage mode, 'lognormal' --abundance-dist would apply 
    # a random multiplier to the single coverage the user asked for. Setting 'even' so
    # that genome gets exactly --median-coverage, unless the user explicitly chose
    # a dist. The auto_even_single flag lets preflight_checks give a clearer message if
    # --spread was also supplied (since the user never typed `even` themselves)
    args.auto_even_single = False
    total_genomes = (args.num_genomes or 0) + _count_accessions(args.accessions)
    if (total_genomes == 1 and args.abundance_mode == "coverage"
            and not getattr(args, "abundance_dist_explicit", False)
            and args.abundance_dist != "even"):
        args.abundance_dist = "even"
        args.auto_even_single = True


def preflight_checks(args):

    if not args.num_genomes and not args.accessions:
        report_message("You must specify `--num-genomes` and/or provide input `--accessions` to enjoy this ride.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.output_dir:
        check_if_output_dir_exists(args.output_dir, force_overwrite=args.force_overwrite)

    # let the input TSV's columns drive abundance/mutation modes (and reject
    # explicit settings that contradict the file) before the checks below, which
    # then see the resolved modes.
    resolve_input_driven_modes(args)

    # ---- mode/param contradictions ----

    if args.total_reads and args.abundance_mode == "coverage":
        report_message("Parameter `--total-reads` is incompatible with `--abundance-mode coverage`.",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.sigma and args.abundance_dist == "even":
        # a single-genome coverage run auto-flips the dist to 'even' (unless explicitly set), 
        # so the generic "spread incompatible with even" message would be
        # confusing here; giving a clearer one instead
        if getattr(args, "auto_even_single", False):
            report_message("Parameter `--spread` has no effect on a single-genome coverage run "
                           "(the lone genome is given exactly `--median-coverage`, so the "
                           "distribution is set to 'even'). Drop `--spread` to proceed.",
                           initial_indent="    ", subsequent_indent="    ")
        else:
            report_message("Parameter `--spread` is incompatible with `--abundance-dist even`.",
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
                       "`--mutation-mode varied`, not `uniform` (which uses `--mutation-rate`).",
                       initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    if args.mutation_mode == "varied" and args.mutation_rate is not None:
        report_message("Parameter `--mutation-rate` applies to `--mutation-mode uniform`, not "
                       "`varied` (which uses `--mutation-rate-min`/`--mutation-rate-max`).",
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

    # coverage / spread (sigma) must be positive
    _check_range(args.median_coverage, "--median-coverage", lo=0, lo_inclusive=False)
    _check_range(args.sigma, "--sigma", lo=0, lo_inclusive=False)
    _check_range(args.ti_tv_ratio, "--ti-tv-ratio", lo=0, lo_inclusive=False)

    # rates are fractions in [0, 1]
    _check_range(args.mutation_rate, "--mutation-rate", lo=0, hi=1)
    _check_range(args.mutation_rate_min, "--mutation-rate-min", lo=0, hi=1)
    _check_range(args.mutation_rate_max, "--mutation-rate-max", lo=0, hi=1)
    _check_range(args.indel_rate, "--indel-rate", lo=0, hi=1)

    # varied bounds must be ordered (min < max); only when both supplied
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
