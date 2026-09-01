import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_help, add_version_arg,
                            wrap_help)
from bit.modules.ncbi.dl_ncbi_assemblies import dl_ncbi_assemblies
from bit.modules.taxonomy.tax_ranks import RANKS
from bit.modules.taxonomy.get_accs_shared import ASSEMBLY_LEVELS
from bit.modules.taxonomy.target_taxon import SOURCE_CHOICES, SECTION_CHOICES


def build_parser(parent_subparsers=None, show_fine=False):
    """
    `show_fine` toggles whether the fine-tuning parameters carry real help strings
    (the `-H` / `--detailed-help` menu) or are hidden but still parsed, so they keep
    working and keep their defaults either way.
    """

    desc = """
        This program downloads assembly files for NCBI genomes. Targets can be pulled
        based on taxonomy (`-t`) from either GTDB or NCBI (also see the --derep-rank parameter),
        or they can be explicitly specified as assembly accessions (`-w`).
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "dl-ncbi-assemblies",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit dl-ncbi-assemblies -t Nitrospirota`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    # hides a fine-tuning arg from the standard menu while still registering it
    def h(text):
        return argparse.SUPPRESS if not show_fine else wrap_help(text)

    required = parser.add_argument_group("Required Parameters (one or both)")
    selection = parser.add_argument_group("Taxon-selection Parameters (used with `-t`)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        action="append",
        default=None,
        help=wrap_help("Target taxon to pull genomes for (a name, or 'all' for every "
                       "domain available in the source). Can be provided multiple times to "
                       "combine taxa; genomes shared between them are only downloaded "
                       "once. Can also be combined with `-w`."),
    )

    required.add_argument(
        "-w",
        "--wanted-accessions",
        metavar="<FILE>",
        help=wrap_help("Input file with wanted accessions, one per line"),
        default=None,
    )

    ## taxon-selection ##
    selection.add_argument(
        "--source",
        choices=list(SOURCE_CHOICES),
        default="gtdb",
        help=wrap_help("Which taxonomy to pull `-t` targets from (default: gtdb)"),
    )

    selection.add_argument(
        "--ncbi-section",
        choices=list(SECTION_CHOICES),
        default="both",
        help=wrap_help("Which part of NCBI to draw from (default: both). With the default "
                       "of `--derep-rank auto`, 'both' is typically fine. Ignored with "
                       "`--source gtdb`."),
    )

    selection.add_argument(
        "--derep-rank",
        choices=["auto", "off"] + list(RANKS),
        default="auto",
        help=wrap_help("Dereplicate down to a single genome per unique value of "
                       "this rank (default: auto). 'auto' is two ranks finer than the "
                       "target (one finer for eukaryotes). 'off' downloads every "
                       "genome under the taxon, so use with care :)"),
    )

    selection.add_argument(
        "-r", "--target-rank",
        choices=list(RANKS),
        default=None,
        help=wrap_help("Target rank (if needed to disambiguate a taxon name that "
                       "exists at multiple ranks)"),
    )

    selection.add_argument(
        "--target-domain",
        metavar="<STR>",
        default=None,
        help=wrap_help("Target domain (if needed to disambiguate a taxon name that "
                       "exists in multiple domains)"),
    )

    selection.add_argument(
        "--dry-run",
        action="store_true",
        help=wrap_help("Run the selection (based on `-t`), but just report how many genomes would be "
                       "would be downloaded"),
    )

    ## general ##
    optional.add_argument(
        "-f",
        "--format",
        help=wrap_help('Format to download (default: "fasta")'),
        default="fasta",
        choices=["fasta", "protein", "genbank", "gff", "nt_cds", "feature_tab", "report", "stats"],
    )

    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        help=wrap_help("Number of downloads you'd like to run concurrently. NCBI can "
                       "become unhappy with many requests, so a max of 20 will be "
                       "used even if more are requested (default: 10)"),
        default=10,
        type=int,
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help=wrap_help("Directory to save output files (default: current directory)"),
        default=".",
    )

    ## fine-tuning (detailed menu only) ##
    fine = parser.add_argument_group("Fine-tuning Parameters") if show_fine else parser

    fine.add_argument(
        "--assembly-level",
        action="append",
        choices=list(ASSEMBLY_LEVELS),
        default=None,
        help=h("Only include genomes (from `-t`) at these assembly levels. "
               "Can be provided multiple times."),
    )

    fine.add_argument(
        "--min-completeness",
        metavar="<FLOAT>",
        type=float,
        default=None,
        help=h("Don't include any genomes (from `-t`) below this checkm completeness "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
    )

    fine.add_argument(
        "--max-contamination",
        metavar="<FLOAT>",
        type=float,
        default=None,
        help=h("Don't include any genomes (from `-t`) above this checkm contamination "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
    )

    fine.add_argument(
        "--representatives-only",
        action="store_true",
        help=h("If `--source gtdb`, only pull GTDB representative genomes. "
               "If `--source ncbi`, only pull NCBI reference genomes. "
               "With the default of `--derep-rank auto`, these are often not needed."),
    )

    optional.add_argument(
        "-H",
        "--detailed-help",
        action="store_true",
        help=wrap_help("Show detailed help, including fine-tuning parameters"),
    )

    add_help(optional)
    add_version_arg(optional)

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

    dl_ncbi_assemblies(args)
