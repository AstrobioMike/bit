import sys
import argparse
from bit.cli.common import (CustomRichHelpFormatter, add_version_arg,
                            wrap_help)
from bit.modules.ncbi.dl_ncbi_assemblies import dl_ncbi_assemblies
from bit.modules.taxonomy.tax_ranks import RANKS
from bit.modules.taxonomy.get_accs_shared import ASSEMBLY_LEVELS
from bit.modules.taxonomy.target_taxon import SOURCE_CHOICES, SECTION_CHOICES
from bit.modules.taxonomy.exclusion_list import exclusion_list_help


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        This program downloads assembly files for NCBI genomes. Targets can be pulled
        based on taxonomy (`-t`) from either GTDB or NCBI (also see the --derep-rank parameter),
        and/or they can be explicitly specified as assembly accessions in a file (passed to `-w`).
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
    # every `-t` knob lives in this one group, including the ones only the
    # detailed menu shows -- they are all taxon-selection parameters, so splitting
    # them into a separate "fine-tuning" section put related flags in two places
    selection = parser.add_argument_group("Taxon-selection Parameters (used with `-t`)")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-t",
        "--target-tax",
        metavar="<STR>",
        dest="target_taxon",
        action="append",
        default=None,
        help=wrap_help("Target tax to pull genomes for (a name, or 'all'). Can be given "
                       "multiple times, and can be combined with `-w`."),
    )

    required.add_argument(
        "-w",
        "--wanted-accessions",
        metavar="<FILE>",
        help=wrap_help("Single-column file of wanted NCBI assembly accessions"),
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
        help=wrap_help("Which part of NCBI to draw from (default: both). You probably only "
                       "need to worry about changing this from 'both' if you are setting "
                       "`--derep-rank off` and/or targeting a single species. Ignored with "
                       "`--source gtdb`."),
    )

    selection.add_argument(
        "--derep-rank",
        choices=["auto", "off"] + list(RANKS),
        default="auto",
        help=wrap_help("Dereplicate down to one genome per unique value of this rank "
                       "(default: auto). 'auto' is two ranks finer than the "
                       "target (one rank finer for eukaryotes). 'off' downloads every "
                       "genome under the target taxon, so use with care :)"),
    )

    selection.add_argument(
        "--target-rank",
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
        "--exclusion-list",
        metavar="<FILE>",
        dest="exclusion_list",
        default=None,
        help=h(exclusion_list_help("-t")),
    )

    selection.add_argument(
        "--representatives-only",
        action="store_true",
        help=h("With `--source gtdb`, only pull GTDB representative genomes; "
               "with `--source ncbi`, only pull NCBI reference genomes. If the goal is "
               "removing redundancy, the `--derep-rank` parameter can handle that while "
               "ensuring the breadth of available diversity is maintained."),
    )

    selection.add_argument(
        "--assembly-level",
        action="append",
        choices=list(ASSEMBLY_LEVELS),
        default=None,
        help=h("Only include genomes (from `-t`) at these assembly levels. "
               "Can be provided multiple times."),
    )

    selection.add_argument(
        "--min-completeness",
        metavar="<FLOAT>",
        type=float,
        default=None,
        help=h("Don't include any genomes (from `-t`) below this checkm completeness "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
    )

    selection.add_argument(
        "--max-contamination",
        metavar="<FLOAT>",
        type=float,
        default=None,
        help=h("Don't include any genomes (from `-t`) above this checkm contamination "
               "(default: None). If set, genomes with no recorded "
               "value are also excluded."),
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
        help=wrap_help("Number of concurrent downloads (default: 10; capped at "
                       "20 to keep NCBI happy)"),
        default=10,
        type=int,
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help=wrap_help("Directory for output files (default: current directory)"),
        default=".",
    )

    optional.add_argument(
        "-H",
        "--detailed-help",
        action="store_true",
        help=wrap_help("Show detailed help, including additional taxon-selection "
                       "parameters"),
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show basic help",
    )

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
