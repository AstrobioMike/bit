import sys
import argparse
from bit.modules.general import report_message
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.ncbi.get_accessions_from_ncbi import get_accessions_from_ncbi
from bit.modules.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS,
                                                  add_common_get_accs_args,
                                                  apply_derep_default,
                                                  is_derep_on)


def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate retrieving assembly accessions from NCBI. "
            "It returns NCBI accessions and metadata subsets based on NCBI-taxonomy or taxid searches, "
            "with optional filtering and dereplication settings.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "get-accs-from-ncbi",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit get-accs-from-ncbi -t Alteromonas`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    add_common_get_accs_args(
        required, optional, "NCBI assembly-summary",
        taxon_help=("Target tax to get accessions for (a name, taxid, or 'all'). Not needed with `--get-rank-counts`.")
    )

    optional.add_argument(
        "--ncbi-section",
        dest="ncbi_section",
        type=str.lower,
        default="both",
        choices=["refseq", "genbank", "both"],
        help=("Which part of NCBI to draw from (default: both). You probably only "
              "need to worry about changing this from 'both' if you are setting "
              "`--derep-rank off` and/or targeting a single species."),
    )

    optional.add_argument(
        "--assembly-level",
        action="append",
        choices=list(ASSEMBLY_LEVELS),
        default=None,
        help=("Only include genomes (from `-w`) at this assembly level. "
               "Can be provided multiple times."),
    )

    add_help(optional)
    add_version_arg(optional)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    preflight_checks(args)

    get_accessions_from_ncbi(args)


def preflight_checks(args):

    if not args.get_rank_counts and not args.get_table and not args.target_taxon:
        report_message("A target must be provided to `-t` (a taxon name or an NCBI "
                       "taxid), unless you're using `--get-rank-counts` or "
                       "`--get-table`.",
                       trailing_newline=True)
        sys.exit()

    check_derep_rank_is_applicable(args)


def check_derep_rank_is_applicable(args):
    """
    A taxid's rank is known only after the lookup, and the taxid path doesn't run
    through the selection core, so there is no resolved rank to group beneath. So it's
    not applicable with --derep-rank

    What happens depends on where the value came from: an explicit `--derep-rank`
    is a request we can't honor when given a taxid, so it's an error, while the default
    is just an inherited setting and quietly steps aside so a plain taxid pull still
    works.
    """
    set_explicitly = apply_derep_default(args)

    target = str(args.target_taxon or "")

    if not is_derep_on(args.derep_rank) or not target.isdigit():
        return

    if not set_explicitly:
        args.derep_rank = "off"
        return

    report_message(
        "`--derep-rank` can't be applied with a taxid. Pass the taxon "
        "name instead if you also want to dereplicate.",
        trailing_newline=True)
    sys.exit()
