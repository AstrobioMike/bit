import sys
import argparse
from bit.modules.general import report_message
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.ncbi.get_accessions_from_ncbi import get_accessions_from_ncbi
from bit.modules.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS,
                                                  add_common_get_accs_args)


def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate getting NCBI accessions and assembly "
            "metadata based on an NCBI-taxonomy search. It has optional filtering "
            "by source (RefSeq/GenBank), assembly level, and/or RefSeq 'reference' genomes "
            "only, plus optional dereplication down to one genome per specified rank.")

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
        taxon_flags=("-t", "--target-taxon"),
        taxon_help=("Target taxon (a name, an NCBI taxid, or 'all' for every domain "
                    "in the table). Not needed with `--get-rank-counts`."))

    optional.add_argument(
        "-s",
        "--source",
        default="refseq",
        choices=["refseq", "genbank", "both"],
        help=("Specify which section of NCBI to pull from (default: refseq)"),
        action="store",
    )

    optional.add_argument(
        "-a",
        "--assembly-level",
        choices=list(ASSEMBLY_LEVELS),
        nargs="+",
        help=("Restrict to one or more assembly levels (can be multiple space-separated)"),
        action="store",
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

    if args.refseq_reference_genomes_only and args.source != "refseq":
        report_message("The `--reference-genomes-only` flag is only compatible with "
                       "`--source refseq`.",
                       trailing_newline=True)
        sys.exit()

    check_derep_rank_is_applicable(args)


def check_derep_rank_is_applicable(args):
    """
    `--derep-rank` is not applicable to a taxid input
    """
    target = str(args.target_taxon or "")
    derep_rank = getattr(args, "derep_rank", "off")

    if derep_rank in (None, "off", "none", "None"):
        return

    if not target.isdigit():
        return

    report_message(
        "`--derep-rank` can't be applied with a taxid. Pass the taxon "
        "name instead if you also want to dereplicate.",
        trailing_newline=True)
    sys.exit()
