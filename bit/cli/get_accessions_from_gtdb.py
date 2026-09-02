import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.general import report_message
from bit.modules.gtdb.get_accessions_from_gtdb import get_accessions_from_gtdb
from bit.modules.taxonomy.get_accs_shared import (add_common_get_accs_args,
                                                  apply_derep_default)


def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate retrieving accessions from the "
            "Genome Taxonomy Database (gtdb.ecogenomic.org). It returns NCBI "
            "accessions and GTDB metadata subsets based on GTDB-taxonomy searches, "
            "with optional filtering and dereplication settings.")

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "get-accs-from-gtdb",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit get-accs-from-gtdb -t Archaea --gtdb-representatives-only`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    add_common_get_accs_args(
        required, optional, "GTDB")

    optional.add_argument(
        "-G",
        "--gtdb-representatives-only",
        action="store_true",
        help=("Pull only genomes designated as GTDB species representatives."),
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

    get_accessions_from_gtdb(args)


def preflight_checks(args):

    apply_derep_default(args)

    if not args.get_rank_counts and not args.get_table and not args.target_taxon:
        report_message("A target must be provided to `-t` (a taxon name), "
                       "unless you're using `--get-rank-counts` or `--get-table`.",
                       trailing_newline=True)
        sys.exit()

