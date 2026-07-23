import sys
import argparse
from bit.modules.general import report_message
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.ncbi.get_accessions_from_ncbi import (get_accessions_from_ncbi,
                                                       _ASSEMBLY_LEVELS)
from bit.modules.taxonomy.tax_ranks import RANKS


def build_parser(parent_subparsers=None):

    desc = ("This is a helper program to facilitate getting NCBI accessions and assembly "
            "metadata based on an NCBI-taxonomy search. It primarily returns NCBI "
            "accessions and metadata subsets based on NCBI-taxonomy searches, with optional "
            "filtering by source (RefSeq/GenBank), assembly level, and/or RefSeq 'reference' genomes "
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

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=("Target taxon (a name, an NCBI taxid, or 'all'). Not needed with "
              "`--get-rank-counts`."),
        action="store",
    )

    optional.add_argument(
        "-r",
        "--target-rank",
        choices=list(RANKS),
        help=("Target rank (if needed to disambiguate a taxon name that exists at multiple ranks)"),
        action="store",
    )

    optional.add_argument(
        "--derep-rank",
        choices=["auto", "off"] + list(RANKS),
        default="off",
        help=("Dereplicate the pulled genomes down to a single best genome per unique "
              "value of this rank (default: off). E.g., '--derep-rank family' keeps one genome per "
              "family within the target taxon). Use 'auto' for two ranks finer than the target. "
              "Only applies to a taxon-name search (not a taxid or 'all')."),
        action="store",
    )

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
        choices=list(_ASSEMBLY_LEVELS),
        nargs="+",
        help=("Restrict to one or more assembly levels (can be multiple space-separated)"),
        action="store",
    )

    optional.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        dest="refseq_reference_genomes_only",
        action="store_true",
        help=("Pull only genomes designated as RefSeq reference genomes."),
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=("Provide this flag along with a specified taxon to `-t` to see how many "
              "genomes match the set parameters (excluding --derep-rank)"),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=("Provide just this flag alone to see counts of how many unique taxa there "
              "are for each rank."),
    )

    optional.add_argument(
        "--get-table",
        action="store_true",
        help=("Provide just this flag alone to write out a tsv of GToTree's "
              "NCBI assembly-summary metadata table."),
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
