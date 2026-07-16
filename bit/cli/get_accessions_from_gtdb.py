import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.general import report_message
from bit.modules.gtdb.get_accessions_from_gtdb import get_accessions_from_gtdb


def build_parser(parent_subparsers=None):

    desc = """
        This is a helper program to facilitate using taxonomy and genomes from
        the Genome Taxonomy Database (gtdb.ecogenomic.org). It primarily returns
        NCBI accessions and GTDB summary tables based on GTDB-taxonomy searches,
        which could then be passed to, e.g., `bit dl-ncbi-assemblies`. bit
        caches the GTDB metadata. If you want to update it, run
        `bit data get gtdb-data -f`.
        """

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

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=wrap_help("Target taxon (enter 'all' for all). "
                       "Not needed with `--get-rank-counts`."),
        action="store",
    )

    optional.add_argument(
        "-r",
        "--target-rank",
        metavar="<STR>",
        help=wrap_help("Target rank (if needed to disambiguate a taxon name that exists at multiple ranks)"),
        action="store",
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=wrap_help("Add this flag along with a specified taxon to the `-t` parameter "
                       "to see how many of that taxon are in the database."),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=wrap_help("Provide just this flag alone to see counts of how many "
                       "unique taxa there are for each rank."),
    )

    optional.add_argument(
        "-G",
        "--gtdb-representatives-only",
        action="store_true",
        help=wrap_help("Add this flag to only pull accessions for genomes "
                       "designated as GTDB species representatives (see, e.g., "
                       "https://gtdb.ecogenomic.org/faq#gtdb_species_clusters)."),
    )

    optional.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        action="store_true",
        help=wrap_help("Add this flag to only pull accessions for genomes designated as "
                       "RefSeq \"reference\" genomes (these used to be called \"representative\" genomes, see, e.g., "
                       "https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes)."),
    )

    optional.add_argument(
        "--get-table",
        action="store_true",
        help=wrap_help("Provide just this flag alone to write out a tsv of bit's GTDB metadata "
                       "table."),
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

    if not args.get_rank_counts and not args.target_taxon:
        report_message("A target must be provided to `-t` (a taxon name), "
                       "unless you're using `--get-rank-counts`.",
                       trailing_newline=True)
        sys.exit()
