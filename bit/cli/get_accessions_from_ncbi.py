import sys
import argparse
from bit.modules.general import report_message
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.ncbi.get_accessions_from_ncbi import get_accessions_from_ncbi


def build_parser(parent_subparsers=None):

    desc = """
        This is a helper program to facilitate getting NCBI accessions and assembly
        metadata based on an NCBI-taxonomy search. It primarily returns NCBI
        accessions and a summary metadata table based on NCBI-taxonomy or taxid searches, 
        which can then be passed to, e.g., `bit dl-ncbi-assemblies`. bit 
        caches the NCBI assembly metadata. If you want to update it, run
        `bit data get ncbi-assembly-data -f`.
        """

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
        help=wrap_help("Target taxon (enter 'all' for all). Can be a taxon name or an "
                       "NCBI taxid (e.g., 'Alteromonas', 'archaea', or '226'). "
                       "Not needed with `--get-rank-counts`."),
        action="store",
    )

    optional.add_argument(
        "-r",
        "--target-rank",
        metavar="<STR>",
        help=wrap_help("If the target taxon name occurs at more than one rank, "
                       "specify which rank is wanted (e.g. 'genus'). Only needed to "
                       "disambiguate a homonym."),
        action="store",
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=wrap_help("Add this flag along with a specified taxon to the `-t` parameter to "
                       "see how many genomes match (under any filters also provided) without "
                       "writing output files."),
    )

    optional.add_argument(
        "--get-rank-counts",
        action="store_true",
        help=wrap_help("Provide just this flag alone to see counts of how many "
                       "unique taxa there are for each rank."),
    )

    optional.add_argument(
        "-s",
        "--source",
        help=wrap_help("Specify source (default: \"refseq\"; note that 'both' or 'genbank' "
                       "will also grab 'suppressed' assembly accessions)"),
        choices=["refseq", "genbank", "both"],
        default="refseq",
    )

    optional.add_argument(
        "-R",
        "--refseq-reference-genomes-only",
        action="store_true",
        help=wrap_help("Add this flag to only pull accessions "
                       "for genomes designated as RefSeq \"reference\" genomes (these "
                       "used to be called \"representative\" genomes, see, e.g., "
                       "https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes)"),
    )

    optional.add_argument(
        "--assembly-level",
        nargs="+",
        choices=["complete", "chromosome", "scaffold", "contig"],
        help=wrap_help("Limit to one or more assembly levels, space-separated (e.g., "
                       "'complete chromosome')."),
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

    if not args.get_rank_counts and not args.target_taxon:
        report_message("A target must be provided to `-t` (a taxon name or an NCBI "
                       "taxid), unless you're using `--get-rank-counts`.",
                       trailing_newline=True)
        sys.exit()

    if args.refseq_reference_genomes_only and args.source != "refseq":
        report_message("The `--reference-genomes-only` flag is only compatible with " 
                       "`--source refseq`.", 
                       trailing_newline=True)
        sys.exit()
