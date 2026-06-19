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
        which can then be passed to, e.g., `bit dl-ncbi-assemblies`.
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
        help=wrap_help("Target to search for. Can be a taxon name or an "
                       "NCBI taxid (e.g. 'Alteromonas', 'archaea', or '226')."),
        action="store",
        required=True
    )

    optional.add_argument(
        "-s",
        "--source",
        help=wrap_help("Specify source (default: \"refseq\"; note that 'both' "
                       "will also grab 'suppressed' assembly accessions)"),
        choices=["refseq", "genbank", "both"],
        default="refseq",
    )

    optional.add_argument(
        "-r",
        "--reference-genomes-only",
        action="store_true",
        help=wrap_help("Add this flag to only pull accessions "
                       "for genomes designated as RefSeq \"reference\" genomes (these "
                       "used to be called \"representative\" genomes, see e.g.: "
                       "https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes)"),
    )

    optional.add_argument(
        "--assembly-level",
        nargs="+",
        choices=["complete", "chromosome", "scaffold", "contig"],
        help=wrap_help("Limit to one or more assembly levels, space-separated (e.g. "
                       "'complete chromosome')."),
        action="store",
    )

    optional.add_argument(
        "--annotated-only",
        action="store_true",
        help=wrap_help("Add this flag to limit results to annotated genomes only."),
    )

    optional.add_argument(
        "--get-taxon-counts",
        action="store_true",
        help=wrap_help("Add this flag along with a specified taxon to the `-t` parameter to "
                       "see how many genomes match (under any filters also provided) without "
                       "writing output files."),
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

    if args.reference_genomes_only and args.source == "genbank":
        report_message("The `--reference-genomes-only` flag is incompatible with " 
                       "`--source genbank` as it restricts the search to refseq only.", 
                       trailing_newline=True)
        sys.exit()
