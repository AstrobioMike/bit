import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, wrap_help, add_version_arg
from bit.modules.ncbi.get_accessions_from_ncbi import get_accessions_from_ncbi


def build_parser(parent_subparsers=None):

    desc = """
        This is a helper program to facilitate getting NCBI accessions and assembly
        metadata based on an NCBI-taxonomy search. It primarily returns NCBI
        accessions and a summary metadata table based on NCBI-taxonomy or taxid searches, 
        which can then be passed to, e.g., `bit dl-ncbi-assemblies`. Filters are available 
        for assembly source (GenBank/RefSeq), RefSeq "reference" genomes only, 
        assembly level, and annotation status.
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
            epilog="Ex. usage: `bit get-accs-from-ncbi -t Alteromonas --assembly-source RefSeq`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-t",
        "--target-taxon",
        metavar="<STR>",
        help=wrap_help("Target taxon to search for. Can be a scientific name or common "
                       "name, or an NCBI taxid (e.g. 'Alteromonas', "
                       "'archaea', or '226')."),
        action="store",
    )

    optional.add_argument(
        "--source",
        help=wrap_help("Specify source (default: \"refseq\"; note that 'both' "
                       "will also grab 'suppressed' assembly accessions)"),
        choices=["genbank", "refseq", "both"],
        default="refseq",
    )

    optional.add_argument(
        "--reference-only",
        action="store_true",
        help=wrap_help("Add this flag to only pull accessions "
                       "for genomes designated as RefSeq \"reference\" genomes (these "
                       "used to be called \"representative\" genomes, see e.g.: "
                       "https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes)"),
    )

    optional.add_argument(
        "--assembly-level",
        metavar="<STR>",
        help=wrap_help("Limit to one or more assembly levels, comma-separated. Options are "
                       "'chromosome', 'complete', 'contig', and 'scaffold' (e.g. "
                       "'complete,chromosome')."),
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

    get_accessions_from_ncbi(args)
