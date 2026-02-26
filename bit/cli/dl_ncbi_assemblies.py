import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.dl_ncbi_assemblies import dl_ncbi_assemblies

def main():

    desc = """
        This program downloads assembly files for NCBI genomes. It takes as input
        assembly accessions (either GCA_* or GCF_*) and optionally a specification of
        which format to download. For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: `bit-dl-ncbi-assemblies -w wanted-accessions.txt -f fasta -j 4`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("Required Parameters")
    optional = parser.add_argument_group("Optional Parameters")

    required.add_argument(
        "-w",
        "--wanted-accessions",
        metavar="<FILE>",
        help='Input file with wanted accessions',
        required=True,
    )

    optional.add_argument(
        "-f",
        "--format",
        help='Format to download (default: "fasta")',
        default="fasta",
        choices=["fasta", "protein", "genbank", "gff", "nt_cds", "feature_tab", "report", "stats"],
    )

    optional.add_argument(
        "-j",
        "--jobs",
        metavar="<INT>",
        help="Number of downloads you'd like to run concurrently. NCBI can become unhappy with many requests,\
             so a max of 20 will be used even if more are requested (default: 10)",
        default=10,
        type=int,
    )

    optional.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help="Directory to save output files (default: current directory)",
        default=".",
    )

    add_help(optional)

    if len(sys.argv) == 1:  # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_ncbi_assemblies(args)
