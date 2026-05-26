import sys
import argparse
from bit.cli.common import add_help, add_version_arg, CustomRichHelpFormatter


PROGRAM_GROUPS = [
    {
        "title": "GTDB/NCBI-related",
        "programs": [
            {
                "name": "bit-dl-ncbi-assemblies",
                "desc": "download NCBI assemblies in different formats given input accessions",
            },
            {
                "name": "bit-get-accs-from-gtdb",
                "desc": "search the GTDB by taxonomy and retrieve NCBI accessions",
            },
        ],
    },
    {
        "title": "Coverage/mapping-related",
        "programs": [
            {
                "name": "bit-cov-analyzer",
                "desc": "analyze coverage patterns from a bam + reference fasta to identify regions of relatively higher or lower coverage",
            },
            {
                "name": "bit-cov-stats",
                "desc": "get detection, coverage, and mean percent ID for single or multiple references given fasta(s) and a bam file",
            },
            {
                "name": "bit-mapped-reads-pid",
                "desc": "get percent ID information for mapped reads in a bam file",
            },
        ],
    },
    {
        "title": "Sequence manipulation / read generation",
        "programs": [
            {
                "name": "bit-gen-reads",
                "desc": "generate reads from fasta files",
            },
            {
                "name": "bit-mutate-seqs",
                "desc": "introduce point mutations (substitutions/indels) into nucleotide or amino-acid fasta files",
            },
            {
                "name": "bit-add-insertion",
                "desc": "add insertions into nucleotide or amino-acid fasta sequences",
            },
        ],
    },
    {
        "title": "Sequence searching",
        "programs": [
            {
                "name": "bit-ez-screen",
                "desc": "",
                "subcommands": [
                    ("assembly",         "run blast-based screening of targets in assemblies"),
                    ("reads",            "run mapping-based screening of targets in reads"),
                ]
            },
        ],
    },
    {
        "title": "Fasta utilities",
        "programs": [
            {
                "name": "bit-fasta",
                "desc": "",
                "subcommands": [
                    ("calc-gc",           "calculate GC content per sequence or for the full file"),
                    ("calc-var-in-msa",   "calculate variation in each column of a multiple-sequence alignment"),
                    ("count",             "count and summarize bases or sequences in a fasta file"),
                    ("extract-by-coords", "extract sequences by genomic coordinates"),
                    ("extract-by-headers","extract sequences by header names"),
                    ("extract-by-primers","extract sequences based on primer sequences"),
                    ("filter-by-length",  "filter sequences by minimum/maximum length"),
                    ("modify-headers",    "rename or reformat sequence headers"),
                    ("remove-wraps",      "remove soft line wraps"),
                    ("to-bed",            "convert fasta to BED format"),
                    ("to-genbank",        "convert fasta to GenBank format"),
                ],
            },
        ],
    },
    {
        "title": "Assembly-related",
        "programs": [
            {
                "name": "bit-assemble",
                "desc": "simple wrapper for assembly with optional quality trimming and normalization",
            },
            {
                "name": "bit-summarize-assembly",
                "desc": "quickly summarize nucleotide assemblies",
            },
        ],
    },
    {
        "title": "GenBank-format utilities",
        "programs": [
            {
                "name": "bit-genbank",
                "desc": "",
                "subcommands": [
                    ("to-fasta",   "extract nucleotide sequences"),
                    ("to-AA-seqs", "extract amino acid sequences"),
                    ("to-cds-tsv", "extract CDS info to a TSV"),
                    ("to-cds-seqs","extract CDS nucleotide sequences"),
                ],
            },
        ],
    },
    {
        "title": "Taxonomy and lineage helpers",
        "programs": [
            {
                "name": "bit-kraken2",
                "desc": "summarize and visualize kraken2/bracken outputs",
                "subcommands": [
                    ("tax-summary","generate summary tables from kraken2/bracken outputs"),
                    ("tax-plots",  "generate standard taxonomy barplots from kraken2/bracken outputs"),
                ],
            },
            {
                "name": "bit-lineage",
                "desc": "",
                "subcommands": [
                    ("from-taxids","get full lineage info from a list of NCBI taxon IDs"),
                    ("to-tsv",     "reformat lineage info to a TSV"),
                ],
            },
        ],
    },
    {
        "title": "Table utilities",
        "programs": [
            {
                "name": "bit-table",
                "desc": "",
                "subcommands": [
                    ("colnames",         "print column names with numbers (handy for cut/awk)"),
                    ("filter",           "filter a table based on wanted IDs"),
                    ("normalize",        "normalize to CPM or with the DESeq2 median-ratio method"),
                    ("summarize-column", "summarize a numeric column"),
                ],
            },
        ],
    },
    {
        "title": "Functional-annotation helpers",
        "programs": [
            {
                "name": "bit-filter-ko-results",
                "desc": "filter KOFamScan results",
            },
            {
                "name": "bit-go",
                "desc": "",
                "subcommands": [
                    ("get-term-info",        "look up GO term info"),
                    ("summarize-annotations","summarize GO annotations"),
                    ("combine-summaries",    "combine GO summary outputs"),
                    ("slim-terms",           "slim GO terms to a specified ontology"),
                ],
            },
        ],
    },
    {
        "title": "iTOL helpers",
        "programs": [
            {
                "name": "bit-itol",
                "desc": "",
                "subcommands": [
                    ("binary-dataset","generate a binary dataset annotation file"),
                    ("colorstrip",    "generate a color strip annotation file"),
                    ("map",           "generate a mapping/connection file"),
                    ("text-dataset",  "generate a text label dataset file"),
                ],
            },
        ],
    },
    {
        "title": "Workflows",
        "programs": [
            {
                "name": "bit-get-workflow",
                "desc": "retrieve a snakemake workflow (available: sra-download, genome-summarize, metagenomics)",
            },
        ],
    },
    {
        "title": "bit-data management",
        "programs": [
            {
                "name": "bit-data",
                "desc": "",
                "subcommands": [
                    ("get",       "download/update bit-utilized databases or get test data"),
                    ("locations", "check or set data-location environment variables"),
                ],
            },
        ],
    },
]


def print_overview():
    from rich.console import Console # type: ignore
    from rich.table import Table # type: ignore
    from importlib.metadata import version
    from datetime import datetime

    console = Console()
    ver = f"v{version('bit')}"

    console.print()
    console.print(f"{'':>33}bit [green]{ver}[/green]")
    console.print(f"{'':>25}github.com/AstrobioMike/bit")
    console.print()
    console.print(f"{'':>28}OVERVIEW OF PROGRAMS")
    console.print()

    # global name-column width — keeps description column aligned across all groups
    name_col_width = max(
        max(
            max(len(p["name"]) for p in g["programs"]),
            max((len(s[0]) + 2 for p in g["programs"] for s in p.get("subcommands", [])), default=0),
        )
        for g in PROGRAM_GROUPS
    ) + 2

    console.rule(style="dim")

    for i, group in enumerate(PROGRAM_GROUPS):
        if i > 0:
            console.print()
            console.rule(style="dim")
        console.print()

        # group header — dark_orange matches argparse.groups style
        console.print(f"  [dark_orange]{group['title']}[/dark_orange]")

        # borderless two-column table: name | description
        # wrapping is handled per-column so descriptions indent cleanly
        table = Table(box=None, show_header=False, padding=(0, 2, 0, 4))
        table.add_column(no_wrap=True, min_width=name_col_width)  # name — fixed width, never wraps
        table.add_column()                                         # description — wraps cleanly

        for prog in group["programs"]:
            name = prog["name"]
            desc = prog["desc"]
            subcommands = prog.get("subcommands")

            table.add_row(
                f"[cyan]{name}[/cyan]",
                desc,
            )

            for sub_name, sub_desc in (subcommands or []):
                table.add_row(
                    f"  [dark_cyan]{sub_name}[/dark_cyan]",
                    sub_desc,
                )

        console.print(table)

    console.print()
    console.rule(style="dim")
    console.print()

    today = datetime.today().strftime('%A')
    signoff = f"Happy {today} :)"
    console.print(f"{'':>51}[green]{signoff}[/green]\n")



def build_parser():

    desc = """
        Print an overview of all available bit programs.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )

    add_help(parser)
    add_version_arg(parser)

    return parser


def main():

    parser = build_parser()

    # print the overview when called with no arguments or with -h/--help
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print_overview()
        sys.exit(0)

    args = parser.parse_args()
