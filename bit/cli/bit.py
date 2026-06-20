import sys
import argparse
import importlib
from bit.cli.common import add_help, add_version_arg, CustomRichHelpFormatter


SUBCOMMAND_MAP = {
    "get-accs-from-gtdb":  "bit.cli.get_accessions_from_gtdb",
    "get-accs-from-ncbi":  "bit.cli.get_accessions_from_ncbi",
    "dl-ncbi-assemblies":  "bit.cli.dl_ncbi_assemblies",
    "gen-metagenome":      "bit.cli.gen_metagenome",
    "gen-reads":           "bit.cli.gen_reads",
    "mutate-seqs":         "bit.cli.mutate_seqs",
    "add-insertion":       "bit.cli.add_insertion",
    "cov-analyzer":        "bit.cli.cov_analyzer",
    "cov-stats":           "bit.cli.cov_stats",
    "mapped-read-stats":   "bit.cli.mapped_read_stats",
    "aa-diff":             "bit.cli.aa_diff",
    "ez-screen":           "bit.cli.ez_screen",
    "fasta":               "bit.cli.fasta",
    "assemble":            "bit.cli.assemble",
    "summarize-assembly":  "bit.cli.summarize_assembly",
    "genbank":             "bit.cli.genbank",
    "kraken2":             "bit.cli.kraken2",
    "lineage":             "bit.cli.lineage",
    "table":               "bit.cli.table",
    "filter-ko-results":   "bit.cli.filter_kofamscan_results",
    "go":                  "bit.cli.go",
    "itol":                "bit.cli.itol",
    "get-workflow":        "bit.cli.get_workflow",
    "data":                "bit.cli.data",
}


PROGRAM_GROUPS = [
    {
        "title": "NCBI/GTDB-related",
        "programs": [
            {
                "name": "get-accs-from-gtdb",
                "desc": "search the GTDB by taxonomy and retrieve NCBI accessions",
            },
            {
                "name": "get-accs-from-ncbi",
                "desc": "search NCBI by taxonomy or taxid and retrieve NCBI accessions",
            },
            {
                "name": "dl-ncbi-assemblies",
                "desc": "download NCBI assemblies in different formats given input accessions",
            },
        ],
    },
    {
        "title": "Simulation / sequence manipulation",
        "programs": [
            {
                "name": "gen-metagenome",
                "desc": "generate an in silico metagenome with ground-truth tables",
            },
            {
                "name": "gen-reads",
                "desc": "generate reads from fasta files",
            },
            {
                "name": "mutate-seqs",
                "desc": "introduce point mutations into nucleotide or amino-acid fasta files",
            },
            {
                "name": "add-insertion",
                "desc": "add insertions into nucleotide or amino-acid fasta sequences",
            },
        ],
    },
    {
        "title": "Coverage/mapping-related",
        "programs": [
            {
                "name": "cov-analyzer",
                "desc": "analyze coverage patterns to identify regions of relatively higher/lower coverage",
            },
            {
                "name": "cov-stats",
                "desc": "get detection, coverage, and mean percent ID for references given fasta(s) and a bam",
            },
            {
                "name": "mapped-read-stats",
                "desc": "get percent ID and other information for mapped reads in a bam",
            },
        ],
    },
    {
        "title": "Sequence searching/comparing",
        "programs": [
            {
                "name": "aa-diff",
                "desc": "compare a query sequence to an amino-acid reference and report differences",
            },
            {
                "name": "ez-screen",
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
                "name": "fasta",
                "desc": "",
                "subcommands": [
                    ("calc-gc",           "calculate GC content per sequence or for the full file"),
                    ("calc-var-in-msa",   "calculate variation in each column of a multiple-sequence alignment"),
                    ("count",             "count and summarize bases or sequences"),
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
                "name": "assemble",
                "desc": "simple wrapper for assembly with optional quality trimming and normalization",
            },
            {
                "name": "summarize-assembly",
                "desc": "quickly summarize nucleotide assemblies",
            },
        ],
    },
    {
        "title": "GenBank-format utilities",
        "programs": [
            {
                "name": "genbank",
                "desc": "",
                "subcommands": [
                    ("to-AA-seqs", "extract amino acid sequences"),
                    ("to-cds-tsv", "extract CDS info to a TSV"),
                    ("to-cds-seqs","extract CDS nucleotide sequences"),
                    ("to-fasta",   "extract nucleotide sequences"),
                ],
            },
        ],
    },
    {
        "title": "Taxonomy and lineage helpers",
        "programs": [
            {
                "name": "kraken2",
                "desc": "",
                "subcommands": [
                    ("tax-summary","generate summary tables from kraken2 or bracken outputs"),
                    ("tax-plots",  "generate standard taxonomy barplots from kraken2 or bracken outputs"),
                ],
            },
            {
                "name": "lineage",
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
                "name": "table",
                "desc": "",
                "subcommands": [
                    ("colnames",         "print column names with numbers (handy for cut/awk)"),
                    ("filter",           "filter a table based on wanted strings"),
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
                "name": "filter-ko-results",
                "desc": "filter KOFamScan results",
            },
            {
                "name": "go",
                "desc": "",
                "subcommands": [
                    ("get-term-info",        "print out GO term info"),
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
                "name": "itol",
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
                "name": "get-workflow",
                "desc": "download a bit-packaged workflow",
            },
        ],
    },
    {
        "title": "bit-data management",
        "programs": [
            {
                "name": "data",
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
    console.print(f"{'':>27}OVERVIEW OF SUBCOMMANDS")
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



def _suppress_help_version_on_group_parsers(parser):
    """
    Suppress -h/--help/-v/--version from argcomplete on parsers that have
    subparsers, so TAB after a group shows only subcommand names.
    Leaf-level parsers are left untouched so their flags still appear
    without requiring a '-' prefix.
    """
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for a in parser._actions:
                if hasattr(a, 'option_strings') and any(
                    s in ('-h', '--help', '-v', '--version')
                    for s in a.option_strings
                ):
                    a.help = argparse.SUPPRESS
            for sub_parser in action.choices.values():
                _suppress_help_version_on_group_parsers(sub_parser)
            break


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

    subparsers = parser.add_subparsers(dest='subcommand')

    # Lazy-load modules during tab completion so only the one relevant module
    # is imported per keypress instead of all 21 at once.
    import os
    comp_line = os.environ.get('COMP_LINE', '')
    if comp_line:
        words = comp_line.split()
        if len(words) >= 2 and words[1] in SUBCOMMAND_MAP:
            # User is completing args/flags for a specific subcommand — only
            # import that one module.
            module = importlib.import_module(SUBCOMMAND_MAP[words[1]])
            module.build_parser(parent_subparsers=subparsers)
        else:
            # User is still completing the subcommand name itself — lightweight
            # stubs are all argcomplete needs to offer the names.
            for name in SUBCOMMAND_MAP:
                subparsers.add_parser(name, add_help=False)
    else:
        # Normal (non-completion) invocation — build the full tree.
        for module_path in SUBCOMMAND_MAP.values():
            module = importlib.import_module(module_path)
            module.build_parser(parent_subparsers=subparsers)

    _suppress_help_version_on_group_parsers(parser)

    return parser


def main():

    parser = build_parser()

    try:
        import argcomplete # type: ignore
        argcomplete.autocomplete(parser)
    except ImportError:
        pass

    # print the overview when called with no arguments or with -h/--help
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print_overview()
        sys.exit(0)

    subcommand = sys.argv[1]

    # handle -v/--version before subcommand dispatch
    if subcommand in ("-v", "--version"):
        parser.parse_args()
        return

    if subcommand not in SUBCOMMAND_MAP:
        from rich.console import Console # type: ignore
        Console(stderr=True).print(f"\n    [yellow]Unknown subcommand:[/yellow] [cyan]{subcommand}[/cyan]\n")
        Console(stderr=True).print(f"  Run [cyan]bit[/cyan] by itself to see available subcommands.\n")
        # print_overview()
        sys.exit(1)

    # rewrite argv so the target module's parser sees the right program name
    sys.argv = [f"bit {subcommand}"] + sys.argv[2:]
    module = importlib.import_module(SUBCOMMAND_MAP[subcommand])
    module.main()
