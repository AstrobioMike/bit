import os
import sys
import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore
from bit.modules.general import color_text, wprint, report_message
from bit.modules.ncbi.get_ncbi_assembly_data import (
    check_ncbi_assembly_info_location_var_is_set,
    get_ncbi_assembly_data,
    read_date_retrieved,
)
from bit.modules.ncbi.build_ncbi_data_parquet import PARQUET_FILENAME
from bit.modules.taxonomy.tax_select import (
    AmbiguousTaxon,
    TaxonNotFound,
    resolve_taxon,
    select,
)


_COLUMNS = [
    "assembly_accession",
    "organism_name",
    "taxid",
    "asm_name",
    "assembly_level",
    "refseq_category",
    "checkm_completeness",
    "checkm_contamination",
    "genome_size",
    "domain", "phylum", "class", "order", "family", "genus", "species",
]

_ASSEMBLY_LEVELS = {
    "complete": "Complete Genome",
    "chromosome": "Chromosome",
    "scaffold": "Scaffold",
    "contig": "Contig",
}


def parse_assembly_levels(value):
    if not value:
        return []
    if isinstance(value, str):
        parts = [v.strip().lower() for v in value.split(",") if v.strip()]
    else:
        parts = [str(v).strip().lower() for v in value]

    unknown = [p for p in parts if p not in _ASSEMBLY_LEVELS]
    if unknown:
        raise ValueError(
            f"unrecognised --assembly-level value(s): {', '.join(unknown)}. "
            f"Choose from: {', '.join(_ASSEMBLY_LEVELS)}")
    return [_ASSEMBLY_LEVELS[p] for p in parts]


def ncbi_table_path(force_update=False, quiet=True):
    get_ncbi_assembly_data(force_update=force_update, quiet=quiet)
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    return os.path.join(ncbi_dir, PARQUET_FILENAME)


def get_accessions_from_ncbi(args):

    table_path = ncbi_table_path()
    _report_ncbi_date(table_path)

    if getattr(args, "get_rank_counts", False):
        report_unique_taxa_counts_of_all_ranks(
            table_path, source=args.source,
            reps_only=args.refseq_reference_genomes_only)
        sys.exit(0)

    if str(args.target_taxon).lower() == "all":
        tab = _select_all(table_path, reps_only=args.refseq_reference_genomes_only)
        label = "all genomes"
    elif str(args.target_taxon).isdigit():
        tab = _select_by_taxid(table_path, str(args.target_taxon))
        label = f"taxid {args.target_taxon}"
    else:
        try:
            canonical, rank = resolve_taxon(table_path, args.target_taxon,
                                            rank=getattr(args, "target_rank", None))
        except AmbiguousTaxon as e:
            report_message(
                f"'{e.taxon}' occurs at more than one rank "
                f"({', '.join(e.ranks_found)}). Specify which with `-r`, or pass "
                f"the NCBI taxid to `-t` instead.", "yellow")
            print("")
            sys.exit(0)
        except TaxonNotFound:
            report_message(f"Input taxon '{args.target_taxon}' doesn't seem to "
                           f"exist at any rank :(", "yellow")
            print("")
            sys.exit(0)

        if canonical != args.target_taxon:
            wprint(color_text(f"Matched input '{args.target_taxon}' to NCBI taxon "
                              f"'{canonical}'.", "yellow"))
            print("")

        tab = select(table_path, "ncbi", rank, canonical,
                     reps_only=args.refseq_reference_genomes_only, columns=_COLUMNS)
        label = f"{rank} '{canonical}'"

    tab = _apply_filters(tab, args)

    if args.get_taxon_counts:
        report_message(f"There are {tab.num_rows:,} genome(s) under {label} with "
                       "the specified filters.", "none",
                       initial_indent="    ", subsequent_indent="    ",
                       trailing_newline=True)
        sys.exit(0)

    if tab.num_rows == 0:
        report_message(f"No genomes were found under {label} with the specified "
                       "filters.", "none",
                       initial_indent="    ", subsequent_indent="    ")
        print("")
        sys.exit(0)

    _write_outputs(tab, args, label)


def report_unique_taxa_counts_of_all_ranks(table_path, source="refseq", reps_only=False):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the NCBI table
    (matching get-accs-from-gtdb's --get-rank-counts), scoped to `source`:
      refseq  -> GCF_ accessions only
      genbank -> GCA_ accessions only
      both    -> no source filter
    If reps_only, also print the counts among RefSeq reference genomes (which are
    RefSeq, so that sub-table is implicitly refseq-scoped).

    Reads the rank columns plus assembly_accession (for the source prefix filter),
    never the whole table. Distinct counts include the "NA" placeholder as one value,
    consistent with the GTDB command.
    """
    from bit.modules.taxonomy.tax_ranks import RANKS

    ranks = list(RANKS)
    tab = pq.read_table(table_path, columns=ranks + ["assembly_accession"])
    tab = _filter_by_source(tab, source)

    label = {"refseq": "RefSeq", "genbank": "GenBank", "both": "all"}.get(source, source)
    print("\n    {:<10} {:}".format("Rank", f"Num. Unique Taxa ({label})"))
    for rank in ranks:
        n = pc.count_distinct(tab.column(rank)).as_py()
        print("    {:<10} {:}".format(rank, str(n)))
    print("")

    if reps_only:
        rep = pq.read_table(table_path, columns=ranks,
                            filters=[("refseq_category", "=", "reference genome")])
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print("    {:<10} {:}".format("Rank", "Num. Unique Ref. Taxa"))
        for rank in ranks:
            n = pc.count_distinct(rep.column(rank)).as_py()
            print("    {:<10} {:}".format(rank, str(n)))
        print("")


def _filter_by_source(tab, source):
    """
    Scope a table to a source by accession prefix (refseq -> GCF_, genbank -> GCA_,
    both -> unfiltered)
    """
    if source in ("refseq", "genbank"):
        prefix = "GCF_" if source == "refseq" else "GCA_"
        tab = tab.filter(pc.starts_with(tab.column("assembly_accession"), prefix))
    return tab


def _report_ncbi_date(table_path):
    date_str = read_date_retrieved(os.path.dirname(table_path))
    print("\n    Date NCBI data retrieved: " + date_str)


def _select_all(table_path, reps_only=False):
    filters = [("refseq_category", "=", "reference genome")] if reps_only else None
    return pq.read_table(table_path, columns=_COLUMNS, filters=filters)


def _select_by_taxid(table_path, taxid):
    from bit.modules.taxonomy.tax_ranks import RANKS

    for rank in RANKS:
        tab = pq.read_table(table_path, columns=_COLUMNS,
                            filters=[(f"{rank}_taxid", "=", str(taxid))])
        if tab.num_rows:
            return tab

    return pq.read_table(table_path, columns=_COLUMNS,
                         filters=[("taxid", "=", str(taxid))])


def _apply_filters(tab, args):
    tab = _filter_by_source(tab, getattr(args, "source", "refseq"))

    wanted = parse_assembly_levels(getattr(args, "assembly_level", None))
    if wanted:
        tab = tab.filter(pc.is_in(tab.column("assembly_level"),
                                  value_set=pa.array(wanted, type=pa.string())))
    return tab


def _write_outputs(tab, args, label):
    taxon_raw = str(args.target_taxon)
    if taxon_raw.lower() == "all":
        taxon_raw = "all"
    taxon_for_filename = taxon_raw.replace(" ", "-").replace("/", "-")

    suffix_bits = []
    if args.refseq_reference_genomes_only:
        suffix_bits.append("refseq-ref")
    elif getattr(args, "source", "refseq") != "both":
        suffix_bits.append(args.source.lower())
    suffix = ("-" + "-".join(suffix_bits)) if suffix_bits else ""

    acc_out = f"NCBI-{taxon_for_filename}{suffix}-accessions.txt"
    tab_out = f"NCBI-{taxon_for_filename}{suffix}-metadata.tsv"

    df = tab.to_pandas()
    df.to_csv(tab_out, sep="\t", index=False)

    accs = df["assembly_accession"].tolist()
    with open(acc_out, "w") as out:
        for acc in accs:
            out.write(acc + "\n")

    print("")
    wprint(f"Wrote {len(accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out))
    print("")
