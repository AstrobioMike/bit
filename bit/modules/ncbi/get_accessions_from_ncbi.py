import os
import sys
import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore
from bit.modules.general import (color_text, wprint, report_message,
                                 write_table_tsv, write_accessions)
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
    find_ranks_for_taxon as _resolve_ranks,
)
from bit.modules.taxonomy.tax_derep import (select_ref_genomes, select_all_domains,
                                            resolve_derep_rank)
from bit.modules.taxonomy.tax_counts import (representatives_filter, count_genomes,
                                             derep_size, rank_counts, read_pool,
                                             render_rank_count_table)
from bit.modules.taxonomy.tax_targets import (is_all_target,
                                              unassigned_domain_summary)
from bit.modules.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS, PoolSpec,
                                                  derep_note as _shared_derep_note,
                                                  is_derep_on,
                                                  parse_assembly_levels as _parse_levels,
                                                  scoped_counts_note,
                                                  source_overlap_note,
                                                  source_prefixes as _shared_prefixes)


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

_ASSEMBLY_LEVELS = ASSEMBLY_LEVELS


def parse_assembly_levels(value):
    """Kept as a name here; the body is shared with the GTDB helper's vocabulary."""
    return _parse_levels(value)


def ncbi_table_path(force_update=False, quiet=True):
    get_ncbi_assembly_data(force_update=force_update, quiet=quiet)
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    return os.path.join(ncbi_dir, PARQUET_FILENAME)


def _source_prefixes(source):
    """Accession prefixes for a --source value (None means no restriction)."""
    return _shared_prefixes(source)


def copy_ncbi_table(table_path):
    """
    Write out the parquet object as a tsv
    """
    out_name = "ncbi-assembly-summary-metadata.tsv"
    write_table_tsv(pq.read_table(table_path), out_name)

    print("")
    wprint("NCBI table written to:")
    print(color_text("    " + out_name + "\n"))


def get_accessions_from_ncbi(args):

    table_path = ncbi_table_path()
    _report_ncbi_date(table_path)
    _report_source_overlap(getattr(args, "source", "refseq"))

    if getattr(args, "get_table", False):
        copy_ncbi_table(table_path)
        sys.exit(0)

    assembly_levels = parse_assembly_levels(getattr(args, "assembly_level", None))

    target = str(args.target_taxon) if args.target_taxon else ""
    named_taxon = bool(target) and not is_all_target(target) and not target.isdigit()

    if getattr(args, "get_rank_counts", False):
        if named_taxon:
            _report_rank_counts_for_taxon_or_exit(table_path, target, args,
                                                  assembly_levels)
        else:
            report_unique_taxa_counts_of_all_ranks(
                table_path, source=args.source,
                reps_only=args.refseq_reference_genomes_only,
                assembly_levels=assembly_levels,
                scoped_to_all=is_all_target(target))
        sys.exit(0)

    if getattr(args, "get_taxon_counts", False) and named_taxon:
        _report_taxon_counts_or_exit(table_path, target, args, assembly_levels)
        sys.exit(0)

    if is_all_target(args.target_taxon):
        if _derep_is_on(args):
            tab, label = _select_all_dereplicated(table_path, args, assembly_levels)
        else:
            tab = _select_all(table_path, args,
                              reps_only=args.refseq_reference_genomes_only,
                              assembly_levels=assembly_levels)
            label = "all genomes"
    elif str(args.target_taxon).isdigit():
        tab = _select_by_taxid(table_path, str(args.target_taxon),
                               assembly_levels=assembly_levels)
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

        try:
            selection = select_ref_genomes(
                table_path, "ncbi", canonical, target_rank=rank,
                derep_rank=getattr(args, "derep_rank", "off"),
                reps_only=args.refseq_reference_genomes_only,
                accession_prefixes=_source_prefixes(args.source),
                assembly_levels=assembly_levels)
        except ValueError as e:
            report_message(str(e), "yellow")
            print("")
            sys.exit(0)

        for w in selection.warnings:
            wprint(color_text(w, "yellow"))
        if selection.warnings:
            print("")

        tab = pa.Table.from_pylist(selection.rows) if selection.rows else \
            select(table_path, "ncbi", rank, canonical,
                   reps_only=args.refseq_reference_genomes_only,
                   columns=_COLUMNS).slice(0, 0)
        label = f"{rank} '{canonical}'"
        if selection.effective_derep_rank:
            label += f" (dereplicated to one genome per {selection.effective_derep_rank})"

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


def _prefix_mask(acc_col, prefixes):
    mask = None
    for p in prefixes:
        m = pc.starts_with(acc_col, p)
        mask = m if mask is None else pc.or_(mask, m)
    return mask


def _count_at_rank(table_path, rank, taxon, prefixes=None, reps_only=False,
                   assembly_levels=None):
    """
    Count rows where column `rank` == `taxon`, with the set POOL filters applied
    (source prefix, RefSeq-reference-only, assembly level)
    """
    return count_genomes(table_path, "ncbi", rank=rank, taxon=taxon,
                         rep_filter=_rep_filter(reps_only),
                         accession_prefixes=prefixes,
                         assembly_levels=assembly_levels)


def _rep_filter(reps_only):
    """The RefSeq-reference-genome predicate, or None."""
    return representatives_filter("ncbi", "refseq" if reps_only else None)


def _derep_count_at_rank(table_path, rank, taxon, derep_rank, prefixes=None,
                         reps_only=False, assembly_levels=None):
    """How many genomes survive dereplication at `derep_rank`, under the same pool."""
    return derep_size(table_path, "ncbi", rank, taxon, derep_rank,
                      rep_filter=_rep_filter(reps_only),
                      accession_prefixes=prefixes,
                      assembly_levels=assembly_levels)


def _derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                reps_only=False):
    """
    The "...dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off

    Returns (line_or_None, warnings).
    """
    pool = PoolSpec(table_path, "ncbi", rep_filter=_rep_filter(reps_only),
                    accession_prefixes=prefixes, assembly_levels=assembly_levels,
                    label="NCBI", taxon_flag="-t")
    return _shared_derep_note(pool, rank, taxon,
                              getattr(args, "derep_rank", "off"))


def _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                       reps_only=False):
    """Print the dereplicated-count line (and any 'auto' advisory) for one rank."""
    line, warnings = _derep_note(table_path, rank, taxon, args, prefixes,
                                 assembly_levels, reps_only=reps_only)
    if line:
        report_message(line, color=None, width=100, initial_indent="      ",
                       subsequent_indent="      ", leading_newline=False,
                       trailing_newline=True)
    for warning in warnings:
        report_message(warning, "yellow", width=100, initial_indent="      ",
                       subsequent_indent="      ", leading_newline=False,
                       trailing_newline=True)


def _report_rank_counts_for_taxon_or_exit(table_path, taxon, args, assembly_levels):
    """
    `--get-rank-counts` scoped to a taxon
    """
    prefixes = _source_prefixes(args.source)
    reps_only = args.refseq_reference_genomes_only
    scope_note = _counts_scope_note(args, assembly_levels)

    try:
        canonical, ranks_found_in = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        sys.exit(0)

    for rank in ranks_found_in:
        total = _count_at_rank(table_path, rank, canonical, prefixes=prefixes,
                               reps_only=reps_only, assembly_levels=assembly_levels)
        rows = rank_counts(table_path, "ncbi", scope_rank=rank, scope_taxon=canonical,
                           rep_filter=_rep_filter(reps_only),
                           accession_prefixes=prefixes,
                           assembly_levels=assembly_levels)
        report_message(f"The rank '{rank}' has {total:,} {canonical} entries{scope_note}.",
                       color=None, width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        print(render_rank_count_table(
            rows, count_header=f"Num. Unique Taxa under '{canonical}'"))

    report_message("Each count above is also how many genomes `--derep-rank <rank>` "
                   "would return, since dereplication keeps one genome per unique "
                   "taxon at that rank.", "yellow", width=100, initial_indent="    ",
                   subsequent_indent="    ", trailing_newline=True)


def _counts_scope_note(args, assembly_levels):
    """short human description of which filters the primary counts block reflects"""
    bits = []
    if args.source == "refseq":
        bits.append("in refseq")
    elif args.source == "genbank":
        bits.append("in genbank")
    if assembly_levels:
        levels = ", ".join(sorted(assembly_levels))
        bits.append(f"at assembly level {levels}")
    if not bits:
        return ""
    return " (" + ", ".join(bits) + ")"


def _report_taxon_counts_or_exit(table_path, taxon, args, assembly_levels):
    """
    Report how many genomes match `taxon` at each rank it occurs at, matching the GTDB
    helper's format

    A primary per-rank block for the base pool (scoped by --source and
    --assembly-level), then if --refseq-reference-genomes-only is set a separate
    "in considering only RefSeq reference genomes" block, like GTDB's reps block.

    The wording is explicit about WHICH filters each block reflects: the primary block
    reflects --source and --assembly-level (but not the reference-genome filter, which
    is applied only in the second block), so the two numbers aren't confused.

    Reporting per-rank (rather than one number for a single resolved rank) also means
    an ambiguous taxon name is informative here instead of an error.
    """
    prefixes = _source_prefixes(args.source)
    scope_note = _counts_scope_note(args, assembly_levels)

    try:
        canonical, ranks_found_in = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        sys.exit(0)

    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                               assembly_levels=assembly_levels)
        report_message(f"The rank '{rank}' has {count:,} {taxon} entries{scope_note}.",
                       color=None, width=100, initial_indent="    ",
                       subsequent_indent="    ", leading_newline=False,
                       trailing_newline=True)
        _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels)

    if args.refseq_reference_genomes_only:
        report_message("Of those, in considering only RefSeq reference genomes:",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", leading_newline=False,
                       trailing_newline=True)
        any_rep = False
        for rank in ranks_found_in:
            count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                                   reps_only=True, assembly_levels=assembly_levels)
            if count:
                any_rep = True
                report_message(f"The rank '{rank}' has {count:,} {taxon} RefSeq "
                               "reference genome entries.", color=None, width=100,
                               initial_indent="    ", subsequent_indent="    ",
                               leading_newline=False, trailing_newline=True)
                _report_derep_note(table_path, rank, taxon, args, prefixes,
                                   assembly_levels, reps_only=True)
        if not any_rep:
            report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank "
                           "as a RefSeq reference genome :(", "yellow", width=100,
                           initial_indent="    ", subsequent_indent="    ",
                           leading_newline=False, trailing_newline=True)
            sys.exit(0)


def _report_unassigned_domains(summary):
    """
    Say what an 'all' pull leaves behind, if anything
    """
    message = summary.message("NCBI") if summary else None
    if message:
        report_message(message, "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)


def _report_source_overlap(source):
    """
    `--source both` pulls most assemblies twice (a RefSeq GCF_ record and the GenBank
    GCA_ original it was derived from), so say so rather than let it look like twice
    as many genomes.
    """
    note = source_overlap_note(source)
    if note:
        report_message(note, "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)


def report_unique_taxa_counts_of_all_ranks(table_path, source="refseq", reps_only=False,
                                           assembly_levels=None, scoped_to_all=False):
    """
    Print, for each of the 7 ranks, how many unique taxa exist in the NCBI table
    (matching get-accs-from-gtdb's --get-rank-counts), scoped to `source`:
      refseq  -> GCF_ accessions only
      genbank -> GCA_ accessions only
      both    -> no source filter
    If reps_only, also print the counts among RefSeq reference genomes (which are
    RefSeq, so that sub-table is implicitly refseq-scoped).

    Counting is done by tax_counts, which EXCLUDES the "NA" placeholder.
    """
    prefixes = _source_prefixes(source)
    label = {"refseq": "RefSeq", "genbank": "GenBank", "both": "all"}.get(source, source)
    domain_assigned = True if scoped_to_all else None

    rows = rank_counts(table_path, "ncbi", accession_prefixes=prefixes,
                       assembly_levels=assembly_levels,
                       domain_assigned=domain_assigned)
    print("")
    print(render_rank_count_table(rows, count_header=f"Num. Unique Taxa ({label})"))
    print("")

    if scoped_to_all:
        report_message(scoped_counts_note("-t"), "yellow", width=100,
                       initial_indent="    ", subsequent_indent="    ",
                       trailing_newline=False)
        _report_unassigned_domains(unassigned_domain_summary(
            table_path, "ncbi", rep_filter=_rep_filter(reps_only),
            accession_prefixes=prefixes, assembly_levels=assembly_levels))

    if reps_only:
        rep_rows = rank_counts(table_path, "ncbi", accession_prefixes=prefixes,
                               rep_filter=_rep_filter(True),
                               assembly_levels=assembly_levels,
                               domain_assigned=domain_assigned)
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print(render_rank_count_table(rep_rows, count_header="Num. Unique Ref. Taxa"))
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


def _derep_is_on(args):
    """True when --derep-rank asks for actual dereplication."""
    return is_derep_on(getattr(args, "derep_rank", "off"))


def _select_all_dereplicated(table_path, args, assembly_levels=None):
    """`-t all` WITH --derep-rank: one selection per domain, merged."""
    try:
        selection = select_all_domains(
            table_path, "ncbi", derep_rank=args.derep_rank,
            reps_only=args.refseq_reference_genomes_only,
            accession_prefixes=_source_prefixes(args.source),
            assembly_levels=assembly_levels)
    except ValueError as e:
        report_message(str(e), "yellow")
        print("")
        sys.exit(0)

    report_message(f"Dereplicating within each domain "
                   f"({', '.join(selection.domains)}).", "yellow",
                   initial_indent="    ", subsequent_indent="    ",
                   trailing_newline=False)
    _report_unassigned_domains(getattr(selection, "unassigned", None))
    for w in selection.warnings:
        wprint(color_text(w, "yellow"))

    tab = (pa.Table.from_pylist(selection.rows) if selection.rows
           else pq.read_table(table_path, columns=_COLUMNS).slice(0, 0))
    label = ("all genomes (dereplicated within each domain: "
             + ", ".join(selection.domains) + ")")
    return tab, label


def _select_all(table_path, args, reps_only=False, assembly_levels=None):
    prefixes = _source_prefixes(getattr(args, "source", "refseq"))
    tab = read_pool(table_path, "ncbi", _COLUMNS,
                    rep_filter=_rep_filter(reps_only),
                    accession_prefixes=prefixes,
                    assembly_levels=assembly_levels,
                    domain_assigned=True)
    _report_unassigned_domains(unassigned_domain_summary(
        table_path, "ncbi", rep_filter=_rep_filter(reps_only),
        accession_prefixes=prefixes, assembly_levels=assembly_levels))
    return tab


def _select_by_taxid(table_path, taxid, assembly_levels=None):
    from bit.modules.taxonomy.tax_ranks import RANKS

    level_filter = ([("assembly_level", "in", set(assembly_levels))]
                    if assembly_levels else [])

    for rank in RANKS:
        tab = pq.read_table(table_path, columns=_COLUMNS,
                            filters=[(f"{rank}_taxid", "=", str(taxid))] + level_filter)
        if tab.num_rows:
            return tab

    return pq.read_table(table_path, columns=_COLUMNS,
                         filters=[("taxid", "=", str(taxid))] + level_filter)


def _apply_filters(tab, args):
    """
    Source scoping only
    """
    return _filter_by_source(tab, getattr(args, "source", "refseq"))


def _write_outputs(tab, args, label):
    taxon_raw = str(args.target_taxon)
    if taxon_raw.lower() == "all":
        taxon_raw = "all"
    taxon_for_filename = taxon_raw.replace(" ", "-").replace("/", "-").lower()

    suffix_bits = []
    if args.refseq_reference_genomes_only:
        suffix_bits.append("refseq-ref")
    elif getattr(args, "source", "refseq") != "both":
        suffix_bits.append(args.source.lower())
    suffix = ("-" + "-".join(suffix_bits)) if suffix_bits else ""

    acc_out = f"ncbi-{taxon_for_filename}{suffix}-accessions.txt"
    tab_out = f"ncbi-{taxon_for_filename}{suffix}-metadata.tsv"

    write_table_tsv(tab, tab_out)

    accs = tab.column("assembly_accession").to_pylist()
    write_accessions(acc_out, accs)

    print("")
    wprint(f"Wrote {len(accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out))
    print("")
