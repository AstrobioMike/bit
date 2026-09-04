import os
import sys
from collections import namedtuple
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
    CrossDomainTaxon,
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
from bit.modules.taxonomy.empty_selection import empty_pull_message
from bit.modules.taxonomy.get_accs_shared import (ASSEMBLY_LEVELS,
                                                  FILTERS_APPLIED_NOTE, PoolSpec,
                                                  apply_derep_default,
                                                  resolved_derep_rank,
                                                  derep_note as _shared_derep_note,
                                                  is_derep_on,
                                                  parse_assembly_levels as _parse_levels,
                                                  pull_count_lines,
                                                  scoped_counts_note,
                                                  source_prefixes as _shared_prefixes,
                                                  with_filters_note)
from bit.modules.taxonomy.exclusion_list import load_exclusion_cores


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

_NcbiSelection = namedtuple(
    "_NcbiSelection",
    ["table", "label", "rank", "taxon", "effective_derep_rank", "ref_selection"],
    defaults=[None, None],
)

_ASSEMBLY_LEVELS = ASSEMBLY_LEVELS


def parse_assembly_levels(value):
    """
    Kept as a name here; the body is shared with the GTDB helper's vocabulary
    """
    return _parse_levels(value)


def ncbi_table_path(force_update=False, quiet=True):
    get_ncbi_assembly_data(force_update=force_update, quiet=quiet)
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    return os.path.join(ncbi_dir, PARQUET_FILENAME)


def _source_prefixes(source):
    """
    Accession prefixes for an --ncbi-section value (None means no restriction)
    """
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

    # normally already done by the CLI's preflight, but this is also called directly,
    # and an unresolved --derep-rank sentinel would blow up further down
    apply_derep_default(args)

    exclude_cores = load_exclusion_cores(getattr(args, "exclusion_list", None))

    table_path = ncbi_table_path()
    _report_ncbi_date(table_path)

    if getattr(args, "get_table", False):
        copy_ncbi_table(table_path)
        sys.exit(0)

    assembly_levels = parse_assembly_levels(getattr(args, "assembly_level", None))

    target = str(args.target_taxon) if args.target_taxon else ""
    named_taxon = bool(target) and not is_all_target(target) and not target.isdigit()

    if getattr(args, "get_rank_counts", False):
        if named_taxon:
            _report_rank_counts_for_taxon_or_exit(table_path, target, args,
                                                  assembly_levels,
                                                  exclude_cores=exclude_cores)
        else:
            report_unique_taxa_counts_of_all_ranks(
                table_path, source=args.ncbi_section,
                reps_only=args.refseq_reference_genomes_only,
                assembly_levels=assembly_levels,
                scoped_to_all=is_all_target(target),
                                            exclude_cores=exclude_cores)
        sys.exit(0)

    if getattr(args, "get_taxon_counts", False) and named_taxon:
        _report_taxon_counts_or_exit(table_path, target, args, assembly_levels,
                                     exclude_cores=exclude_cores)
        sys.exit(0)

    selection = _select_rows(table_path, args, assembly_levels)
    tab, label = selection.table, selection.label

    if args.get_taxon_counts:
        report_message(with_filters_note(
                           f"There are {tab.num_rows:,} genome(s) under {label}."),
                       "none",
                       initial_indent="    ", subsequent_indent="    ",
                       trailing_newline=True)
        sys.exit(0)

    if tab.num_rows == 0:
        report_message(empty_pull_message(
                           f"No genomes were found under {label}.",
                           selection.ref_selection,
                           assembly_levels=assembly_levels,
                           ncbi_section=getattr(args, "ncbi_section", None),
                           reps_only_requested=bool(
                               args.refseq_reference_genomes_only),
                           emoticon=":("),
                       "none",
                       initial_indent="    ", subsequent_indent="    ")
        print("")
        sys.exit(0)

    _report_pull_counts(table_path, args, assembly_levels, selection, tab.num_rows,
                        exclude_cores=exclude_cores)
    _write_outputs(tab, args, selection)


def _report_pull_counts(table_path, args, assembly_levels, selection, kept_n,
                        exclude_cores=None):
    """
    The "The rank 'X' has N <taxon> entries." (+ dereplicated) block a taxon-name
    pull prints just above its "Wrote N accession(s)" lines.

    Only the taxon-name path has a single resolved rank to report against; the `all`
    and taxid paths pass rank=None and are left alone.
    """
    rank, taxon = selection.rank, selection.taxon
    if rank is None:
        return

    pool = PoolSpec(
        table_path, "ncbi",
        rep_filter=_rep_filter(args.refseq_reference_genomes_only),
        accession_prefixes=_source_prefixes(getattr(args, "ncbi_section", "refseq")),
        assembly_levels=assembly_levels or None, label="NCBI", taxon_flag="-t",
        exclude_cores=exclude_cores)
    header, derep = pull_count_lines(pool, rank, taxon,
                                    selection.effective_derep_rank, kept_n)
    print("")
    wprint("  " + header)
    if derep:
        wprint("    " + derep)


def _prefix_mask(acc_col, prefixes):
    mask = None
    for p in prefixes:
        m = pc.starts_with(acc_col, p)
        mask = m if mask is None else pc.or_(mask, m)
    return mask


def _count_at_rank(table_path, rank, taxon, prefixes=None, reps_only=False,
                   assembly_levels=None, exclude_cores=None):
    """
    Count rows where column `rank` == `taxon`, with the set POOL filters applied
    (source prefix, RefSeq-reference-only, assembly level)
    """
    return count_genomes(table_path, "ncbi", rank=rank, taxon=taxon,
                         rep_filter=_rep_filter(reps_only),
                         accession_prefixes=prefixes,
                         assembly_levels=assembly_levels,
                         exclude_cores=exclude_cores)


def _rep_filter(reps_only):
    """
    The RefSeq-reference-genome predicate, or None
    """
    return representatives_filter("ncbi", "refseq" if reps_only else None)


def _derep_count_at_rank(table_path, rank, taxon, derep_rank, prefixes=None,
                         reps_only=False, assembly_levels=None, exclude_cores=None):
    """
    How many genomes survive dereplication at `derep_rank`, under the same pool
    """
    return derep_size(table_path, "ncbi", rank, taxon, derep_rank,
                      rep_filter=_rep_filter(reps_only),
                      accession_prefixes=prefixes,
                      assembly_levels=assembly_levels,
                      exclude_cores=exclude_cores)


def _derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                reps_only=False, exclude_cores=None):
    """
    The "...dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off

    Returns (line_or_None, warnings).
    """
    pool = PoolSpec(table_path, "ncbi", rep_filter=_rep_filter(reps_only),
                    accession_prefixes=prefixes, assembly_levels=assembly_levels,
                    label="NCBI", taxon_flag="-t",
                    exclude_cores=exclude_cores)
    return _shared_derep_note(pool, rank, taxon,
                              resolved_derep_rank(args))


def _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                       reps_only=False, exclude_cores=None):
    """
    Print the dereplicated-count line (and any 'auto' advisory) for one rank
    """
    line, warnings = _derep_note(table_path, rank, taxon, args, prefixes,
                                 assembly_levels, reps_only=reps_only,
                                 exclude_cores=exclude_cores)
    if line:
        report_message(line, color=None, width=100, initial_indent="      ",
                       subsequent_indent="      ", leading_newline=False,
                       trailing_newline=True)
    for warning in warnings:
        report_message(warning, "yellow", width=100, initial_indent="      ",
                       subsequent_indent="      ", leading_newline=False,
                       trailing_newline=True)


def _report_rank_counts_for_taxon_or_exit(table_path, taxon, args, assembly_levels, exclude_cores=None):
    """
    `--get-rank-counts` scoped to a taxon
    """
    prefixes = _source_prefixes(args.ncbi_section)
    reps_only = args.refseq_reference_genomes_only

    try:
        canonical, ranks_found_in, _domains_found = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        sys.exit(0)

    for rank in ranks_found_in:
        total = _count_at_rank(table_path, rank, canonical, prefixes=prefixes,
                               reps_only=reps_only, assembly_levels=assembly_levels,
                               exclude_cores=exclude_cores)
        rows = rank_counts(table_path, "ncbi", scope_rank=rank, scope_taxon=canonical,
                           rep_filter=_rep_filter(reps_only),
                           accession_prefixes=prefixes,
                           assembly_levels=assembly_levels,
                           exclude_cores=exclude_cores)
        report_message(with_filters_note(
                           f"The rank '{rank}' has {total:,} {canonical} entries."),
                       color=None, width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        print(render_rank_count_table(
            rows, count_header=f"Num. Unique Taxa under '{canonical}'"))

    report_message("Each count above is also how many genomes `--derep-rank <rank>` "
                   "would return, since dereplication keeps one genome per unique "
                   "taxon at that rank.", "yellow", width=90, initial_indent="  ",
                   subsequent_indent="  ", trailing_newline=True)


def _report_taxon_counts_or_exit(table_path, taxon, args, assembly_levels, exclude_cores=None):
    """
    Report how many genomes match `taxon` at each rank it occurs at, matching the GTDB
    helper's format

    A primary per-rank block for the base pool (scoped by every pool filter --
    --ncbi-section, --assembly-level, --exclusion-list), then if
    --refseq-reference-genomes-only is set a separate "in considering only RefSeq
    reference genomes" block, like GTDB's reps block.

    Each primary line carries the generic "(after any specified filters)" tag rather
    than enumerating which filters are set, so it can't drift from what the pool
    actually applies. The reference-genome filter is called out separately because it
    is applied only in the second block.

    Reporting per-rank (rather than one number for a single resolved rank) also means
    an ambiguous taxon name is informative here instead of an error.
    """
    prefixes = _source_prefixes(args.ncbi_section)

    try:
        canonical, ranks_found_in, _domains_found = _resolve_ranks(table_path, taxon)
    except TaxonNotFound:
        report_message(f"Input taxon '{taxon}' doesn't seem to exist at any rank :(",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", trailing_newline=True)
        sys.exit(0)

    taxon = canonical

    print("")
    for rank in ranks_found_in:
        count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                               assembly_levels=assembly_levels,
                               exclude_cores=exclude_cores)
        report_message(with_filters_note(
                           f"The rank '{rank}' has {count:,} {taxon} entries."),
                       color=None, width=100, initial_indent="    ",
                       subsequent_indent="    ", leading_newline=False,
                       trailing_newline=False)
        _report_derep_note(table_path, rank, taxon, args, prefixes, assembly_levels,
                           exclude_cores=exclude_cores)

    if args.refseq_reference_genomes_only:
        report_message("Of those, in considering only RefSeq reference genomes:",
                       "yellow", width=100, initial_indent="    ",
                       subsequent_indent="    ", leading_newline=False,
                       trailing_newline=True)
        any_rep = False
        for rank in ranks_found_in:
            count = _count_at_rank(table_path, rank, taxon, prefixes=prefixes,
                                   reps_only=True, assembly_levels=assembly_levels,
                                   exclude_cores=exclude_cores)
            if count:
                any_rep = True
                report_message(f"The rank '{rank}' has {count:,} {taxon} RefSeq "
                               "reference genome entries.", color=None, width=100,
                               initial_indent="    ", subsequent_indent="    ",
                               leading_newline=False, trailing_newline=True)
                _report_derep_note(table_path, rank, taxon, args, prefixes,
                                   assembly_levels, reps_only=True,
                           exclude_cores=exclude_cores)
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


def report_unique_taxa_counts_of_all_ranks(table_path, source="refseq", reps_only=False,
                                           assembly_levels=None, scoped_to_all=False, exclude_cores=None):
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
                       domain_assigned=domain_assigned,
                       exclude_cores=exclude_cores)
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
                               domain_assigned=domain_assigned,
                               exclude_cores=exclude_cores)
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
    print("\n    Using NCBI assembly-data retrieved: " + date_str)


def _derep_is_on(args):
    """
    True when --derep-rank asks for actual dereplication.
    """
    return is_derep_on(resolved_derep_rank(args))


def _select_all_dereplicated(table_path, args, assembly_levels=None):
    """
    `-t all` WITH --derep-rank: one selection per domain, merged.
    """
    try:
        selection = select_all_domains(
            table_path, "ncbi", derep_rank=args.derep_rank,
            reps_only=args.refseq_reference_genomes_only,
            accession_prefixes=_source_prefixes(args.ncbi_section),
            assembly_levels=assembly_levels,
            exclude_cores=load_exclusion_cores(
                getattr(args, "exclusion_list", None)))
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
    prefixes = _source_prefixes(getattr(args, "ncbi_section", "refseq"))
    tab = read_pool(table_path, "ncbi", _COLUMNS,
                    rep_filter=_rep_filter(reps_only),
                    accession_prefixes=prefixes,
                    assembly_levels=assembly_levels,
                    domain_assigned=True)
    _report_unassigned_domains(unassigned_domain_summary(
        table_path, "ncbi", rep_filter=_rep_filter(reps_only),
        accession_prefixes=prefixes, assembly_levels=assembly_levels))
    return tab


def _select_rows(table_path, args, assembly_levels=None):
    """
    Resolve the target and return an _NcbiSelection(table, label, rank, taxon) where
    `rank` is the resolved rank for a taxon-name search.

    Mirrors GToTree's function of the same name. Every branch ends in one
    _NcbiSelection so the caller has a single shape to work with, and an empty result
    is NOT reported here, it goes back to the one check in
    get_accessions_from_ncbi(), so all three paths produce one message built one way.
    """
    if is_all_target(args.target_taxon):
        if _derep_is_on(args):
            tab, label = _select_all_dereplicated(table_path, args, assembly_levels)
        else:
            tab = _select_all(table_path, args,
                              reps_only=args.refseq_reference_genomes_only,
                              assembly_levels=assembly_levels)
            label = "all genomes"
        return _NcbiSelection(_apply_filters(tab, args), label, None, "all")

    if str(args.target_taxon).isdigit():
        tab = _select_by_taxid(table_path, str(args.target_taxon),
                               assembly_levels=assembly_levels)
        label = f"taxid {args.target_taxon}"
        return _NcbiSelection(_apply_filters(tab, args), label, None,
                              f"taxid-{args.target_taxon}")

    try:
        canonical, rank, domain = resolve_taxon(
            table_path, args.target_taxon,
            rank=getattr(args, "target_rank", None),
            domain=getattr(args, "target_domain", None))
    except AmbiguousTaxon as e:
        report_message(
            f"'{e.taxon}' occurs at more than one rank "
            f"({', '.join(e.ranks_found)}). Specify which with `-r`, or pass "
            f"the NCBI taxid to `-t` instead.", "yellow")
        print("")
        sys.exit(0)
    except CrossDomainTaxon as e:
        report_message(
            f"'{e.taxon}' occurs in more than one domain "
            f"({', '.join(e.domains_found)}). Specify which with "
            f"`--target-domain`.", "yellow")
        print("")
        sys.exit(0)
    except TaxonNotFound:
        report_message(f"Input taxon '{args.target_taxon}' doesn't seem to "
                       f"exist at any rank :(", "yellow")
        print("")
        sys.exit(0)

    try:
        selection = select_ref_genomes(
            table_path, "ncbi", canonical, target_rank=rank,
            diagnose_empty=True,
            derep_rank=resolved_derep_rank(args),
            reps_only=args.refseq_reference_genomes_only,
            accession_prefixes=_source_prefixes(args.ncbi_section),
            assembly_levels=assembly_levels,
            exclude_cores=load_exclusion_cores(
                getattr(args, "exclusion_list", None)),
            target_domain=domain)
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

    return _NcbiSelection(_apply_filters(tab, args), label, rank, canonical,
                          selection.effective_derep_rank,
                          ref_selection=selection)


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
    return _filter_by_source(tab, getattr(args, "ncbi_section", "refseq"))


def _write_outputs(tab, args, selection):
    """
    Names match GToTree's NCBI helper and both GTDB helpers:
    ncbi-<taxon>-<rank><suffix>-accs.txt

    The taxon comes from the SELECTION, not the raw `-t` value, so the same organism
    lands in the same filename however it was typed or aliased. The rank is included
    because a name can exist at more than one, and is omitted for `all` and taxid
    pulls, which resolve to no single rank.
    """
    taxon_for_filename = selection.taxon.replace(" ", "-").replace("/", "-").lower()
    rank_bit = f"-{selection.rank}" if selection.rank else ""

    suffix_bits = []
    if args.refseq_reference_genomes_only:
        suffix_bits.append("refseq-ref")
    elif getattr(args, "ncbi_section", "refseq") != "both":
        suffix_bits.append(args.ncbi_section.lower())
    suffix = ("-" + "-".join(suffix_bits)) if suffix_bits else ""

    acc_out = f"ncbi-{taxon_for_filename}{rank_bit}{suffix}-accs.txt"
    tab_out = f"ncbi-{taxon_for_filename}{rank_bit}{suffix}-metadata.tsv"

    write_table_tsv(tab, tab_out)

    accs = tab.column("assembly_accession").to_pylist()
    write_accessions(acc_out, accs)

    print("")
    wprint("  " + f"Wrote {len(accs):,} accession(s) to:")
    wprint("    " + color_text(acc_out))
    print("")
    wprint("  " + "Associated taxonomy and metadata of these targets written to:")
    wprint("    " + color_text(tab_out))
    print("")
