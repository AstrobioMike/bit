import sys
import os
import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore
from bit.modules.general import (color_text, report_message, wprint, write_table_tsv,
                                 write_accessions)
from bit.modules.gtdb.get_gtdb_data import (get_gtdb_data,
                                            report_gtdb_version_info as _read_gtdb_version_info)
from bit.modules.taxonomy.tax_ranks import RANKS
from bit.modules.taxonomy.tax_select import (resolve_taxon, select,
                                             TaxonNotFound, AmbiguousTaxon,
                                             find_ranks_for_taxon as _resolve_ranks)
from bit.modules.taxonomy.tax_derep import (select_ref_genomes, select_all_domains,
                                            resolve_derep_rank)
from bit.modules.taxonomy.tax_counts import (representatives_filter, count_genomes,
                                             derep_size, rank_counts,
                                             render_rank_count_table)
from bit.modules.taxonomy.tax_targets import (domains_in_asset, is_all_target,
                                              unassigned_domain_summary)
from bit.modules.taxonomy.get_accs_shared import (PoolSpec, all_derep_size,
                                                  derep_note as _shared_derep_note,
                                                  is_derep_on, scoped_counts_note)


_RANK_COLUMNS = list(RANKS)


def _all_columns(gtdb_path):
    return pq.ParquetFile(gtdb_path).schema_arrow.names


def get_accessions_from_gtdb(args):

    gtdb_path = get_gtdb_data(quiet=True)

    if args.get_table:
        copy_gtdb_table(gtdb_path)
        sys.exit(0)

    if args.get_taxon_counts and not args.target_taxon:
        print("")
        wprint(color_text("A specific taxon needs to also be provided to the `-t` flag in order to use `--get-taxon-counts`.", "yellow"))
        print("")
        wprint("  E.g.: bit get-accs-from-gtdb --get-taxon-counts -t Alteromonas")
        print("")
        sys.exit(0)

    if args.gtdb_representatives_only and args.refseq_reference_genomes_only:
        print("")
        wprint(color_text("Only one of `--gtdb-representatives-only` or `--refseq-reference-genomes-only` can be provided.", "yellow"))
        print("")
        sys.exit(0)

    _report_gtdb_version(gtdb_path)

    if args.gtdb_representatives_only:
        representatives_source = "gtdb"
    elif args.refseq_reference_genomes_only:
        representatives_source = "refseq"
    else:
        representatives_source = None

    named_taxon = bool(args.target_taxon) and not is_all_target(args.target_taxon)

    if args.get_rank_counts:
        if named_taxon:
            report_rank_counts_for_taxon(gtdb_path, args.target_taxon,
                                         representatives_source)
        else:
            report_unique_taxa_counts_of_all_ranks(
                gtdb_path, representatives_source=representatives_source,
                scoped_to_all=is_all_target(args.target_taxon))
        sys.exit(0)

    if not args.target_taxon:
        return

    target_taxon, resolved_rank = _resolve_or_exit(gtdb_path, args.target_taxon, args.target_rank)

    if args.get_taxon_counts:
        report_taxon_counts(gtdb_path, target_taxon, representatives_source, args)
        sys.exit(0)

    if is_all_target(target_taxon):
        if _derep_is_on(args):
            _write_all_dereplicated(gtdb_path, args, representatives_source)
        else:
            _write_all(gtdb_path, representatives_source)
        sys.exit(0)

    table, kept_rank = _select_rows(gtdb_path, args, target_taxon, resolved_rank,
                                    representatives_source)
    _write_outputs(table, target_taxon, kept_rank, representatives_source)
    sys.exit(0)


def _select_rows(gtdb_path, args, target_taxon, resolved_rank, representatives_source):
    """
    The taxon slice to write, as a pyarrow Table carrying every column of the asset
    """
    reps_only = representatives_source is not None
    cols = _all_columns(gtdb_path)

    table = select(gtdb_path, "gtdb", resolved_rank, target_taxon,
                   reps_only=reps_only, columns=cols)

    derep_rank = getattr(args, "derep_rank", "off")
    if derep_rank in (None, "off", "none", "None"):
        return table, resolved_rank

    try:
        selection = select_ref_genomes(
            gtdb_path, "gtdb", target_taxon, target_rank=resolved_rank,
            derep_rank=derep_rank, reps_only=reps_only)
    except ValueError as e:
        print("")
        wprint(color_text(str(e), "yellow"))
        print("")
        sys.exit(0)

    for w in selection.warnings:
        wprint(color_text(w, "yellow"))
    if selection.warnings:
        print("")

    table = _restrict_to_accessions(table, selection.accessions)

    if selection.effective_derep_rank:
        wprint(f"Dereplicated to one genome per {selection.effective_derep_rank}: "
               f"{table.num_rows:,} genome(s) kept.")
        print("")

    return table, resolved_rank


def _restrict_to_accessions(table, accessions):
    """
    Keep only the rows whose accession is in `accessions`
    """
    wanted = pa.array(sorted(set(accessions)), type=pa.string())
    col = table.column("ncbi_genbank_assembly_accession")
    return table.filter(pc.is_in(col, value_set=wanted))


def _write_outputs(table, taxon, rank, representatives_source=None):
    """Write the accession list + metadata TSV for a resolved taxon."""
    taxon_for_filename = taxon.replace(" ", "-").replace("/", "-").lower()

    if representatives_source:
        stem = f"gtdb-{taxon_for_filename}-{rank}-{representatives_source}-rep"
    else:
        stem = f"gtdb-{taxon_for_filename}-{rank}"

    acc_out_filename = stem + "-accs.txt"
    tab_out_filename = stem + "-metadata.tsv"

    write_table_tsv(table, tab_out_filename)

    accs = table.column("ncbi_genbank_assembly_accession").to_pylist()
    write_accessions(acc_out_filename, accs)

    print("")
    wprint(f"Wrote {len(accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated taxonomy and metadata of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def _derep_is_on(args):
    """True when --derep-rank asks for actual dereplication."""
    return is_derep_on(getattr(args, "derep_rank", "off"))


def _write_all_dereplicated(gtdb_path, args, representatives_source=None):
    """`-t all` WITH --derep-rank: one selection per domain, merged."""
    reps_only = representatives_source is not None

    try:
        selection = select_all_domains(gtdb_path, "gtdb", derep_rank=args.derep_rank,
                                       reps_only=reps_only)
    except ValueError as e:
        print("")
        wprint(color_text(str(e), "yellow"))
        print("")
        sys.exit(0)

    print("")
    wprint(color_text(f"Dereplicating within each domain "
                      f"({', '.join(selection.domains)}).", "yellow"))
    _report_unassigned_domains(getattr(selection, "unassigned", None))
    for w in selection.warnings:
        wprint(color_text(w, "yellow"))
    print("")

    stem = (f"gtdb-arc-and-bac-{representatives_source}-rep"
            if representatives_source else "gtdb-arc-and-bac")
    acc_out_filename = stem + "-accessions.txt"
    tab_out_filename = stem + "-metadata.tsv"

    table = _restrict_to_accessions(
        pq.read_table(gtdb_path, columns=_all_columns(gtdb_path)),
        selection.accessions)
    write_table_tsv(table, tab_out_filename)
    write_accessions(acc_out_filename, selection.accessions)

    print("")
    wprint(f"Wrote {len(selection.accessions):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated taxonomy and metadata written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def _report_unassigned_domains(summary):
    """
    Say what an 'all' pull leaves behind, if anything
    """
    message = summary.message("GTDB") if summary else None
    if message:
        wprint(color_text(message, "yellow"))
        print("")


def _write_all(gtdb_path, representatives_source=None):
    """
    Bulk dump of every accession (optionally representatives-only)
    """
    rep_filter = _rep_filter_for(representatives_source)
    filters = [(rep_filter[0], "=", rep_filter[1])] if rep_filter else None

    table = pq.read_table(gtdb_path, columns=_all_columns(gtdb_path), filters=filters)

    if representatives_source:
        stem = f"gtdb-arc-and-bac-{representatives_source}-rep"
        tab_out_filename = stem + "-metadata.tsv"
        write_table_tsv(table, tab_out_filename)
    else:
        stem = "gtdb-arc-and-bac"
        tab_out_filename = None

    acc_out_filename = stem + "-accessions.txt"
    accs = table.column("ncbi_genbank_assembly_accession").to_pylist()
    write_accessions(acc_out_filename, accs)

    print("")
    wprint(f"Wrote {len(accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    if tab_out_filename:
        wprint("Associated taxonomy and metadata written to:")
        wprint("  " + color_text(tab_out_filename))
        print("")
    else:
        wprint("The full GTDB table holds all of these; `--get-table` writes it out as "
               + color_text("gtdb-arc-and-bac-metadata.tsv"))
        print("")


def _resolve_or_exit(gtdb_path, taxon, rank=None):
    if is_all_target(taxon):
        return "all", None
    try:
        canonical, resolved_rank = resolve_taxon(gtdb_path, taxon, rank)
    except AmbiguousTaxon:
        wprint(color_text("Since '" + str(taxon) + "' occurs at more than 1 rank, we'll need to specify "
               "which rank is wanted as well before we pull the accessions. This can be specified with the `-r` flag.", "yellow"))
        print("")
        sys.exit(0)
    except TaxonNotFound:
        wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
        print("")
        sys.exit(0)
    # if canonical != taxon:
    #     wprint(color_text("Matched input '" + taxon + "' to GTDB taxon '" + canonical + "'.", "yellow"))
    #     print("")
    return canonical, resolved_rank


def _report_gtdb_version(gtdb_path):
    report_gtdb_version_info(os.path.dirname(gtdb_path))


def report_gtdb_version_info(location):
    gtdb_version, gtdb_release_date = _read_gtdb_version_info(location)
    print("\n    Using GTDB " + gtdb_version + ": " + gtdb_release_date)


def copy_gtdb_table(gtdb_path):
    """
    Write the parquet object as a tsv
    """
    report_gtdb_version_info(os.path.dirname(gtdb_path))

    out_name = "gtdb-arc-and-bac-metadata.tsv"
    write_table_tsv(pq.read_table(gtdb_path), out_name)

    print("")
    wprint("GTDB table written to:")
    print(color_text("    " + out_name + "\n"))


def _rep_filter_for(representatives_source):
    """
    The Parquet predicate for a representatives source (or None for no filter).
    Delegates to tax_counts so the column names live with the SourceSpec.
    """
    if representatives_source == "gtdb":
        return representatives_filter("gtdb", "source")
    if representatives_source == "refseq":
        return representatives_filter("gtdb", "refseq")
    return None


def _derep_note(gtdb_path, rank, taxon, args, rep_filter=None):
    """
    The "...dereplicated at X, that would be N" line for one rank, or None when
    dereplication is off
    """
    pool = PoolSpec(gtdb_path, "gtdb", rep_filter=rep_filter, label="GTDB",
                    taxon_flag="-t")
    return _shared_derep_note(pool, rank, taxon, getattr(args, "derep_rank", "off"))


def _report_derep_note(gtdb_path, rank, taxon, args, rep_filter=None):
    """Print the dereplicated-count line (and any 'auto' advisory) for one rank."""
    if args is None:
        return
    line, warnings = _derep_note(gtdb_path, rank, taxon, args, rep_filter=rep_filter)
    if line:
        wprint("    " + line)
    for warning in warnings:
        wprint(color_text("      " + warning, "yellow"))


def report_taxon_counts(gtdb_path, taxon, representatives_source=None, args=None):
    """
    Counts for a specific taxon (or 'all'), read straight from the Parquet.

    When --derep-rank is set, each per-rank line is followed by how many genomes would
    survive dereplication
    """
    rep_filter = _rep_filter_for(representatives_source)

    if taxon == "all":
        count = count_genomes(gtdb_path, "gtdb")
        print("")
        wprint(f"  There are {count:,} total genomes in the database.")
        if args is not None and _derep_is_on(args):
            wprint(f"    Dereplicated within each domain, that would be "
                   f"{_all_derep_size(gtdb_path, 'gtdb', args):,} genome(s).")
        print("")
        if representatives_source:
            rep_source = "GTDB" if representatives_source == "gtdb" else "RefSeq"
            wprint(color_text(f"In considering only {rep_source} representative genomes:",
                              "yellow"))
            print("")
            wprint(f"  There are {count_genomes(gtdb_path, 'gtdb', rep_filter=rep_filter):,} "
                   f"total representative genomes in the database.")
            print("")
        return

    try:
        canonical, ranks_found_in = _resolve_ranks(gtdb_path, taxon)
    except TaxonNotFound:
        wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
        print("")
        sys.exit(0)

    taxon = canonical

    for rank in ranks_found_in:
        count = count_genomes(gtdb_path, "gtdb", rank=rank, taxon=taxon)
        wprint("  The rank '" + rank + "' has " + f"{count:,}" + " " + taxon + " entries.")
        _report_derep_note(gtdb_path, rank, taxon, args)

    print("")

    if representatives_source:
        rep_source = "GTDB" if representatives_source == "gtdb" else "RefSeq"
        wprint(color_text(f"In considering only {rep_source} representative genomes:", "yellow"))
        print("")
        any_rep = False
        for rank in ranks_found_in:
            count = count_genomes(gtdb_path, "gtdb", rank=rank, taxon=taxon,
                                  rep_filter=rep_filter)
            if count:
                any_rep = True
                wprint("  The rank '" + rank + "' has " + f"{count:,}" + " " + taxon +
                       " representative genome entries.")
                _report_derep_note(gtdb_path, rank, taxon, args, rep_filter=rep_filter)
                print("")
        if not any_rep:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any "
                              "rank as a representative genome :(", "yellow"))
            print("")
            sys.exit(0)


def _all_derep_size(path, source, args, rep_filter=None):
    """How many genomes `-t all --derep-rank X` would return."""
    pool = PoolSpec(path, source, rep_filter=rep_filter, label="GTDB", taxon_flag="-t")
    return all_derep_size(pool, args.derep_rank)


def report_rank_counts_for_taxon(gtdb_path, taxon, representatives_source=None):
    """
    `--get-rank-counts` scoped to a taxon
    """
    rep_filter = _rep_filter_for(representatives_source)

    try:
        canonical, ranks_found_in = _resolve_ranks(gtdb_path, taxon)
    except TaxonNotFound:
        wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
        print("")
        sys.exit(0)

    for rank in ranks_found_in:
        total = count_genomes(gtdb_path, "gtdb", rank=rank, taxon=canonical,
                              rep_filter=rep_filter)
        rows = rank_counts(gtdb_path, "gtdb", scope_rank=rank, scope_taxon=canonical,
                           rep_filter=rep_filter)
        print("")
        wprint(f"  The rank '{rank}' has {total:,} {canonical} entries.")
        print("")
        print(render_rank_count_table(
            rows, count_header=f"Num. Unique Taxa under '{canonical}'"))

    report_message("Each count above is also how many genomes `--derep-rank <rank>` "
                   "would return, since dereplication keeps one genome per unique "
                   "taxon at that rank.", "yellow", width=90, initial_indent="  ",
                   subsequent_indent="  ", trailing_newline=True)


def report_unique_taxa_counts_of_all_ranks(gtdb_path, representatives_source=None,
                                           scoped_to_all=False):
    """
    Counts of unique taxa at each rank across the whole table.

    Counting is done by tax_counts, which excludes the "NA" placeholder.

    scoped_to_all restricts the counts to genomes with an assigned domain, i.e., what
    `-t all` can actually pull, so each count is also how many accessions
    `-t all --derep-rank <that rank>` returns.
    """
    domain_assigned = True if scoped_to_all else None

    print("")
    print(render_rank_count_table(rank_counts(gtdb_path, "gtdb",
                                              domain_assigned=domain_assigned)))
    print("")

    if scoped_to_all:
        report_message(scoped_counts_note("-t"), "yellow", width=100,
                       initial_indent="    ", subsequent_indent="    ",
                       trailing_newline=True)
        _report_unassigned_domains(unassigned_domain_summary(gtdb_path, "gtdb"))

    if representatives_source == "gtdb":
        report_message("(The `--gtdb-representatives-only` flag doesn't change these counts: "
                       "every GTDB taxon has a representative genome, so the number of unique "
                       "taxa per rank is the same with or without it.)",
                       "yellow", initial_indent="    ", subsequent_indent="    ",
                        width=100, leading_newline=False, trailing_newline=True)
    elif representatives_source == "refseq":
        rep_rows = rank_counts(gtdb_path, "gtdb",
                               rep_filter=_rep_filter_for("refseq"))
        wprint(color_text("  In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print(render_rank_count_table(rep_rows, count_header="Num. Unique Ref. Taxa"))
        print("")
