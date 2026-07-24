import sys
import os
import pyarrow.parquet as pq # type: ignore
from bit.modules.general import color_text, wprint
from bit.modules.gtdb.get_gtdb_data import (get_gtdb_data,
                                            report_gtdb_version_info as _read_gtdb_version_info)
from bit.modules.taxonomy.tax_ranks import RANKS
from bit.modules.taxonomy.tax_select import (resolve_taxon, select,
                                             TaxonNotFound, AmbiguousTaxon)
from bit.modules.taxonomy.tax_derep import select_ref_genomes


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

    if args.get_rank_counts:
        full = _read_rank_columns(gtdb_path)
        rep = _apply_reps_filter(gtdb_path, representatives_source)
        get_unique_taxa_counts_of_all_ranks(full, rep, representatives_source=representatives_source)
        sys.exit(0)

    if not args.target_taxon:
        return

    target_taxon, resolved_rank = _resolve_or_exit(gtdb_path, args.target_taxon, args.target_rank)

    if args.get_taxon_counts:
        if target_taxon == "all":
            full = _read_rank_columns(gtdb_path)
            rep = _apply_reps_filter(gtdb_path, representatives_source)
        else:
            full = _read_rank_columns(gtdb_path)
            rep = _apply_reps_filter(gtdb_path, representatives_source) if representatives_source else None
        get_unique_taxon_counts(target_taxon, full, rep, representatives_source=representatives_source)
        sys.exit(0)

    working, rank = _read_taxon_full_slice(gtdb_path, target_taxon, resolved_rank,
                                           representatives_source)

    working = _apply_derep(gtdb_path, args, target_taxon, rank, working,
                           representatives_source)

    _write_accessions(target_taxon, working, rank,
                      representatives_source=representatives_source)
    sys.exit(0)


def _apply_derep(gtdb_path, args, target_taxon, rank, working, representatives_source):
    """
    Dereplicate `working` (a full-slice DataFrame) down to one best genome per unique
    value of --derep-rank, using the shared selection core so this matches
    get-accs-from-ncbi and GToTree.

    Derep only applies to a resolved taxon name, 'all' has no single rank to
    dereplicate relative to, so it is passed through untouched. Mike, remember, this is
    not gen-mg's dereplication (see modules/gen_mg/selection.py); that one samples
    randomly toward a target count and is deliberately a different approach.
    """
    derep_rank = getattr(args, "derep_rank", "off")
    if derep_rank in (None, "off") or target_taxon == "all":
        return working

    try:
        selection = select_ref_genomes(
            gtdb_path, "gtdb", target_taxon, target_rank=rank,
            derep_rank=derep_rank,
            reps_only=representatives_source is not None)
    except ValueError as e:
        print("")
        wprint(color_text(str(e), "yellow"))
        print("")
        sys.exit(0)

    for w in selection.warnings:
        wprint(color_text(w, "yellow"))
    if selection.warnings:
        print("")

    keep = set(selection.accessions)
    out = working[working["ncbi_genbank_assembly_accession"].isin(keep)]

    if selection.effective_derep_rank:
        wprint(f"Dereplicated to one genome per {selection.effective_derep_rank}: "
               f"{len(out):,} genome(s) kept.")
        print("")

    return out


def _read_rank_columns(gtdb_path):
    return pq.read_table(gtdb_path, columns=_RANK_COLUMNS).to_pandas()


def _apply_reps_filter(gtdb_path, representatives_source):
    if not representatives_source:
        return None
    if representatives_source == "gtdb":
        filt = [("gtdb_representative", "=", "t")]
    else:
        filt = [("ncbi_refseq_category", "=", "reference genome")]
    return pq.read_table(gtdb_path, columns=_RANK_COLUMNS, filters=filt).to_pandas()


def _read_taxon_full_slice(gtdb_path, taxon, resolved_rank, representatives_source):
    reps_only = representatives_source is not None
    cols = _all_columns(gtdb_path)

    if taxon == "all":
        if reps_only:
            filt = [("gtdb_representative", "=", "t")] if representatives_source == "gtdb" \
                else [("ncbi_refseq_category", "=", "reference genome")]
            tab = pq.read_table(gtdb_path, columns=cols, filters=filt).to_pandas()
        else:
            tab = pq.read_table(gtdb_path, columns=cols).to_pandas()
        return tab, None

    tbl = select(gtdb_path, "gtdb", resolved_rank, taxon,
                 reps_only=reps_only, columns=cols)
    return tbl.to_pandas(), resolved_rank


def _resolve_or_exit(gtdb_path, taxon, rank=None):
    if taxon.lower() == "all":
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
    if canonical != taxon:
        wprint(color_text("Matched input '" + taxon + "' to GTDB taxon '" + canonical + "'.", "yellow"))
        print("")
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
    df = pq.read_table(gtdb_path).to_pandas()
    df.to_csv(out_name, sep="\t", index=False)

    print("")
    wprint("GTDB table written to:")
    print(color_text("    " + out_name + "\n"))


def find_ranks_for_taxon(taxon, tab):
    """ ranks (of the 7) whose column contains `taxon`, within a DataFrame """
    return [rank for rank in RANKS if taxon in tab[rank].unique()]


def get_accessions(taxon, gtdb_tab, gtdb_rep_tab=None, rank=None, representatives_source=None):
    """ write accessions (+ metadata TSV) for a taxon, from a DataFrame """
    working_tab = gtdb_rep_tab if representatives_source else gtdb_tab

    if taxon == "all":
        if representatives_source:
            tab_out_filename = "gtdb-arc-and-bac-refseq-rep-metadata.tsv"
            acc_out_filename = "gtdb-arc-and-bac-refseq-rep-accessions.txt"
            working_tab.to_csv(tab_out_filename, sep="\t", index=False)
        else:
            acc_out_filename = "gtdb-arc-and-bac-accessions.txt"
            tab_out_filename = None
    else:
        ranks_found_in = find_ranks_for_taxon(taxon, working_tab)

        if len(ranks_found_in) == 0:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
            print("")
            sys.exit(0)

        if not rank:
            if len(ranks_found_in) > 1:
                wprint(color_text("Since '" + str(taxon) + "' occurs at more than 1 rank, we'll need to specify "
                       "which rank is wanted as well before we pull the accessions. This can be specified with the `-r` flag.", "yellow"))
                print("")
                sys.exit(0)
            else:
                rank = ranks_found_in[0]
        else:
            rank = rank.lower()

        working_tab = working_tab[working_tab[rank] == taxon]

        taxon_for_filename = taxon.replace(" ", "-").replace("/", "-").lower()

        if representatives_source:
            tab_out_filename = "gtdb-" + taxon_for_filename + "-" + rank + "-" + representatives_source + "-rep-metadata.tsv"
            acc_out_filename = "gtdb-" + taxon_for_filename + "-" + rank + "-" + representatives_source + "-rep-accs.txt"
        else:
            tab_out_filename = "gtdb-" + taxon_for_filename + "-" + rank + "-metadata.tsv"
            acc_out_filename = "gtdb-" + taxon_for_filename + "-" + rank + "-accs.txt"

        working_tab.to_csv(tab_out_filename, sep="\t", index=False)

    target_accs = working_tab["ncbi_genbank_assembly_accession"].tolist()

    with open(acc_out_filename, "w") as out:
        for acc in target_accs:
            out.write(acc + "\n")

    print("")
    wprint(f"Wrote {len(target_accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    if tab_out_filename:
        wprint("Associated taxonomy and metadata of these targets written to:")
        wprint("  " + color_text(tab_out_filename))
        print("")
    else:
        wprint("The GTDB table that already exists holds all of these: " + color_text("gtdb-arc-and-bac-metadata.tsv"))
        print("")


def _write_accessions(taxon, working_df, rank, representatives_source=None):
    """ thin orchestrator->worker adapter: working_df is already the taxon slice """
    if representatives_source:
        get_accessions(taxon, None, gtdb_rep_tab=working_df, rank=rank,
                       representatives_source=representatives_source)
    else:
        get_accessions(taxon, working_df, rank=rank)


def get_unique_taxon_counts(taxon, gtdb_tab, gtdb_rep_tab=None, representatives_source=None):
    """ counts of a specific taxon (or all) per rank, from a DataFrame """
    if taxon == "all":
        count = len(gtdb_tab.index)
        print("")
        wprint("  There are " + str(count) + " total genomes in the database.")
        print("")
        if representatives_source:
            count = len(gtdb_rep_tab.index)
            rep_source = "GTDB" if representatives_source == "gtdb" else "RefSeq"
            wprint(color_text("In considering only " + rep_source + " representative genomes:", "yellow"))
            print("")
            wprint("  There are " + str(count) + " total representative genomes in the database.")
            print("")
    else:
        ranks_found_in = find_ranks_for_taxon(taxon, gtdb_tab)

        for rank in ranks_found_in:
            count = len(gtdb_tab[gtdb_tab[rank] == taxon].index)
            wprint("  The rank '" + rank + "' has " + str(count) + " " + taxon + " entries.")

        if len(ranks_found_in) == 0:
            wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
            print("")
            sys.exit(0)

        print("")

        if representatives_source:
            ranks_found_in_rep = find_ranks_for_taxon(taxon, gtdb_rep_tab)
            rep_source = "GTDB" if representatives_source == "gtdb" else "RefSeq"
            wprint(color_text("In considering only " + rep_source + " representative genomes:", "yellow"))
            print("")
            for rank in ranks_found_in_rep:
                count = len(gtdb_rep_tab[gtdb_rep_tab[rank] == taxon].index)
                wprint("  The rank '" + rank + "' has " + str(count) + " " + taxon + " representative genome entries.")
                print("")
            if len(ranks_found_in_rep) == 0:
                wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank as a representative genome :(", "yellow"))
                print("")
                sys.exit(0)


def get_unique_taxa_counts_of_all_ranks(gtdb_tab, gtdb_rep_tab=None, representatives_source=None):
    """ counts of unique taxa at each rank, from a DataFrame """
    print("\n    {:<10} {:}".format("Rank", "Num. Unique Taxa"))
    for rank in RANKS:
        print("    {:<10} {:}".format(rank, str(gtdb_tab[rank].nunique())))
    print("")

    if representatives_source == "gtdb":
        wprint(color_text("(The `--gtdb-representatives-only` flag doesn't change these "
                          "counts: every GTDB taxon has a representative genome, so the "
                          "number of unique taxa per rank is the same with or without it.)",
                          "yellow"))
        print("")
    elif representatives_source == "refseq":
        wprint(color_text("In considering only RefSeq reference genomes:", "yellow"))
        print("")
        print("    {:<10} {:}".format("Rank", "Num. Unique Ref. Taxa"))
        for rank in RANKS:
            print("    {:<10} {:}".format(rank, str(gtdb_rep_tab[rank].nunique())))
        print("")
