import sys
import os
import shutil
import pandas as pd # type: ignore
from bit.modules.general import color_text, wprint
from bit.modules.gtdb.get_gtdb_data import get_gtdb_data


def get_accessions_from_gtdb(args):

    GTDB_dir = get_gtdb_data(quiet=True)

    if args.get_table:
        copy_gtdb_table(GTDB_dir)
        sys.exit(0)

    # some checks to prevent things that either should be provided together, or should be mutually exclusive
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

    # reading in the GTDB table
    gtdb_tab = read_gtdb_tab(GTDB_dir)

    gtdb_rep_tab = None
    representatives_source = None

    if args.gtdb_representatives_only:
        gtdb_rep_tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]
        representatives_source = "GTDB"

    if args.refseq_reference_genomes_only:
        gtdb_rep_tab = gtdb_tab[gtdb_tab["ncbi_refseq_category"] == "reference genome"]
        representatives_source = "RefSeq"

    # resolve taxon name case-insensitively if provided
    target_taxon = resolve_taxon_name(args.target_taxon, gtdb_tab) if args.target_taxon else None

    if args.get_rank_counts:
        get_unique_taxa_counts_of_all_ranks(gtdb_tab, gtdb_rep_tab, representatives_source=representatives_source)
        sys.exit(0)

    if args.get_taxon_counts:
        get_unique_taxon_counts(target_taxon, gtdb_tab, gtdb_rep_tab, representatives_source=representatives_source)
        sys.exit(0)

    if target_taxon:
        get_accessions(target_taxon, gtdb_tab, gtdb_rep_tab, rank=args.target_rank, representatives_source=representatives_source)
        sys.exit(0)


def report_gtdb_version_info(location):
    """ reporting GTDB version info """

    version_info = []

    with open(os.path.join(location, "GTDB-version-info.txt")) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line)

    gtdb_version = version_info[0]
    gtdb_release_date = version_info[1]

    print("    Using GTDB " + gtdb_version + ": " + gtdb_release_date + "\n")


def copy_gtdb_table(GTDB_dir):
    """ copies the GTDB metadata table to the current directory """

    report_gtdb_version_info(GTDB_dir)

    shutil.copy(os.path.join(GTDB_dir, "GTDB-arc-and-bac-metadata.tsv"), "GTDB-arc-and-bac-metadata.tsv")
    print("")
    report_gtdb_version_info(GTDB_dir)
    wprint("GTDB table copied to:")
    print(color_text("    GTDB-arc-and-bac-metadata.tsv\n"))


def read_gtdb_tab(GTDB_dir):
    """ reads in file with gtdb info """

    print("")
    wprint(color_text("Reading in the GTDB info table...", "yellow"))

    report_gtdb_version_info(GTDB_dir)

    gtdb_tab = pd.read_csv(os.path.join(GTDB_dir, "GTDB-arc-and-bac-metadata.tsv"), sep="\t", low_memory=False)

    return gtdb_tab


def resolve_taxon_name(taxon, gtdb_tab):
    """ resolves user-provided taxon to canonical GTDB name via case-insensitive matching """

    if taxon.lower() == "all":
        return "all"

    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # build a lookup of lowercase -> canonical name across all rank columns
    canonical_names = {}
    for rank in ranks:
        for name in gtdb_tab[rank].unique():
            canonical_names[name.lower()] = name

    taxon_lower = taxon.lower()

    if taxon_lower in canonical_names:
        canonical = canonical_names[taxon_lower]
        if canonical != taxon:
            wprint(color_text("Matched input '" + taxon + "' to GTDB taxon '" + canonical + "'.", "yellow"))
            print("")
        return canonical

    wprint(color_text("Input taxon '" + taxon + "' doesn't seem to exist at any rank :(", "yellow"))
    print("")
    sys.exit(0)


def find_ranks_for_taxon(taxon, tab):
    """ returns list of ranks that contain the given taxon """

    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    return [rank for rank in ranks if taxon in tab[rank].unique()]


def get_accessions(taxon, gtdb_tab, gtdb_rep_tab=None, rank=None, representatives_source=None):
    """ get accessions based on specified taxon, if the provided taxon is in more than one rank, will require specified rank """

    # figure out which table to pull from
    working_tab = gtdb_rep_tab if representatives_source else gtdb_tab

    if taxon == "all":

        if representatives_source:
            tab_out_filename = "GTDB-arc-and-bac-refseq-rep-metadata.tsv"
            acc_out_filename = "GTDB-arc-and-bac-refseq-rep-accessions.txt"
            working_tab.to_csv(tab_out_filename, sep="\t", index=False)
        else:
            acc_out_filename = "GTDB-arc-and-bac-accessions.txt"
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

        # swapping spaces for dashes in case it's a species taxon
        taxon_for_filename = taxon.replace(" ", "-")

        if representatives_source:
            tab_out_filename = "GTDB-" + taxon_for_filename + "-" + rank + "-" + representatives_source + "-rep-metadata.tsv"
            acc_out_filename = "GTDB-" + taxon_for_filename + "-" + rank + "-" + representatives_source + "-rep-accs.txt"
        else:
            tab_out_filename = "GTDB-" + taxon_for_filename + "-" + rank + "-metadata.tsv"
            acc_out_filename = "GTDB-" + taxon_for_filename + "-" + rank + "-accs.txt"

        working_tab.to_csv(tab_out_filename, sep="\t", index=False)

    # write accessions
    target_accs = working_tab["ncbi_genbank_assembly_accession"].tolist()

    with open(acc_out_filename, "w") as out:
        for acc in target_accs:
            out.write(acc + "\n")

    # report results
    print("")

    wprint(f"Wrote {len(target_accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    if tab_out_filename:
        wprint("Associated taxonomy and metadata of these targets written to:")
        wprint("  " + color_text(tab_out_filename))
        print("")
    else:
        wprint("The GTDB table that already exists holds all of these: " + color_text("GTDB-arc-and-bac-metadata.tsv"))
        print("")


def get_unique_taxon_counts(taxon, gtdb_tab, gtdb_rep_tab=None, representatives_source=None):
    """ get counts of specific taxa """

    if taxon == "all":
        count = len(gtdb_tab.index)
        print("")
        wprint("  There are " + str(count) + " total genomes in the database.")
        print("")

        if representatives_source:
            count = len(gtdb_rep_tab.index)
            wprint(color_text("In considering only " + representatives_source + " representative genomes:", "yellow"))
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

            wprint(color_text("In considering only " + representatives_source + " representative genomes:", "yellow"))
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
    """ get counts of unique taxa at each rank """

    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    print("\n    {:<10} {:}".format("Rank", "Num. Unique Taxa"))
    for rank in ranks:
        print("    {:<10} {:}".format(rank, str(gtdb_tab[rank].nunique())))

    print("")

    # below only needed for RefSeq, because if it is GTDB, it is equivalent (i.e., they have 1 rep genome for each unique rank in their system)
    if representatives_source == "RefSeq":

        wprint(color_text("In considering only " + representatives_source + " representative genomes:", "yellow"))
        print("")

        print("    {:<10} {:}".format("Rank", "Num. Unique Rep. Taxa"))
        for rank in ranks:
            print("    {:<10} {:}".format(rank, str(gtdb_rep_tab[rank].nunique())))

        print("")
