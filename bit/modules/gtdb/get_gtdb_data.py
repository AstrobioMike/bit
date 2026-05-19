#!/usr/bin/env python

import sys
import os
import pandas as pd
import urllib
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)


################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="This is a bit helper program to download or update the stored GTDB metadata.",
        epilog="Ex. usage: `get-gtdb-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    parser.add_argument("-q", "--quiet", help="Exit silently if GTDB data already present", action = "store_true")
    parser.add_argument("-f", "--force-update", help='Update the stored GTDB metadata', action = "store_true")

    add_help(parser)

    args = parser.parse_args()

    get_gtdb_data(force_update=args.force_update, quiet=args.quiet)

################################################################################

def get_gtdb_data(force_update=False, quiet=False):
    GTDB_dir = check_gtdb_location_var_is_set()
    data_present = check_if_gtdb_data_present(GTDB_dir)

    if data_present and not force_update:
        if not quiet:
            report_message("GTDB data already present at:")
            print(f"        {GTDB_dir}")
            report_message("Run `get-gtdb-data -f` if you want ot re-download/update the GTDB data.")
            print("")
            return GTDB_dir
        return GTDB_dir
    else:
        gen_gtdb_tab(GTDB_dir)
        return GTDB_dir


def check_gtdb_location_var_is_set():
    try:
        gtdb_data_dir = os.environ['GTDB_DIR']
    except KeyError:
        wprint(color_text("The environment variable 'GTDB_DIR' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `bit-data-locations check`.")
        sys.exit(1)
    return gtdb_data_dir


def check_if_gtdb_data_present(location):

    metadata_path = os.path.join(str(location), "GTDB-arc-and-bac-metadata.tsv")
    version_info_path = os.path.join(str(location), "GTDB-version-info.txt")

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(metadata_path) or not is_nonempty_file(version_info_path):

        for p in (metadata_path, version_info_path):
            if os.path.exists(p):
                if os.path.isfile(p):
                    os.remove(p)
        return False
    return True


def gen_gtdb_tab(location):
    """ downloads and parses the GTDB info tables """

    print(color_text("\n    Downloading and parsing metadata tables from GTDB (only needs to be done once)...\n", "yellow"))

    arc_path = os.path.join(location, "ar53_metadata.tsv.gz")
    bac_path = os.path.join(location, "bac120_metadata.tsv.gz")
    metadata_path = os.path.join(location, "GTDB-arc-and-bac-metadata.tsv")
    version_info_path = os.path.join(location, "GTDB-version-info.txt")

    # getting archaea
    arc_link = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz"
    arc_tsv_gz = download_with_tqdm(arc_link, "        GTDB archaeal data", arc_path)
    arc_tab = pd.read_csv(arc_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    arc_tab.rename(columns={arc_tab.columns[0]:"accession"}, inplace=True)
    arc_tab.dropna(inplace=True, how="all")

    # getting bacteria
    bac_link = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz"
    bac_tsv_gz = download_with_tqdm(bac_link, "        GTDB bacterial data", bac_path)
    bac_tab = pd.read_csv(bac_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    bac_tab.rename(columns={bac_tab.columns[0]:"accession"}, inplace=True)
    bac_tab.dropna(inplace=True, how="all")

    # combining
    gtdb_tab = pd.concat([arc_tab, bac_tab])

    print("")

    # splitting gtdb taxonomy column into 7 and dropping the single column
    domain, phylum, rclass, order, family, genus, species = [], [], [], [], [], [], []

    for index, row in gtdb_tab.iterrows():
        curr_acc = row["accession"]
        tax_list = row["gtdb_taxonomy"].split(";")

        if len(tax_list) != 7:
            wprint(color_text("GTDB entry " + curr_acc + " doesn't seem to have 7-column lineage info. Something is likely wrong :(", "yellow"))
            print("")
            wprint("If this continues to happen, please file an issue at github.com/AstrobioMike/GToTree/issues")
            print("")
            wprint("Aborting for now.")
            print("")
            sys.exit(0)

        else:
            domain.append(tax_list[0][3:])
            phylum.append(tax_list[1][3:])
            rclass.append(tax_list[2][3:])
            order.append(tax_list[3][3:])
            family.append(tax_list[4][3:])
            genus.append(tax_list[5][3:])
            species.append(tax_list[6][3:])

    gtdb_tab.insert(1, "species", species)
    gtdb_tab.insert(1, "genus", genus)
    gtdb_tab.insert(1, "family", family)
    gtdb_tab.insert(1, "order", order)
    gtdb_tab.insert(1, "class", rclass)
    gtdb_tab.insert(1, "phylum", phylum)
    gtdb_tab.insert(1, "domain", domain)

    # writing out
    gtdb_tab.to_csv(metadata_path, index=False, sep="\t")

    urllib.request.urlretrieve("https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION.txt", version_info_path)


if __name__ == "__main__":
    main()
