#!/usr/bin/env python

"""
This is a helper program of bit, that I initially wrote for GToTree (https://github.com/AstrobioMike/GToTree/wiki).
It is for setting up reference files for the glorious Genome Taxonomy Database (gtdb.ecogenomic.org/).

For examples, please visit the GToTree wiki here: https://github.com/AstrobioMike/GToTree/wiki/example-usage
"""

import sys
import os
import urllib.request
import pandas as pd
import textwrap
import argparse

parser = argparse.ArgumentParser(description = "This is a helper program to facilitate setting up the reference files for the \
                                              glorious Genome Taxonomy Database (gtdb.ecogenomic.org). It's really meant for internal \
                                              use only by other bit programs.")

args = parser.parse_args()

################################################################################

def main():

    ## checking env variable is set and writable
    check_location_var_is_set_and_writable("GTDB_DIR")

    ## setting up ref GTDB files if needed
    check_and_or_get_gtdb_files(os.environ["GTDB_DIR"])

################################################################################


# setting some colors
tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


### functions ###
def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    print(textwrap.fill(text, width=80, initial_indent="  ",
          subsequent_indent="  ", break_on_hyphens=False))


def check_location_var_is_set_and_writable(variable):

    # making sure there is an env variable
    try:
        path = os.environ[variable]

        if path == "":
            raise

    except:
        print()
        wprint(color_text("The environment variable '" + str(variable) + "' does not seem to be set :(", "red"))
        print()
        wprint("Try to set it with `bit-data-locations set`, then try again.")
        print("\nExiting for now.\n")
        sys.exit(1)

    # making sure path is writable for the user
    path_writable = os.access(path, os.W_OK)

    if not path_writable:
        print()
        wprint(color_text("The environment variable '" + str(variable) + "' does not seem to be writable :(", "red"))
        print()
        wprint("Try to set it somewhere else with `bit-data-locations set`, then try again.")
        print("\nExiting for now.\n")
        sys.exit(1)

    return()


def gen_gtdb_tab(location):
    """ downloads and parses the GTDB info tables """

    # getting archaea
    arc_tar_gz = urllib.request.urlopen("https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz")
    arc_tab = pd.read_csv(arc_tar_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    arc_tab.rename(columns={arc_tab.columns[0]:"accession"}, inplace=True)
    arc_tab.dropna(inplace=True, how="all")

    # getting bacteria
    bac_tar_gz = urllib.request.urlopen("https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz")
    bac_tab = pd.read_csv(bac_tar_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
    bac_tab.rename(columns={bac_tab.columns[0]:"accession"}, inplace=True)
    bac_tab.dropna(inplace=True, how="all")

    # combining
    gtdb_tab = pd.concat([arc_tab, bac_tab])

    # splitting gtdb taxonomy column into 7 and dropping the single column
    domain, phylum, rclass, order, family, genus, species = [], [], [], [], [], [], []

    for index, row in gtdb_tab.iterrows():
        curr_acc = row["accession"]
        tax_list = row["gtdb_taxonomy"].split(";")

        if len(tax_list) != 7:
            wprint(color_text("GTDB entry " + curr_acc + " doesn't seem to have 7-column lineage info. Something is likely wrong :(", "yellow"))
            print("")
            wprint("If this continues to happen, please file an issue at github.com/AstrobioMike/bit/issues")
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
    gtdb_tab.to_csv(location + "GTDB-arc-and-bac-metadata.tsv", index=False, sep="\t")

    gtdb_version_info = urllib.request.urlretrieve("https://data.gtdb.ecogenomic.org/releases/latest/VERSION", location + "GTDB-version-info.txt")


def check_and_or_get_gtdb_files(GTDB_DIR):
    """ checks for and sets up ref GTDB files if needed """

    if os.path.exists(GTDB_DIR + "GTDB-arc-and-bac-metadata.tsv") and os.path.exists(GTDB_DIR + "GTDB-version-info.txt"):

        sys.exit(0)

    # generating when table doesn't exist yet
    else:
        print("")
        wprint(color_text("Downloading and parsing archaeal and bacterial metadata tables from GTDB (only needs to be done once, or when a new version is available)...", "yellow"))
        print("")

        gen_gtdb_tab(GTDB_DIR)


if __name__ == "__main__":
    main()
