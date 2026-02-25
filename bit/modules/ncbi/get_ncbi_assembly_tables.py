#!/usr/bin/env python

"""
This is a helper module of bit taken from my GToTree package (https://github.com/AstrobioMike/GToTree/wiki)
to download the NCBI assembly summary tables if they are not present, or to update them.
"""

import sys
import os
import argparse
import shutil
from datetime import date
from bit.cli.common import CustomRichHelpFormatter
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)



################################################################################

def main():

    parser = argparse.ArgumentParser(description="This is a helper program to download and setup the NCBI assembly summary tables if they are \
                                              not present, or to update them.", \
                                 epilog="Ex. usage: `get-ncbi-assembly-tables`",
                                 formatter_class=CustomRichHelpFormatter)

    parser.add_argument("-f", "--force-update", help='Update the stored NCBI assembly tables', action = "store_true")

    args = parser.parse_args()

    get_ncbi_assembly_data(force_update = args.force_update)

################################################################################


def check_ncbi_assembly_info_location_var_is_set():

    # making sure there is a NCBI_assembly_data_dir env variable
    try:
        NCBI_data_dir = os.environ['NCBI_assembly_data_dir']
    except:
        wprint(color_text("The environment variable 'NCBI_assembly_data_dir' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `bit-data-locations check`.")
        print("")
        sys.exit(0)

    return(NCBI_data_dir)


def check_if_data_present(location):

    # seeing if present already
    table_path = os.path.join(str(location), "ncbi-assembly-info.tsv")
    date_retrieved_path = os.path.join(str(location), "date-retrieved.txt")

    # if either file is missing, we are going to download, we also package the date-retrieved file empty with conda to retain directory, so checking it's not empty as well
    if not os.path.isfile(table_path) or not os.path.isfile(date_retrieved_path) or not os.path.getsize(date_retrieved_path) > 0:

        if os.path.exists(table_path):
            os.remove(table_path)
        if os.path.isdir(date_retrieved_path):
            shutil.rmtree(date_retrieved_path)

        return(False)

    return(True)


def get_NCBI_assembly_summary_data(location):

    """ downloads the needed ncbi assembly summary tables and combines them """

    # setting links
    genbank_link = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
    refseq_link = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

    table_path = os.path.join(str(location), "ncbi-assembly-info.tsv")
    refseq_temp_path = os.path.join(str(location), "refseq-assembly-info.tmp")

    print(color_text("\n    Downloading NCBI assembly summaries (only done once)...\n", "yellow"))

    try:
        download_with_tqdm(genbank_link, "        Genbank assemblies summary", table_path)
        download_with_tqdm(refseq_link, "        RefSeq assemblies summary", refseq_temp_path)
        print("")
    except Exception as e:
        report_message(f"Downloading the NCBI assembly summary tables failed with the following error:\n{e}", "red")
        notify_premature_exit()

    # combining
    with open (table_path, "a") as final_table:
        with open(refseq_temp_path, "r") as refseq:
            final_table.write(refseq.read())

    # removing temp
    if os.path.exists(refseq_temp_path):
        os.remove(refseq_temp_path)

    # storing date retrieved
    date_retrieved = str(date.today()).replace("-", ",")
    date_retrieved.replace("-", ",")

    date_retrieved_path = os.path.join(str(location), "date-retrieved.txt")

    with open(date_retrieved_path, "w") as outfile:
        outfile.write(date_retrieved + "\n")

def get_ncbi_assembly_data(force_update=False):
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        return
    else:
        get_NCBI_assembly_summary_data(ncbi_dir)


################################################################################

if __name__ == "__main__":
    main()