#!/usr/bin/env python

"""
This is a helper module of bit taken from my GToTree package (https://github.com/AstrobioMike/GToTree/wiki)
to download NCBI tax data (largely for using TaxonKit (https://bioinf.shenwei.me/taxonkit/)).
"""

import sys
import os
import argparse
import tarfile
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import (wprint, color_text, report_message,
                                 notify_premature_exit, download_with_tqdm)

################################################################################

def main():

    parser = argparse.ArgumentParser(
        description="This is a bit helper program to download NCBI taxonomy data.",
        epilog="Example usage: `get-ncbi-tax-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    parser.add_argument("-q", "--quiet", help="Exit silently if tax data already present", action = "store_true")
    parser.add_argument("-f", "--force-update", help='Update the stored NCBI taxonomy data', action = "store_true")

    add_help(parser)

    args = parser.parse_args()

    get_ncbi_tax_data(force_update=args.force_update, quiet=args.quiet)

################################################################################

def check_tax_location_var_is_set():
    try:
        ncbi_tax_data_dir = os.environ['TAXONKIT_DB']
    except KeyError:
        wprint(color_text("The environment variable 'TAXONKIT_DB' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `bit-data-locations check`.")
        print("")
        sys.exit(1)

    return ncbi_tax_data_dir


def check_if_data_present(location):

    names_path = os.path.join(location, "names.dmp")
    nodes_path = os.path.join(location, "nodes.dmp")

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not (is_nonempty_file(names_path) and is_nonempty_file(nodes_path)):

        # clean up any partial or empty files so the helper can redownload cleanly
        for p in (names_path, nodes_path):
            if os.path.exists(p):
                if os.path.isfile(p):
                    os.remove(p)
        return False
    return True


def download_ncbi_tax_data(location):

    taxdump_path = os.path.join(location, "taxdump.tar.gz")
    taxdump_link = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    print(color_text("\n    Downloading required NCBI taxonomy data (only needs to be done once)...\n", "yellow"))

    try:
        download_with_tqdm(taxdump_link, "        NCBI Taxonomy data", taxdump_path)
        print("")
    except Exception as e:
        report_message(f"Downloading the NCBI taxonomy data failed with the following error:\n{e}", "red")
        notify_premature_exit()

    with tarfile.open(taxdump_path) as tarball:
        tarball.extractall(location)

    os.remove(taxdump_path)


def get_ncbi_tax_data(force_update=False, quiet = False):

    ncbi_data_dir = check_tax_location_var_is_set()
    data_present = check_if_data_present(ncbi_data_dir)

    if data_present and not force_update:
        if not quiet:
            report_message(f"Tax data already present at:")
            print(f"        {ncbi_data_dir}")
            report_message(f"Run `get-ncbi-tax-data -f` if you want to re-download/update it.")
            print()
            return
        return
    else:
        download_ncbi_tax_data(ncbi_data_dir)


################################################################################


if __name__ == "__main__":
    main()