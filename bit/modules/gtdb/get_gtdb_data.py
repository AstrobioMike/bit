import sys
import os
import socket
import pandas as pd
import urllib
import urllib.error
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)


GTDB_BASE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"


def get_gtdb_data(force_update=False, quiet=False):
    GTDB_dir = check_gtdb_location_var_is_set()
    data_present = check_if_gtdb_data_present(GTDB_dir)

    if data_present and not force_update:
        if not quiet:
            report_message("GTDB data already present at:")
            print(f"        {GTDB_dir}")
            report_message("Run `bit data get gtdb-data -f` if you want to re-download/update it.")
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
        wprint("This shouldn't happen, check on things with `bit data-locations check`.")
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


def report_gtdb_unreachable(err):

    print("")
    wprint(color_text("  The GTDB data server could not be reached :(", "red"))
    report_message(f"While trying to download from: {GTDB_BASE_URL}", color = "none",
                   initial_indent = "    ", subsequent_indent = "      ")
    report_message("This is typically a network/connectivity problem on the system running "
           "bit (e.g., no internet access, a firewall or proxy blocking the "
           "connection, or the GTDB server being temporarily unavailable), rather "
           "than a problem with bit itself.", color = "none", initial_indent = "    ",
           subsequent_indent = "    ")
    report_message(("Things to try:"), color = "none", initial_indent = "    ", subsequent_indent = "    ")
    print("          - confirm the system has internet access")
    print(f"          - confirm {GTDB_BASE_URL} is reachable (e.g., in a browser or with curl)")
    print("          - if behind a proxy/firewall, check that outbound HTTPS is allowed")
    print("          - wait and try again later in case the server is temporarily down")
    report_message(f"Underlying error: {err}", initial_indent = "    ", subsequent_indent = "    ")
    print("")
    notify_premature_exit()
    sys.exit(1)


def gen_gtdb_tab(location):
    """ downloads and parses the GTDB info tables """

    print(color_text("\n    Downloading and parsing metadata tables from GTDB (only needs to be done once)...\n", "yellow"))

    arc_path = os.path.join(location, "ar53_metadata.tsv.gz")
    bac_path = os.path.join(location, "bac120_metadata.tsv.gz")
    metadata_path = os.path.join(location, "GTDB-arc-and-bac-metadata.tsv")
    version_info_path = os.path.join(location, "GTDB-version-info.txt")

    # fail fast (instead of hanging) if the server can't be reached
    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)

    try:
        # getting archaea
        arc_link = f"{GTDB_BASE_URL}/ar53_metadata.tsv.gz"
        arc_tsv_gz = download_with_tqdm(arc_link, "        GTDB archaeal data", arc_path)
        arc_tab = pd.read_csv(arc_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
        arc_tab.rename(columns={arc_tab.columns[0]:"accession"}, inplace=True)
        arc_tab.dropna(inplace=True, how="all")

        # getting bacteria
        bac_link = f"{GTDB_BASE_URL}/bac120_metadata.tsv.gz"
        bac_tsv_gz = download_with_tqdm(bac_link, "        GTDB bacterial data", bac_path)
        bac_tab = pd.read_csv(bac_tsv_gz, sep="\t", compression="gzip", on_bad_lines = 'skip', header=0, low_memory=False)
        bac_tab.rename(columns={bac_tab.columns[0]:"accession"}, inplace=True)
        bac_tab.dropna(inplace=True, how="all")
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as err:
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)

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

    try:
        socket.setdefaulttimeout(30)
        urllib.request.urlretrieve(f"{GTDB_BASE_URL}/VERSION.txt", version_info_path)
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as err:
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)
