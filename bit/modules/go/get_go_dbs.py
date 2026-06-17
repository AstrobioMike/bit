import sys
import os
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)


def check_go_data_location_var_is_set():

    try:
        go_data_dir = os.environ['GO_DB_DIR']
    except:
        wprint(color_text("The environment variable 'GO_DB_DIR' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `bit data locations check`.")
        print("")
        sys.exit(1)

    return go_data_dir


def check_if_data_present(location):

    go_basic_path = os.path.join(str(location), "go-basic.obo")
    goslim_metagenomics_path = os.path.join(str(location), "goslim_metagenomics.obo")

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(go_basic_path) or not is_nonempty_file(goslim_metagenomics_path):

        for p in (go_basic_path, goslim_metagenomics_path):
            if os.path.exists(p):
                if os.path.isfile(p):
                    os.remove(p)
        return False
    return True


def download_go_data(location):

    # setting links
    go_basic_link = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    goslim_metagenomics_link = "http://current.geneontology.org/ontology/subsets/goslim_metagenomics.obo"

    go_basic_path = os.path.join(str(location), "go-basic.obo")
    goslim_metagenomics_path = os.path.join(str(location), "goslim_metagenomics.obo")

    print(color_text("\n    Downloading GO databases...\n", "yellow"))

    try:
        download_with_tqdm(go_basic_link, "        GO basic", go_basic_path)
        download_with_tqdm(goslim_metagenomics_link, "        GO slim metagenomics", goslim_metagenomics_path)
        print("")
    except Exception as e:
        report_message(f"Downloading the GO databases failed with the following error:\n{e}", "red")
        notify_premature_exit()


def get_go_data(force_update=False, quiet=False):

    go_db_dir = check_go_data_location_var_is_set()
    data_present = check_if_data_present(go_db_dir)

    if data_present and not force_update:
        if not quiet:
            report_message(f"GO data already present at:")
            print(f"        {go_db_dir}")
            report_message(f"Run `bit data get go-dbs -f` if you want to re-download/update it.")
            print()
            return
        return
    else:
        download_go_data(go_db_dir)
