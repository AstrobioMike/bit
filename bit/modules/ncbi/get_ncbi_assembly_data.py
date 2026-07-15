import os
import sys
import socket
import urllib
import urllib.error
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)



_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/ncbi-assembly-info-latest"

TABLE_FILENAME = "ncbi-data.parquet"
DATE_FILENAME = "date-retrieved.txt"

NCBI_DATA_URL = f"{_RELEASE_BASE}/{TABLE_FILENAME}"
NCBI_DATE_URL = f"{_RELEASE_BASE}/{DATE_FILENAME}"


def check_ncbi_assembly_info_location_var_is_set():

    # making sure there is a NCBI_assembly_data_dir env variable
    try:
        ncbi_assembly_data_dir = os.environ['NCBI_assembly_data_dir']
    except KeyError:
        wprint(color_text("The environment variable 'NCBI_assembly_data_dir' does not seem to be set :(", "yellow"))
        wprint("This shouldn't happen, check on things with `bit data locations check`.")
        print("")
        sys.exit(0)

    return ncbi_assembly_data_dir


def ncbi_data_table_path(location=None):
    if location is None:
        location = check_ncbi_assembly_info_location_var_is_set()
    return os.path.join(str(location), TABLE_FILENAME)


def check_if_data_present(location):

    table_path = os.path.join(str(location), TABLE_FILENAME)
    date_retrieved_path = os.path.join(str(location), DATE_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(table_path) or not is_nonempty_file(date_retrieved_path):
        for p in (table_path, date_retrieved_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False
    return True


def _report_unavailable(err):
    print("")
    wprint(color_text("Couldn't download the prepared NCBI assembly-info table :(", "yellow"))
    report_message(f"Underlying issue: {err}", initial_indent="    ",
                   subsequent_indent="    ")
    print("")
    wprint("This is usually a transient network problem, and trying again in a few minutes "
           "often works. If it persists, the table can be fetched manually from:")
    print(f"        {color_text(NCBI_DATA_URL)}")
    print(f"        {color_text(NCBI_DATE_URL)}")
    wprint(f"and placed (as '{TABLE_FILENAME}' and '{DATE_FILENAME}') in the directory "
           "shown by `bit data locations check`.")
    print("")


def get_slim_ncbi_assembly_data(location, quiet=False):
    table_path = os.path.join(location, TABLE_FILENAME)
    date_path = os.path.join(location, DATE_FILENAME)

    print(color_text("\n    Downloading the prepared NCBI assembly-info table (only needs to be done once)...\n", "yellow"))

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        download_with_tqdm(NCBI_DATA_URL, "        NCBI prepared data", table_path,
                           speed_gate=True)

        # confirm the file is a readable Parquet before we trust it
        _verify_parquet(table_path)

        _download_date_file(date_path)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            ValueError, OSError) as err:
        for p in (table_path, date_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        if not quiet:
            _report_unavailable(err)
        notify_premature_exit()
        return
    finally:
        socket.setdefaulttimeout(default_timeout)

    print("")


def _verify_parquet(path):
    """
    Cheap integrity check: open the Parquet metadata (footer) and confirm the file
    has a schema and at least one row group. This reads only the footer, not the
    whole table, and turns a truncated/corrupt download into a clean failure.
    """
    import pyarrow.parquet as pq # type: ignore
    md = pq.ParquetFile(path).metadata
    if md.num_columns == 0 or md.num_row_groups == 0:
        raise ValueError("downloaded NCBI table has no data (truncated download?)")


def _download_date_file(date_path):
    tmp = date_path + ".part"
    try:
        download_with_tqdm(NCBI_DATE_URL, "        date stamp", tmp, leave=False)
        _validate_date_file(tmp)
        os.replace(tmp, date_path)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def _validate_date_file(path):
    with open(path) as fh:
        first = fh.readline().strip()
    parts = first.split(",")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise ValueError(f"date-retrieved.txt is not a 'YYYY,MM,DD' stamp: {first!r}")


def get_ncbi_assembly_data(force_update=False, quiet=False):
    """
    Ensure the NCBI Parquet table is present locally, and return its path
    """
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        if not quiet:
            report_message("Assembly data already present at:")
            print(f"        {ncbi_dir}")
            report_message("Run `bit data get ncbi-assembly-data -f` if you want to re-download/update it.")
            print()
    else:
        get_slim_ncbi_assembly_data(ncbi_dir, quiet=quiet)

    return ncbi_data_table_path(ncbi_dir)
