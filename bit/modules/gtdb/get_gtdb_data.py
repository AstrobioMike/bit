import os
import sys
import socket
import urllib
import urllib.error
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)
from bit.modules.gtdb.build_gtdb_data_parquet import PARQUET_FILENAME, VERSION_FILENAME

_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest"

GTDB_DATA_URL = f"{_RELEASE_BASE}/{PARQUET_FILENAME}"
GTDB_VERSION_URL = f"{_RELEASE_BASE}/{VERSION_FILENAME}"


def check_gtdb_location_var_is_set():
    try:
        gtdb_data_dir = os.environ['GTDB_DIR']
    except KeyError:
        wprint(color_text("The environment variable 'GTDB_DIR' does not seem to be set :(", "red"))
        wprint("This shouldn't happen, check on things with `bit data-locations check`.")
        sys.exit(1)
    return gtdb_data_dir


def gtdb_data_table_path(location=None):
    if location is None:
        location = check_gtdb_location_var_is_set()
    return os.path.join(str(location), PARQUET_FILENAME)


def check_if_gtdb_data_present(location):

    table_path = os.path.join(str(location), PARQUET_FILENAME)
    version_info_path = os.path.join(str(location), VERSION_FILENAME)

    def is_nonempty_file(p):
        return os.path.isfile(p) and os.path.getsize(p) > 0

    if not is_nonempty_file(table_path) or not is_nonempty_file(version_info_path):
        for p in (table_path, version_info_path):
            if os.path.exists(p) and os.path.isfile(p):
                os.remove(p)
        return False
    return True


def _report_unavailable(err):
    print("")
    wprint(color_text("Couldn't download the prepared GTDB table :(", "yellow"))
    report_message(f"Underlying issue: {err}", initial_indent="    ",
                   subsequent_indent="    ")
    print("")
    wprint("This is usually a transient network problem, and trying again in a few minutes "
           "often works. If it persists, the table can be fetched manually from:")
    print(f"        {color_text(GTDB_DATA_URL)}")
    print(f"        {color_text(GTDB_VERSION_URL)}")
    wprint(f"and placed (as '{PARQUET_FILENAME}' and '{VERSION_FILENAME}') in the directory "
           "shown by `bit data-locations check`.")
    print("")


def get_slim_gtdb_tab(location, quiet=False):
    table_path = os.path.join(location, PARQUET_FILENAME)
    version_path = os.path.join(location, VERSION_FILENAME)

    print(color_text("\n    Downloading the prepared GTDB table (only needs to be done once)...\n", "yellow"))

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        download_with_tqdm(GTDB_DATA_URL, "        GTDB prepared data", table_path,
                           speed_gate=True)

        # confirm the file is a readable Parquet before we trust it
        _verify_parquet(table_path)

        _download_version_file(version_path)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            ValueError, OSError) as err:
        for p in (table_path, version_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        if not quiet:
            _report_unavailable(err)
        notify_premature_exit()
        sys.exit(1)
    finally:
        socket.setdefaulttimeout(default_timeout)

    print("")


def _verify_parquet(path):
    """
    Cheap integrity check: open the Parquet footer and confirm the file has a schema
    and at least one row group. Reads only the footer, not the whole table.
    """
    import pyarrow.parquet as pq # type: ignore
    md = pq.ParquetFile(path).metadata
    if md.num_columns == 0 or md.num_row_groups == 0:
        raise ValueError("downloaded GTDB table has no data (truncated download?)")


def _download_version_file(version_path):
    tmp = version_path + ".part"
    try:
        download_with_tqdm(GTDB_VERSION_URL, "        version info", tmp, leave=False)
        _validate_version_file(tmp)
        os.replace(tmp, version_path)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def _validate_version_file(path):
    lines = [ln.strip() for ln in open(path) if ln.strip()]
    if len(lines) < 2:
        raise ValueError("GTDB version file doesn't have the expected version + date lines")


def report_gtdb_version_info(location):
    version_info = []
    with open(os.path.join(location, VERSION_FILENAME)) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line)
    return version_info[0], version_info[1]


def get_gtdb_data(force_update=False, quiet=False):
    """
    Ensure the GTDB Parquet table is present locally, and return its path
    """
    gtdb_dir = check_gtdb_location_var_is_set()
    data_present = check_if_gtdb_data_present(gtdb_dir)

    if data_present and not force_update:
        if not quiet:
            report_message("GTDB data already present at:")
            print(f"        {gtdb_dir}")
            report_message("Run `bit data get gtdb-data -f` if you want to re-download/update it.")
            print("")
    else:
        get_slim_gtdb_tab(gtdb_dir, quiet=quiet)

    return gtdb_data_table_path(gtdb_dir)
