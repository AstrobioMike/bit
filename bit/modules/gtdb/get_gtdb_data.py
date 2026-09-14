#!/usr/bin/env python

"""
Ensures the prepared GTDB metadata table is present, downloading bit's hosted
Parquet asset (gtdb-data.parquet) if it isn't. The Parquet asset is already slimmed
to the columns used here and already has the 7 taxonomic ranks split into their own
columns.

The download/verify/cleanup machinery is shared with the NCBI table, see
bit/modules/hosted_parquet_asset.py. What lives here is the GTDB-specific
configuration and the names the rest of the codebase imports.
"""

import os

from bit.modules.general import report_message
from bit.modules.hosted_parquet_asset import (HostedParquetAsset,
                                              validate_version_lines)
from bit.modules.gtdb.build_gtdb_data_parquet import PARQUET_FILENAME, VERSION_FILENAME


# the prepared GTDB assets live on a rolling "-latest" GitHub release. The Parquet
# file is the slimmed, rank-split metadata table; VERSION.txt carries the GTDB
# release + date lines.
_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest"

GTDB_ASSET = HostedParquetAsset(
    env_var="GTDB_DIR",
    release_base=_RELEASE_BASE,
    parquet_filename=PARQUET_FILENAME,
    sidecar_filename=VERSION_FILENAME,
    sidecar_validator=validate_version_lines,
    display_name="GTDB table",
    download_label="GTDB prepared data",
    sidecar_label="version info",
)

GTDB_DATA_URL = GTDB_ASSET.data_url
GTDB_VERSION_URL = GTDB_ASSET.sidecar_url


def get_gtdb_data(force_update=False, quiet=False):
    """
    Ensure the GTDB Parquet table is present locally, and return its path.

    `quiet` silences the "already present" note only. A failed download always
    explains itself, since it's fatal and there's nothing else to go on.
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
        get_slim_gtdb_tab(gtdb_dir)

    return gtdb_data_table_path(gtdb_dir)


def check_gtdb_location_var_is_set():
    return GTDB_ASSET.location()


def gtdb_data_table_path(location=None):
    return GTDB_ASSET.table_path(location)


def check_if_gtdb_data_present(location):
    """
    True if both the Parquet table and version-info file are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False so a
    fresh copy is pulled.
    """
    return GTDB_ASSET.is_present(location)


def get_slim_gtdb_tab(location):
    """
    Download the prepared GTDB Parquet asset and its version-info file into
    `location`. The Parquet footer is verified before we trust the table, and the
    version file is written atomically. On any network/integrity failure the partial
    artifacts are cleaned up and we exit with a helpful message.
    """
    GTDB_ASSET.download(location)


def report_gtdb_version_info(location):
    """
    Return (version, release_date) from the local VERSION.txt (first two lines)
    """
    version_info = []
    with open(os.path.join(location, VERSION_FILENAME)) as version_info_file:
        for line in version_info_file:
            line = line.strip()
            if line != "":
                version_info.append(line.replace("Released ", ""))
    return version_info[0], version_info[1]
