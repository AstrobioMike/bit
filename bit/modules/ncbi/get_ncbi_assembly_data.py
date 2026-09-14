#!/usr/bin/env python

"""
Ensures the prepared NCBI assembly-info table is present, downloading bit's hosted
Parquet asset (ncbi-data.parquet) if it isn't. The asset is a combined, slimmed
GenBank + RefSeq assembly summary with taxonomy resolved.

The download/verify/cleanup machinery is shared with the GTDB table, see
bit/modules/hosted_parquet_asset.py. What lives here is the NCBI-specific
configuration and the names the rest of the codebase imports.
"""

import os

from bit.modules.general import report_message
from bit.modules.hosted_parquet_asset import (HostedParquetAsset,
                                              validate_date_stamp)
from bit.modules.ncbi.build_ncbi_data_parquet import PARQUET_FILENAME, DATE_FILENAME


_RELEASE_BASE = "https://github.com/AstrobioMike/bit/releases/download/ncbi-assembly-info-latest"

NCBI_ASSET = HostedParquetAsset(
    env_var="NCBI_assembly_data_dir",
    release_base=_RELEASE_BASE,
    parquet_filename=PARQUET_FILENAME,
    sidecar_filename=DATE_FILENAME,
    sidecar_validator=validate_date_stamp,
    display_name="NCBI assembly-info table",
    download_label="NCBI prepared data",
    sidecar_label="date stamp",
)

NCBI_DATA_URL = NCBI_ASSET.data_url
NCBI_DATE_URL = NCBI_ASSET.sidecar_url


def get_ncbi_assembly_data(force_update=False, quiet=False):
    """
    Ensure the NCBI Parquet table is present locally, and return its path.

    `quiet` silences the "already present" note only. A failed download always
    explains itself, since it's fatal and there's nothing else to go on.
    """
    ncbi_dir = check_ncbi_assembly_info_location_var_is_set()
    data_present = check_if_data_present(ncbi_dir)

    if data_present and not force_update:
        if not quiet:
            report_message("Assembly data already present at:")
            print(f"        {ncbi_dir}")
            report_message("Run `bit data get ncbi-assembly-data -f` if you want to re-download/update it.")
            print("")
    else:
        get_slim_ncbi_assembly_data(ncbi_dir)

    return ncbi_data_table_path(ncbi_dir)


def check_ncbi_assembly_info_location_var_is_set():
    return NCBI_ASSET.location()


def ncbi_data_table_path(location=None):
    """Path to the local NCBI Parquet asset (resolving the location if not given)."""
    return NCBI_ASSET.table_path(location)


def check_if_data_present(location):
    """
    True if both the Parquet table and date-retrieved.txt are present and non-empty.
    If either is missing/empty, any stray copy is cleaned up and we return False so a
    fresh copy is pulled.
    """
    return NCBI_ASSET.is_present(location)


def get_slim_ncbi_assembly_data(location):
    """
    Download the prepared NCBI Parquet asset and its date-retrieved file into
    `location`. The Parquet footer is verified before we trust the table, and the
    date file is written atomically. On any network/integrity failure the partial
    artifacts are cleaned up and we exit with a helpful message -- there is no
    NCBI-rebuild fallback, since the hosted asset is the prepared table.
    """
    NCBI_ASSET.download(location)


def read_date_retrieved(location):
    """
    Read date-retrieved.txt (a 'YYYY,MM,DD' stamp) from `location` and return it
    formatted like 'Jan 05, 2026'. Returns the raw string if it can't be parsed.
    """
    import datetime

    with open(os.path.join(location, DATE_FILENAME)) as fh:
        stamp = fh.readline().strip()
    try:
        y, m, d = (int(p) for p in stamp.split(","))
        return datetime.date(y, m, d).strftime("%b %d, %Y")
    except (ValueError, TypeError):
        return stamp
