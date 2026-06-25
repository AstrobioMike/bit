import sys
import os
import socket
import shutil
import tarfile
import time
import pandas as pd # type: ignore
import urllib
import urllib.error
from bit.modules.general import (wprint, color_text,
                                 report_message, notify_premature_exit,
                                 download_with_tqdm)


GTDB_BASE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"

# pre-slimmed table + version-info, packaged as a single .tar.gz and hosted on a
# bit GitHub release for a fast default download. The archive contains exactly
# two files at its root:
#   GTDB-arc-and-bac-metadata.tsv   (already slimmed to GTDB_KEPT_COLUMNS)
#   GTDB-version-info.txt
# Used by default; `-f/--force-update` bypasses this and rebuilds from
# GTDB_BASE_URL instead (gen_gtdb_tab). If this is empty or the download fails,
# the default path falls back to that same upstream rebuild.
GTDB_SLIM_TARBALL_URL = "https://github.com/AstrobioMike/bit/releases/download/gtdb-metadata-latest/GTDB-arc-and-bac-metadata.tar.gz"

# the stored GTDB-arc-and-bac-metadata.tsv is slimmed to only the columns bit
# uses anywhere (see gen_metagenome.GTDB_USED_COLUMNS, which this must stay in
# sync with). The upstream GTDB metadata tables carry ~120 columns; keeping only
# these shrinks the stored table ~6x and speeds every downstream read. The order
# here is the on-disk column order. Any column absent from a given GTDB release
# is silently skipped at write time.
GTDB_KEPT_COLUMNS = [
    "accession", "ncbi_genbank_assembly_accession", "ncbi_taxid",
    "gtdb_representative", "ncbi_refseq_category",
    "domain", "phylum", "class", "order", "family", "genus", "species",
    "checkm2_completeness", "checkm2_contamination",
    "checkm_completeness", "checkm_contamination",
    "genome_size", "contig_count", "gc_count", "gc_percentage", "ambiguous_bases",
    "coding_bases", "coding_density",
]


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
        # both default and -f/--force-update fetch the prepared slim table from
        # bit's host; -f just forces a re-download even when local data already
        # exists (to pick up a newer GTDB release). If the asset download fails or
        # no tarball URL is configured, get_slim_gtdb_tab falls back to rebuilding
        # from the upstream GTDB release.
        get_slim_gtdb_tab(GTDB_dir, quiet=quiet)
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


def _download_with_retries(url, label, dest, quiet=False,
                           attempts=4, retry_wait=3):
    """
    download `url` to `dest`, retrying up to `attempts` times on transient
    failures (timeouts, connection resets, SSL handshake stalls, transient 5xx),
    with a short wait between tries. A 404 is not retried -- it's raised
    immediately so the caller can treat a missing asset distinctly. Raises the
    last error if all attempts fail.
    """
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            download_with_tqdm(url, label, dest)
            return
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise          # missing asset -> not transient, don't retry
            last_err = err     # transient 5xx etc. -> retry
        except (urllib.error.URLError, socket.timeout, TimeoutError,
                ConnectionError, OSError) as err:
            last_err = err     # timeout / reset / SSL stall -> retry

        if attempt < attempts:
            if not quiet:
                wprint(color_text(
                    f"    download failed (attempt {attempt}/{attempts}); retrying...",
                    "yellow"))
            time.sleep(retry_wait)

    raise last_err


def get_slim_gtdb_tab(location, quiet=False):
    """
    fast default path: download bit's pre-slimmed GTDB tarball and extract its
    two files (the slimmed metadata table and the version-info file) into
    `location`. Falls back to rebuilding from the upstream GTDB release
    (gen_gtdb_tab) if no tarball URL is configured or the download/extract fails,
    so the user always ends up with a usable table.
    """
    metadata_path = os.path.join(location, "GTDB-arc-and-bac-metadata.tsv")
    version_info_path = os.path.join(location, "GTDB-version-info.txt")

    if not GTDB_SLIM_TARBALL_URL:
        # no host configured -> rebuild from upstream
        gen_gtdb_tab(location)
        return

    if not quiet:
        print(color_text("\n    Downloading the prepared GTDB table (only needs to be done once)...\n", "yellow"))

    tarball_path = os.path.join(location, "GTDB-slim.tar.gz")
    expected = {"GTDB-arc-and-bac-metadata.tsv", "GTDB-version-info.txt"}

    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(30)
    try:
        _download_with_retries(
            GTDB_SLIM_TARBALL_URL, "        GTDB prepared data", tarball_path,
            quiet=quiet)

        # a truncated/corrupt download fails here (full-stream read via
        # getmembers), which trips the fallback below rather than writing a
        # corrupt table.
        with tarfile.open(tarball_path, "r:gz") as tar:
            members = tar.getmembers()
            # take only the two expected files, matched by basename, and guard
            # against path traversal / nested dirs by extracting to flat names.
            wanted = {}
            for m in members:
                base = os.path.basename(m.name)
                if base in expected and m.isfile():
                    wanted[base] = m
            missing = expected - set(wanted)
            if missing:
                raise ValueError(
                    f"prepared GTDB archive is missing expected file(s): {', '.join(sorted(missing))}")

            for base, m in wanted.items():
                src = tar.extractfile(m)
                if src is None:
                    raise ValueError(f"could not read '{base}' from prepared GTDB archive")
                dest = os.path.join(location, base)
                with src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)

    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
            tarfile.TarError, ValueError, OSError) as err:
        # clean up partial artifacts, then fall back to the upstream rebuild
        for p in (tarball_path, metadata_path, version_info_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        if not quiet:
            print("")
            wprint(color_text("  Couldn't get the prepared GTDB table; "))
            wprint(color_text("  rebuilding from the upstream GTDB release instead.", "yellow"))
            report_message(f"Underlying issue: {err}", initial_indent="    ",
                           subsequent_indent="    ")
            print("")
        gen_gtdb_tab(location)
        return
    finally:
        socket.setdefaulttimeout(default_timeout)
        if os.path.exists(tarball_path):
            try:
                os.remove(tarball_path)
            except OSError:
                pass
    print("")


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

    # writing out, slimmed to only the columns bit uses (intersect with present
    # columns so a release missing one of them doesn't error)
    keep = [c for c in GTDB_KEPT_COLUMNS if c in gtdb_tab.columns]
    gtdb_tab = gtdb_tab[keep]
    gtdb_tab.to_csv(metadata_path, index=False, sep="\t")

    try:
        socket.setdefaulttimeout(30)
        urllib.request.urlretrieve(f"{GTDB_BASE_URL}/VERSION.txt", version_info_path)
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as err:
        report_gtdb_unreachable(err)
    finally:
        socket.setdefaulttimeout(default_timeout)


# names of the two files packaged inside the distributed tarball
GTDB_TABLE_FILENAME = "GTDB-arc-and-bac-metadata.tsv"
GTDB_VERSION_FILENAME = "GTDB-version-info.txt"


def fetch_upstream_version(timeout=30):
    """
    return the current upstream GTDB VERSION.txt content as a string. Used by the
    refresh Action to decide whether a new release is available (by comparing
    against the version currently published alongside the hosted asset). Raises on
    network error so the caller can decide how to handle it.
    """
    import urllib.request
    req_url = f"{GTDB_BASE_URL}/VERSION.txt"
    default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(timeout)
    try:
        with urllib.request.urlopen(req_url) as resp:
            return resp.read().decode("utf-8", errors="replace")
    finally:
        socket.setdefaulttimeout(default_timeout)


def main(argv=None):
    """
    build the slim GTDB bundle for the refresh GitHub Action: download + slim the
    metadata via gen_gtdb_tab, then package the table + version file into a single
    .tar.gz, and (optionally) leave the version file alongside it so the Action can
    upload it as a separate small asset for cheap weekly version checks.

    Output in <out-dir>:
      <archive-name>            (default GTDB-arc-and-bac-metadata.tar.gz)
      GTDB-version-info.txt     (left in place for separate upload and rapid version check)
    """
    import argparse
    import tarfile

    ap = argparse.ArgumentParser(description="Build bit's slim GTDB bundle.")
    ap.add_argument("-o", "--out-dir", default=".", help="output directory")
    ap.add_argument("--archive-name", default="GTDB-arc-and-bac-metadata.tar.gz",
                    help="name of the produced .tar.gz")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    table_path = os.path.join(args.out_dir, GTDB_TABLE_FILENAME)
    version_path = os.path.join(args.out_dir, GTDB_VERSION_FILENAME)
    archive_path = os.path.join(args.out_dir, args.archive_name)

    # gen_gtdb_tab writes both GTDB-arc-and-bac-metadata.tsv (slimmed) and
    # GTDB-version-info.txt into out-dir.
    gen_gtdb_tab(args.out_dir)

    print(f"Packaging -> {archive_path}", flush=True)
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(table_path, arcname=GTDB_TABLE_FILENAME)
        tar.add(version_path, arcname=GTDB_VERSION_FILENAME)

    # remove the large uncompressed table; keep the version file loose so the
    # workflow can upload it as a separate small asset.
    if os.path.exists(table_path):
        os.remove(table_path)

    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
