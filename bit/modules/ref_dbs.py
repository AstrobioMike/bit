"""
The large reference databases the workflows use: where they are kept, whether they are set up, and setting them up.

They all live under one directory, each in its own subdirectory: the one given with -d/--ref-dbs-dir where a command
takes it, otherwise the one in the BIT_REF_DBS_DIR environment variable (set with `bit data locations set`). Each is
set up with `bit data get <name>`, or by the workflows before they start when one they need is missing.
"""

import json
import os
import platform
import shutil
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path

# kept light (only the standard library at import), as the command-line parsers read REF_DBS; the downloading
# parts are imported where they're used

REF_DBS_VAR = "BIT_REF_DBS_DIR"

# GTDB-Tk can't be installed on macOS (its pplacer dependency is only built for Linux); there the genomes and MAGs get
# their taxonomy from bit-tax-db instead, and this is said when a run starts
GTDBTK_RUNS_HERE = platform.system() != "Darwin"
GTDBTK_NOT_HERE_NOTE = ("GTDB-Tk can't be run on macOS (one of the programs it depends on, pplacer, isn't available for it), so bit is "
               "assigning the genomes' taxonomy itself, labelling each genome as a whole with bit-tax-db (as `bit classify-seqs` "
               "does). GTDB-Tk is still the best way to assign taxonomy to MAGs and genomes, so it's worth running this on "
               "Linux when you can.")

TAX_DB_VERSION = "2026-11"
TAX_DB_NAME = "bit-tax-db-" + TAX_DB_VERSION
TAX_DB_URL = f"https://github.com/AstrobioMike/bit-tax-db/releases/download/{TAX_DB_VERSION}/"

CHECKM2_DL_LINK = "https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content"
CHECKM2_RECORD = "https://zenodo.org/api/records/14897628"   # the Zenodo record, which gives the file's md5
CHECKM2_DB_FILENAME = "uniref100.KO.1.dmnd"

GTDB_DL_LINK = ("https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/auxillary_files/gtdbtk_package/"
                "full_package/gtdbtk_r232_data.tar.gz")
GTDB_MD5_LIST = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/MD5SUM.txt"
GTDB_MD5_ENTRY = "./auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"

# KEGG's KOfam files as re-hosted (and refreshed monthly) on the GToTree repo: README, ko_list, profiles.tar.gz and
# date-retrieved.txt
KOFAMSCAN_DL_LINK = "https://github.com/AstrobioMike/gtotree/releases/download/kofamscan-data-latest/kofamscan-data.tar.gz"

EUKCC_DL_LINK = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz"

# the amplicon workflow's taxonomy references: a DECIPHER training set per target region
AMPLICON_REFS_DL_LINK = "https://github.com/AstrobioMike/bit/releases/download/27-May-2026-DECIPHER-dbs/"
AMPLICON_REF_FILES = {"16S": "SILVA_SSU_r138.2_v2.RData", "ITS": "UNITE_v2025.RData", "18S": "PR2_v4_13_species_level.RData"}

# what each database is downloaded from: the sites checked before a download starts (bit.modules.connectivity)
DOWNLOAD_URLS = {
    "tax-db": (TAX_DB_URL,),
    "checkm2-db": (CHECKM2_RECORD, CHECKM2_DL_LINK),
    "gtdbtk-db": (GTDB_MD5_LIST, GTDB_DL_LINK),
    "kofamscan-db": (KOFAMSCAN_DL_LINK,),
    "eukcc-db": (EUKCC_DL_LINK,),
    "amplicon-refs": (AMPLICON_REFS_DL_LINK,),
}


@dataclass(frozen=True)
class RefDb:
    name: str        # as in `bit data get <name>`
    title: str
    subdir: str      # its directory under the reference-database directory
    ready: tuple     # files (in subdir) whose presence means it is set up; the last one written is the last listed
    size: str        # what it takes, for the messages
    used_by: str

    def path(self, root):
        return Path(root) / self.subdir

    def is_ready(self, root):
        return all((self.path(root) / f).exists() for f in self.ready)


REF_DBS = {db.name: db for db in (
    RefDb("tax-db", f"taxonomy reference database ({TAX_DB_NAME})", TAX_DB_NAME,
          ("proteins.dmnd", "protein-lineages.tsv", "genome-kmers.hashes.u64", "genome-kmers.genomes.tsv",
           "genome-lineages.tsv", "TAX_DB_SETUP"),
          "~97 GB to download, ~160 GB once unpacked (~260 GB free needed while unpacking)",
          "bit classify-seqs, the metagenomics workflow, characterize-genomes with --euk, and (where GTDB-Tk doesn't run) "
          "magrec and characterize-genomes"),
    RefDb("checkm2-db", "CheckM2 reference database", "checkm2-ref-db-1.1.0", (CHECKM2_DB_FILENAME, "CHECKM2_DB_SETUP"),
          "~1.7 GB to download, ~3 GB once unpacked", "the metagenomics and magrec workflows, and characterize-genomes"),
    RefDb("gtdbtk-db", "GTDB-Tk reference database (R232)", "GTDB-tk-ref-db-r232", ("GTDBTK_DB_SETUP",),
          "~61 GB to download, plus its unpacked size", "the metagenomics and magrec workflows, and characterize-genomes (Linux only)"),
    RefDb("kofamscan-db", "KOFamScan reference database", "kofamscan-db", ("ko_list", "KO_DB_SETUP"),
          "~1.6 GB to download", "the metagenomics workflow"),
    RefDb("eukcc-db", "EukCC reference database", "eukcc2_db_ver_1.1", ("eukcc_DB_SETUP",),
          "~6.3 GB to download", "characterize-genomes with --euk"),
    RefDb("amplicon-refs", "amplicon taxonomy references (DECIPHER training sets for 16S, 18S and ITS)",
          "DECIPHER-dbs-27-May-2026", tuple(AMPLICON_REF_FILES.values()), "~575 MB to download", "the amplicon workflow"),
)}


# the -d/--ref-dbs-dir setting of the workflows
REF_DBS_DIR_HELP = (f"Directory holding the reference databases (default: the one in the {REF_DBS_VAR} variable, set with "
                    "`bit data locations set`)")


class RefDbsDirNotSet(Exception):
    pass


def ref_dbs_root(given=None):
    """ the reference-database directory: the one given (-d/--ref-dbs-dir), otherwise BIT_REF_DBS_DIR """

    root = given or os.environ.get(REF_DBS_VAR, "")
    if not root:
        raise RefDbsDirNotSet(
            f"No reference-database directory is set: give one with -d/--ref-dbs-dir, or set the {REF_DBS_VAR} "
            "variable with `bit data locations set` (the databases can take a few hundred GB, so somewhere with room)")
    return Path(root).expanduser().resolve()


def missing(root, names):
    """ those of the named databases that aren't set up under root """

    return [name for name in names if not REF_DBS[name].is_ready(root)]


def download_urls(names):
    """ the URLs the named databases are downloaded from """

    return [url for name in names for url in DOWNLOAD_URLS[name]]


def missing_dbs_note(root, names, in_output_dir=False, for_a_run=True):
    """
    The lines saying which databases aren't set up under root and how to get each on another machine, for the
    messages given when they can't be downloaded here (bit.modules.connectivity).
    """

    one = len(names) == 1
    these = "this reference database" if one else "these reference databases"
    where = f"{root}/" if not in_output_dir else \
            f"{root}/ (no reference-database directory is set, so {'it' if one else 'they'} would go in the output directory)"
    first = f"This run needs {these}, which {'isn' if one else 'aren'}'t set up yet in {where}:" if for_a_run else \
            f"{these.capitalize()} {'isn' if one else 'aren'}'t set up yet in {where}:"
    width = max(len(name) for name in names)
    return [
        f"    {first}",
        *[f"        {name:<{width}}   {REF_DBS[name].size}" for name in names],
        "",
        f"    {'It' if one else 'They'} can be downloaded on another machine with internet access and then copied here",
        "    (or into a shared location both can see), using:",
        *[f"        bit data get {name} -d <dir>" for name in names],
        f"    then point {REF_DBS_VAR} (or -d/--ref-dbs-dir) at that directory.",
    ]


def _info(db_dir, source):
    (Path(db_dir) / "DB-INFO.txt").write_text(f"source\t{source}\nretrieved\t{date.today().isoformat()}\n")


def _setup_tax_db(root, db_dir):
    from bit.modules.downloads import fetch, read_url, unpack
    md5_list = db_dir / f"{TAX_DB_NAME}.parts.md5"
    md5_list.write_bytes(read_url(TAX_DB_URL + md5_list.name))
    parts = []
    for line in md5_list.read_text().splitlines():
        if line.strip():
            md5, name = line.split()
            parts.append((TAX_DB_URL + name, str(db_dir / name), md5))
    notes = fetch(parts, TAX_DB_NAME)
    unpack([path for _, path, _ in parts], db_dir, f"Unpacking {TAX_DB_NAME}", ["zstd", "-dc"], strip_components=1)
    for _, path, _ in parts:
        os.remove(path)
    return notes


def _setup_checkm2_db(root, db_dir):
    from bit.modules.downloads import fetch, read_url, unpack
    record = json.loads(read_url(CHECKM2_RECORD))
    md5 = next(f["checksum"].split(":", 1)[1] for f in record["files"] if f["key"] == "checkm2_database.tar.gz")
    archive = db_dir / "checkm2_database.tar.gz"
    notes = fetch([(CHECKM2_DL_LINK, str(archive), md5)], "CheckM2 database")
    # the archive holds CheckM2_database/<the diamond db>, which is put directly in the db dir
    unpack([archive], db_dir, "Unpacking CheckM2 database", ["pigz", "-dc"], strip_components=1)
    archive.unlink()
    _info(db_dir, CHECKM2_DL_LINK)
    return notes


def _setup_gtdbtk_db(root, db_dir):
    from bit.modules.downloads import fetch, read_url, unpack
    md5 = next(line.split()[0] for line in read_url(GTDB_MD5_LIST).decode().splitlines()
               if line.split()[1:] == [GTDB_MD5_ENTRY])
    archive = db_dir / "gtdbtk_data.tar.gz"
    notes = fetch([(GTDB_DL_LINK, str(archive), md5)], "GTDB-Tk database", connections=16)
    unpack([archive], db_dir, "Unpacking GTDB-Tk database", ["pigz", "-dc"], strip_components=1)
    archive.unlink()
    _info(db_dir, GTDB_DL_LINK)
    return notes


def _setup_kofamscan_db(root, db_dir):
    from bit.modules.downloads import DownloadError, fetch, unpack
    archive = db_dir / "kofamscan-data.tar.gz"
    notes = fetch([(KOFAMSCAN_DL_LINK, str(archive), None)], "KOFamScan data")
    unpack([archive], db_dir, "Unpacking KOFamScan data", ["pigz", "-dc"])
    archive.unlink()
    profiles = db_dir / "profiles.tar.gz"
    unpack([profiles], db_dir, "Unpacking KOFamScan profiles", ["pigz", "-dc"])
    profiles.unlink()
    if not (db_dir / "profiles").is_dir():
        raise DownloadError("the KOFamScan data held no profiles")
    retrieved = (db_dir / "date-retrieved.txt").read_text().strip() if (db_dir / "date-retrieved.txt").exists() else "?"
    _info(db_dir, KOFAMSCAN_DL_LINK)
    with open(db_dir / "DB-INFO.txt", "a") as f:
        f.write(f"KOfam files retrieved from KEGG (year,month,day)\t{retrieved}\n")
    return notes


def _setup_eukcc_db(root, db_dir):
    from bit.modules.downloads import fetch, unpack
    # the archive holds eukcc2_db_ver_1.1/, this database's directory
    archive = Path(root) / "eukcc_db.tar.gz"
    notes = fetch([(EUKCC_DL_LINK, str(archive), None)], "EukCC database")
    unpack([archive], root, "Unpacking EukCC database", ["pigz", "-dc"])
    archive.unlink()
    _info(db_dir, EUKCC_DL_LINK)
    return notes


def _setup_amplicon_refs(root, db_dir):
    from bit.modules.downloads import fetch
    return fetch([(AMPLICON_REFS_DL_LINK + name, str(db_dir / name), None) for name in AMPLICON_REF_FILES.values()],
                 "amplicon references")


_SETUP = {
    "tax-db": _setup_tax_db,
    "checkm2-db": _setup_checkm2_db,
    "gtdbtk-db": _setup_gtdbtk_db,
    "kofamscan-db": _setup_kofamscan_db,
    "eukcc-db": _setup_eukcc_db,
    "amplicon-refs": _setup_amplicon_refs,
}


def setup(name, root, force=False, quiet=False, check_internet=True):
    """
    Sets up the named database under root (nothing to do if it already is, unless force: then it is removed and set
    up again). Exits with a message if it fails; running the same command again carries on from where it stopped.
    Before anything is removed or downloaded, the sites it comes from are checked (unless the caller already has).
    """

    from bit.modules.downloads import DownloadError
    from bit.modules.general import color_text, wprint

    db = REF_DBS[name]
    db_dir = db.path(root)

    if db.is_ready(root) and not force:
        if not quiet:
            print(color_text(f"\n    The {db.title} is already set up in {db_dir}/\n", "green"))
        return

    if name == "gtdbtk-db" and not GTDBTK_RUNS_HERE:
        if force and db_dir.exists():
            shutil.rmtree(db_dir)
        wprint(color_text("GTDB-Tk only runs on Linux, so its reference isn't set up on this macOS machine.", "yellow"))
        return

    if check_internet:
        from bit.modules.connectivity import require_internet
        require_internet(download_urls([name]), what="the db", missing_dbs=[name], ref_dbs_dir=root, for_a_run=False)

    if force and db_dir.exists():
        shutil.rmtree(db_dir)

    print(color_text(f"\n    Setting up the {db.title} in {db_dir}/", "yellow"))
    print(f"    ({db.size})\n", flush=True)
    db_dir.mkdir(parents=True, exist_ok=True)
    marker = db_dir / db.ready[-1]
    try:
        notes = _SETUP[name](Path(root), db_dir)
    except (DownloadError, OSError, StopIteration, KeyError, ValueError) as e:
        print()
        wprint(color_text(f"Setting up the {db.title} failed: {e}", "red"))
        wprint("Running the same command again carries on from where it stopped.")
        print()
        sys.exit(1)
    if not marker.exists():
        marker.touch()
    if not db.is_ready(root):
        wprint(color_text(f"The {db.title} was downloaded but doesn't hold what it should; please report this.", "red"))
        sys.exit(1)
    for note in notes:
        print(f"      {note}")
    print(color_text(f"\n    The {db.title} is set up.\n", "green"))


def set_up_missing(root, names, check_internet=True):
    """
    Sets up those of the named databases that aren't yet (for the workflows, before they start). The sites of all of
    them are checked first, so a long download isn't started when a later one couldn't be (unless the caller already
    has checked).
    """

    to_set_up = missing(root, names)
    if to_set_up and check_internet:
        from bit.modules.connectivity import require_internet
        require_internet(download_urls(to_set_up), what="the db", missing_dbs=to_set_up, ref_dbs_dir=root)
    for name in to_set_up:
        setup(name, root, check_internet=False)
