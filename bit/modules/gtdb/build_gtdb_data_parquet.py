"""
Build the slim GTDB metadata table as a single sorted, zstd-compressed Parquet
file, on the same lineage schema as the NCBI asset.

This runs on a GitHub Actions runner (see .github/workflows/check-gtdb-metadata.yaml),
NOT on user machines.
"""

import argparse
import gzip
import os
import shutil
import sys
import urllib.request
import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore
from bit.modules.taxonomy.tax_ranks import RANKS, NA


GTDB_BASE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"
ARC_FILENAME = "ar53_metadata.tsv.gz"
BAC_FILENAME = "bac120_metadata.tsv.gz"
VERSION_URL = f"{GTDB_BASE_URL}/VERSION.txt"

PARQUET_FILENAME = "gtdb-data.parquet"
VERSION_FILENAME = "GTDB-version-info.txt"

GTDB_KEPT_COLUMNS = [
    "accession", "ncbi_genbank_assembly_accession", "ncbi_taxid",
    "gtdb_representative", "ncbi_refseq_category",
    "checkm2_completeness", "checkm2_contamination",
    "checkm_completeness", "checkm_contamination",
    "genome_size", "contig_count", "gc_count", "gc_percentage",
    "ambiguous_bases", "coding_bases", "coding_density",
]

OUT_COLUMNS = GTDB_KEPT_COLUMNS[:5] + list(RANKS) + GTDB_KEPT_COLUMNS[5:]

RANK_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]


# ---------------------------------------------------------------------------
# taxonomy string splitting
# ---------------------------------------------------------------------------

def split_gtdb_taxonomy(taxonomy_string):
    parts = [p.strip() for p in str(taxonomy_string).split(";")]
    if len(parts) != 7:
        raise ValueError(
            f"expected 7 semicolon-delimited ranks in GTDB taxonomy, got "
            f"{len(parts)}: {taxonomy_string!r}")

    out = []
    for part, prefix in zip(parts, RANK_PREFIXES):
        if not part.startswith(prefix):
            raise ValueError(
                f"expected GTDB rank prefix '{prefix}' in {part!r} "
                f"(full lineage: {taxonomy_string!r})")
        name = part[len(prefix):].strip()
        out.append(name if name else NA)
    return out


# ---------------------------------------------------------------------------
# upstream metadata parsing
# ---------------------------------------------------------------------------

def _open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_gtdb_rows(path, keep_columns=None):
    """
    Yield rows as dicts of {kept_column: value} plus the 7 split rank columns,
    from one upstream GTDB metadata table.

    Columns are located BY NAME in the file's own header. A kept column absent
    from a given release is filled with NA rather than failing the build (GTDB
    has churned column names across releases -- e.g. checkm -> checkm2) -- but
    the two columns we cannot proceed without (`accession`, `gtdb_taxonomy`)
    raise.
    """
    keep_columns = list(keep_columns) if keep_columns else list(GTDB_KEPT_COLUMNS)

    with _open_maybe_gzip(path) as fh:
        header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")

        if header and header[0] != "accession":
            header[0] = "accession"

        idx = {name: i for i, name in enumerate(header)}

        for required in ("accession", "gtdb_taxonomy"):
            if required not in idx:
                raise ValueError(
                    f"{path} is missing the required column '{required}'")

        tax_i = idx["gtdb_taxonomy"]
        max_needed = max([tax_i] + [idx[c] for c in keep_columns if c in idx])

        for line in fh:
            s = line.rstrip("\n")
            if not s:
                continue
            fields = s.split("\t")
            if len(fields) <= max_needed:
                continue

            row = {}
            for c in keep_columns:
                row[c] = fields[idx[c]].strip() if c in idx else NA
                if row[c] == "":
                    row[c] = NA

            ranks = split_gtdb_taxonomy(fields[tax_i])
            for rank_name, value in zip(RANKS, ranks):
                row[rank_name] = value

            yield row


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

def build_table(arc_path, bac_path):
    """Combine the archaeal + bacterial metadata into a lineage-sorted table."""
    cols = {c: [] for c in OUT_COLUMNS}
    n = 0
    for src in (arc_path, bac_path):
        for row in iter_gtdb_rows(src):
            for c in OUT_COLUMNS:
                cols[c].append(row.get(c, NA))
            n += 1

    if n == 0:
        raise ValueError("no GTDB rows parsed; refusing to build an empty table")

    print(f"    parsed {n:,} GTDB rows", flush=True)

    table = pa.table({c: pa.array(cols[c], type=pa.string()) for c in OUT_COLUMNS})
    del cols

    # sort by lineage: identical values become contiguous -> far better
    # dictionary+zstd compression, and Parquet row-group stats can then skip
    # whole row groups on a rank predicate
    table = table.sort_by([(r, "ascending") for r in RANKS]
                          + [("accession", "ascending")])
    return table


def write_parquet(table, out_path, compression_level=9, row_group_size=256_000):
    pq.write_table(
        table,
        out_path,
        compression="zstd",
        compression_level=compression_level,
        use_dictionary=True,
        write_statistics=True,
        row_group_size=row_group_size,
    )
    return os.path.getsize(out_path)


def download(url, dest):
    print(f"    downloading {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def resolve_version_file(out_dir, work_dir, version_file=None):
    """
    Produce the GTDB version file that gets published alongside the Parquet.

    This file is load-bearing: next week's workflow compares upstream's VERSION.txt
    against the published one to decide whether to rebuild.
    """
    dest = os.path.join(out_dir, VERSION_FILENAME)

    if version_file:
        shutil.copyfile(version_file, dest)
        print(f"    version file taken from {version_file}", flush=True)
        return dest

    download(VERSION_URL, dest)
    return dest


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build the GTDB metadata Parquet asset.")
    ap.add_argument("-o", "--out-dir", default="build")
    ap.add_argument("--work-dir", default="work")
    ap.add_argument("--arc", help="local ar53_metadata.tsv.gz (skips download)")
    ap.add_argument("--bac", help="local bac120_metadata.tsv.gz (skips download)")
    ap.add_argument("--version-file",
                    help="local VERSION.txt (skips download; the workflow passes the "
                         "one it already fetched)")
    ap.add_argument("--no-version-file", action="store_true",
                    help="don't produce a version file at all (local/dev builds)")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    if not args.no_version_file:
        resolve_version_file(args.out_dir, args.work_dir, args.version_file)

    arc = args.arc or download(f"{GTDB_BASE_URL}/{ARC_FILENAME}",
                               os.path.join(args.work_dir, ARC_FILENAME))
    bac = args.bac or download(f"{GTDB_BASE_URL}/{BAC_FILENAME}",
                               os.path.join(args.work_dir, BAC_FILENAME))

    table = build_table(arc, bac)

    out_path = os.path.join(args.out_dir, PARQUET_FILENAME)
    size = write_parquet(table, out_path)
    print(f"    wrote {out_path}  ({size/1048576:.1f} MB, {table.num_rows:,} rows)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
