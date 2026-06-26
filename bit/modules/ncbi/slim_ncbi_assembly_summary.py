"""
Combine the NCBI GenBank + RefSeq assembly_summary tables into a single slim
TSV holding only the columns bit uses, with a clean (no leading '#') header row.

This is the projection used both by:
  - the scheduled GitHub Action that rebuilds and re-hosts the slim table weekly
    (run as a module / via build_slim_assembly_summary), and
  - any local rebuild fallback.

NCBI's assembly_summary_*.txt format:
  - line 1 is a '##' provenance comment
  - line 2 is the real header, beginning with '#assembly_accession\t...'
  - remaining lines are tab-delimited data rows in a fixed 38-column layout

The eight columns bit reads (parse_ncbi_assembly_summary.py and
gen_mg/taxonomy.py) are selected BY NAME here, so the stored file can be slimmed
and even reordered upstream without breaking positional reads downstream -- the
readers parse this file's header into a name->index map.

KEPT_COLUMNS is the single source of truth for the slim layout and column order.
"""

import gzip
import os

# the columns bit actually reads from the assembly summary, in the order they
# are written to the slim file. Keep in sync with the readers, which look these
# up by name from the header row.
KEPT_COLUMNS = [
    "assembly_accession",
    "taxid",
    "organism_name",
    "infraspecific_name",
    "version_status",
    "assembly_level",
    "asm_name",
    "ftp_path",
]

# NCBI's real header line starts with this token (with the leading '#')
_NCBI_HEADER_LEAD = "#assembly_accession"


def _open_maybe_gzip(path, mode="rt"):
    """open a path as text, transparently handling .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_ncbi_header(line):
    """
    turn an NCBI assembly_summary header line into a list of clean column names.
    Strips a single leading '#' from the first field if present. Raises
    ValueError if the line doesn't look like the expected header.
    """
    stripped = line.rstrip("\n")
    if not stripped.startswith("#"):
        raise ValueError("expected an NCBI assembly_summary header line "
                         "(beginning with '#assembly_accession')")
    # drop only the leading '#', keep the rest of that first field
    stripped = stripped[1:]
    names = stripped.split("\t")
    if not names or names[0] != "assembly_accession":
        raise ValueError(
            "unexpected NCBI assembly_summary header; first column is "
            f"'{names[0] if names else ''}', expected 'assembly_accession'")
    return names


def _read_summary(path):
    """
    yield (header_names, data_line_iterator) for one NCBI assembly_summary file.
    Skips the '##' provenance line, parses the '#...'' header, and yields the
    remaining raw data lines (newline-stripped, non-empty).
    """
    fh = _open_maybe_gzip(path, "rt")
    header = None
    data_lines = []
    for line in fh:
        if header is None:
            if line.startswith("##"):
                continue                    # provenance comment -> skip
            if line.startswith("#"):
                header = parse_ncbi_header(line)
                continue
            # a file with no header at all -> not the expected format
            raise ValueError(f"no header found in {path} before data rows")
        s = line.rstrip("\n")
        if s:
            data_lines.append(s)
    fh.close()
    if header is None:
        raise ValueError(f"no header found in {path}")
    return header, data_lines


def _project_rows(header, data_lines, keep_idx):
    """
    yield slimmed rows (lists of values) selecting keep_idx positions from each
    data line. Rows with too few fields to cover the needed indices are skipped
    (defensive against malformed lines).
    """
    max_needed = max(keep_idx)
    out = []
    for s in data_lines:
        fields = s.split("\t")
        if len(fields) <= max_needed:
            continue
        out.append([fields[i] for i in keep_idx])
    return out


def build_slim_assembly_summary(genbank_path, refseq_path, out_path,
                                keep_columns=None):
    """
    combine the GenBank and RefSeq assembly_summary files into a single slim TSV
    at out_path, keeping only keep_columns (default KEPT_COLUMNS) BY NAME, with a
    clean header row. Returns the number of data rows written.

    Both inputs must share the kept column names (they do -- GenBank and RefSeq
    use the same assembly_summary schema). Columns are located by name in each
    file's own header, so the two files are allowed to differ in column order.
    """
    keep_columns = list(keep_columns) if keep_columns else list(KEPT_COLUMNS)

    gb_header, gb_data = _read_summary(genbank_path)
    rs_header, rs_data = _read_summary(refseq_path)

    def _indices(header, source):
        missing = [c for c in keep_columns if c not in header]
        if missing:
            raise ValueError(
                f"{source} assembly_summary is missing expected column(s): "
                f"{', '.join(missing)}")
        return [header.index(c) for c in keep_columns]

    gb_idx = _indices(gb_header, "GenBank")
    rs_idx = _indices(rs_header, "RefSeq")

    rows = _project_rows(gb_header, gb_data, gb_idx)
    rows += _project_rows(rs_header, rs_data, rs_idx)

    tmp_path = str(out_path) + ".tmp"
    with open(tmp_path, "w") as out:
        out.write("\t".join(keep_columns) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    os.replace(tmp_path, out_path)        # atomic swap into place

    return len(rows)


GENBANK_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

# names of the two files packaged inside the distributed tarball (kept flat at
# the archive root, mirroring the GTDB slim bundle).
TABLE_FILENAME = "ncbi-assembly-info.tsv"
DATE_FILENAME = "date-retrieved.txt"


def write_date_retrieved(path, when=None):
    """
    write a date-retrieved.txt in bit's stored format (YYYY,MM,DD). Defaults to
    today (the build date when called from the Action). gen_metagenome reads this
    back to report how fresh the assembly data is.
    """
    from datetime import date as _date
    when = when or _date.today()
    with open(path, "w") as fh:
        fh.write(when.strftime("%Y,%m,%d") + "\n")


def main(argv=None): # pragma: no cover
    """
    download the GenBank + RefSeq assembly summaries, build the slim combined
    table, stamp a build-date date-retrieved.txt, and package both into a single
    .tar.gz. Intended for the scheduled GitHub Action and for a manual rebuild.
    Network download uses urllib so this has no third-party deps beyond the
    standard library.

    Output: <out_dir>/<archive-name> (default ncbi-assembly-info.tar.gz), holding
    ncbi-assembly-info.tsv and date-retrieved.txt at the archive root.
    """
    import argparse
    import tarfile
    import urllib.request

    ap = argparse.ArgumentParser(description="Build bit's slim NCBI assembly-info bundle.")
    ap.add_argument("-o", "--out-dir", default=".", help="output directory")
    ap.add_argument("--archive-name", default="ncbi-assembly-info.tar.gz",
                    help="name of the produced .tar.gz")
    ap.add_argument("--genbank-url", default=GENBANK_URL)
    ap.add_argument("--refseq-url", default=REFSEQ_URL)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    gb_path = os.path.join(args.out_dir, "assembly_summary_genbank.txt")
    rs_path = os.path.join(args.out_dir, "assembly_summary_refseq.txt")
    table_path = os.path.join(args.out_dir, TABLE_FILENAME)
    date_path = os.path.join(args.out_dir, DATE_FILENAME)
    archive_path = os.path.join(args.out_dir, args.archive_name)

    print(f"Downloading GenBank summary -> {gb_path}", flush=True)
    urllib.request.urlretrieve(args.genbank_url, gb_path)
    print(f"Downloading RefSeq summary -> {rs_path}", flush=True)
    urllib.request.urlretrieve(args.refseq_url, rs_path)

    print("Building slim combined table...", flush=True)
    n = build_slim_assembly_summary(gb_path, rs_path, table_path)
    print(f"  wrote {n:,} rows, {len(KEPT_COLUMNS)} columns", flush=True)

    write_date_retrieved(date_path)
    print(f"Stamped build date -> {date_path}", flush=True)

    print(f"Packaging -> {archive_path}", flush=True)
    with tarfile.open(archive_path, "w:gz") as tar:
        # arcname keeps the files flat at the archive root (no out_dir prefix)
        tar.add(table_path, arcname=TABLE_FILENAME)
        tar.add(date_path, arcname=DATE_FILENAME)

    # tidy intermediates; leave only the archive
    for p in (gb_path, rs_path, table_path, date_path):
        if os.path.exists(p):
            os.remove(p)

    print("Done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
