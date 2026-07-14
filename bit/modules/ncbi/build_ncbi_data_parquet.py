"""
Build the slim NCBI assembly-info table WITH full 7-rank lineages, as a single
sorted, compressed Parquet file.

This runs on a GitHub Actions runner, not on user machines.

Output: ncbi-data.parquet
  assembly_accession, taxid, organism_name, infraspecific_name, version_status,
  assembly_level, asm_name, ftp_path,
  domain, phylum, class, order, family, genus, species,
  domain_taxid, ... , species_taxid

The table is sorted by lineage before writing so that (a) Parquet row-group
statistics allow whole-row-group skipping on rank predicates, and (b) the
dictionary-encoded lineage columns compress far better (identical values become
contiguous).
"""

import argparse
import gzip
import os
import sys
import tarfile
import urllib.request

import pyarrow as pa # type: ignore
import pyarrow.parquet as pq # type: ignore

from bit.modules.taxonomy.tax_ranks import (RANKS, NA, accession_core,
                                            LINEAGE_NAME_COLUMNS,
                                            LINEAGE_TAXID_COLUMNS)


GENBANK_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
CHECKM_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/CheckM_report_prokaryotes.txt"

ASSEMBLY_COLUMNS = [
    "assembly_accession",
    "taxid",
    "organism_name",
    "infraspecific_name",
    "version_status",
    "assembly_level",
    "refseq_category",
    "asm_name",
    "ftp_path",
    "genome_size",
    "genome_size_ungapped",
    "contig_count",
]

CHECKM_COLUMNS = ["checkm_completeness", "checkm_contamination"]

RANK_ALIASES = {
    "domain": "domain",
    "superkingdom": "domain",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}

OUT_COLUMNS = (ASSEMBLY_COLUMNS + CHECKM_COLUMNS
               + LINEAGE_NAME_COLUMNS + LINEAGE_TAXID_COLUMNS)

PARQUET_FILENAME = "ncbi-data.parquet"
DATE_FILENAME = "date-retrieved.txt"

_NCBI_HEADER_LEAD = "#assembly_accession"

_EXCLUDED_VERSION_STATUS = {"suppressed", "replaced"}


# ---------------------------------------------------------------------------
# taxdump parsing
# ---------------------------------------------------------------------------

def parse_nodes_dmp(path):
    """
    nodes.dmp -> {taxid: (parent_taxid, rank)}

    Format is pipe-delimited with surrounding tabs:
        taxid\t|\tparent\t|\trank\t|\t...
    """
    nodes = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.split("\t|\t")
            if len(fields) < 3:
                continue
            taxid = fields[0].strip()
            parent = fields[1].strip()
            rank = fields[2].strip().lower()
            if taxid:
                nodes[taxid] = (parent, rank)
    return nodes


def parse_names_dmp(path):
    """
    names.dmp -> {taxid: scientific_name}

    Only 'scientific name' rows are kept; synonyms/common names are ignored.
    """
    names = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            # trailing "\t|\n" then pipe-delimited
            fields = line.rstrip("\n").rstrip("\t|").split("\t|\t")
            if len(fields) < 4:
                continue
            if fields[3].strip() != "scientific name":
                continue
            taxid = fields[0].strip()
            name = fields[1].strip()
            if taxid:
                names[taxid] = name
    return names


def resolve_lineage(taxid, nodes, names):
    """
    Walk the parent chain from `taxid` to the root, collecting the 7 standard
    ranks. Returns (names_tuple, taxids_tuple), each length 7, coarse -> fine,
    with "NA" filled for any rank absent from the lineage (very common in NCBI:
    candidate phyla, unclassified/environmental clades, etc.).

    Unknown/dangling taxids resolve to all-NA rather than raising, so one bad
    row can never fail the whole build.
    """
    found_names = {r: NA for r in RANKS}
    found_taxids = {r: NA for r in RANKS}

    seen = set()
    current = str(taxid).strip()

    while current and current not in seen and current in nodes:
        seen.add(current)  # cycle guard; taxdump shouldn't have any, but be safe
        parent, rank = nodes[current]

        slot = RANK_ALIASES.get(rank)
        # first (deepest) occurrence of a rank wins; walking upward means we
        # would otherwise overwrite species with a parent that shares a rank
        if slot is not None and found_taxids[slot] == NA:
            found_taxids[slot] = current
            found_names[slot] = names.get(current, NA)

        if current == parent:      # root (taxid 1 is its own parent)
            break
        current = parent

    return (
        tuple(found_names[r] for r in RANKS),
        tuple(found_taxids[r] for r in RANKS),
    )


def resolve_lineages_for(taxids, nodes, names):
    """{taxid: (names_tuple, taxids_tuple)} for an iterable of unique taxids."""
    return {t: resolve_lineage(t, nodes, names) for t in taxids}


# ---------------------------------------------------------------------------
# assembly_summary parsing
# ---------------------------------------------------------------------------

def _open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def parse_ncbi_header(line):
    """'#assembly_accession\\t...' -> clean column-name list."""
    stripped = line.rstrip("\n")
    if not stripped.startswith("#"):
        raise ValueError("expected an NCBI assembly_summary header line "
                         "(beginning with '#assembly_accession')")
    names = stripped[1:].split("\t")
    if not names or names[0] != "assembly_accession":
        raise ValueError(
            "unexpected NCBI assembly_summary header; first column is "
            f"'{names[0] if names else ''}', expected 'assembly_accession'")
    return names


def iter_assembly_rows(path, keep_columns=None):
    """
    Yield slimmed rows (lists of str, ordered as keep_columns) from one
    assembly_summary file. Columns are located BY NAME in that file's own header,
    so GenBank and RefSeq are allowed to differ in column order.
    """
    keep_columns = list(keep_columns) if keep_columns else list(ASSEMBLY_COLUMNS)

    with _open_maybe_gzip(path) as fh:
        header = None
        for line in fh:
            if header is None:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = parse_ncbi_header(line)
                    missing = [c for c in keep_columns if c not in header]
                    if missing:
                        raise ValueError(
                            f"{path} is missing expected column(s): {', '.join(missing)}")
                    idx = [header.index(c) for c in keep_columns]
                    vs_i = header.index("version_status") if "version_status" in header else None
                    max_needed = max(idx + ([vs_i] if vs_i is not None else []))
                    continue
                raise ValueError(f"no header found in {path} before data rows")

            s = line.rstrip("\n")
            if not s:
                continue
            fields = s.split("\t")
            if len(fields) <= max_needed:
                continue
            if vs_i is not None and fields[vs_i].strip().lower() in _EXCLUDED_VERSION_STATUS:
                continue

            yield [fields[i] for i in idx]

        if header is None:
            raise ValueError(f"no header found in {path}")


# ---------------------------------------------------------------------------
# CheckM report
# ---------------------------------------------------------------------------

_CHECKM_ACC_KEYS = ("genbank_accession", "refseq_accession", "assembly_accession",
                    "accession")
_CHECKM_COMP_KEYS = ("checkm_completeness", "completeness")
_CHECKM_CONT_KEYS = ("checkm_contamination", "contamination")


def _norm_header_name(name):
    """'#genbank-accession' / 'CheckM Completeness' -> 'genbank_accession' / 'checkm_completeness'"""
    return (str(name).strip().lstrip("#").strip().lower()
            .replace("-", "_").replace(" ", "_"))


def _find_col(header, candidates):
    normed = [_norm_header_name(h) for h in header]
    for cand in candidates:
        key = _norm_header_name(cand)
        if key in normed:
            return normed.index(key)
    return None


def _find_all_cols(header, candidates):
    """All matching column indices, in candidate order (for the 2 accession cols)."""
    normed = [_norm_header_name(h) for h in header]
    out = []
    for cand in candidates:
        key = _norm_header_name(cand)
        if key in normed:
            i = normed.index(key)
            if i not in out:
                out.append(i)
    return out


def parse_checkm_report(path):
    out = {}
    try:
        with _open_maybe_gzip(path) as fh:
            header = None
            for line in fh:
                if header is None:
                    if line.startswith("##") or not line.strip():
                        continue
                    header = line.rstrip("\n").split("\t")
                    acc_idx = _find_all_cols(header, _CHECKM_ACC_KEYS)
                    ci = _find_col(header, _CHECKM_COMP_KEYS)
                    xi = _find_col(header, _CHECKM_CONT_KEYS)
                    if not acc_idx or ci is None or xi is None:
                        print(f"    WARNING: CheckM report header not understood "
                              f"({header[:6]}...); checkm will be NA", flush=True)
                        return {}
                    need = max(acc_idx + [ci, xi])
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= need:
                    continue
                comp = fields[ci].strip() or NA
                cont = fields[xi].strip() or NA
                for ai in acc_idx:
                    core = accession_core(fields[ai])
                    if core:
                        out[core] = (comp, cont)
    except OSError as e:
        print(f"    WARNING: could not read CheckM report ({e}); checkm will be NA",
              flush=True)
        return {}
    return out




# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

def build_table(genbank_path, refseq_path, nodes_path, names_path,
                checkm_path=None):
    """
    Join the GenBank + RefSeq assembly summaries against the taxdump (and, if
    given, NCBI's bulk CheckM report) and return a lineage-sorted pyarrow.Table.
    """
    rows = []
    for src in (genbank_path, refseq_path):
        rows.extend(iter_assembly_rows(src))

    if not rows:
        raise ValueError("no assembly rows parsed; refusing to build an empty table")

    taxid_pos = ASSEMBLY_COLUMNS.index("taxid")
    unique_taxids = {r[taxid_pos] for r in rows}

    print(f"    parsed {len(rows):,} assembly rows "
          f"({len(unique_taxids):,} unique taxids)", flush=True)

    nodes = parse_nodes_dmp(nodes_path)
    names = parse_names_dmp(names_path)
    print(f"    taxdump: {len(nodes):,} nodes, {len(names):,} scientific names", flush=True)

    lineages = resolve_lineages_for(unique_taxids, nodes, names)
    del nodes, names

    unresolved = sum(1 for t in unique_taxids if lineages[t][1][0] == NA)
    print(f"    resolved lineages for {len(lineages):,} taxids "
          f"({unresolved:,} with no domain)", flush=True)

    checkm = parse_checkm_report(checkm_path) if checkm_path else {}
    if checkm_path:
        print(f"    checkm: {len(checkm):,} assemblies with completeness/contamination",
              flush=True)

    acc_pos = ASSEMBLY_COLUMNS.index("assembly_accession")

    # build columnar arrays directly (avoids a second full row-wise copy)
    cols = {c: [] for c in OUT_COLUMNS}
    n_checkm_hit = 0
    for r in rows:
        for c, v in zip(ASSEMBLY_COLUMNS, r):
            cols[c].append(v)

        comp, cont = checkm.get(accession_core(r[acc_pos]), (NA, NA))
        if comp != NA:
            n_checkm_hit += 1
        cols["checkm_completeness"].append(comp)
        cols["checkm_contamination"].append(cont)

        lin_names, lin_taxids = lineages[r[taxid_pos]]
        for c, v in zip(LINEAGE_NAME_COLUMNS, lin_names):
            cols[c].append(v)
        for c, v in zip(LINEAGE_TAXID_COLUMNS, lin_taxids):
            cols[c].append(v)
    del rows

    if checkm:
        print(f"    checkm joined onto {n_checkm_hit:,} rows "
              f"({n_checkm_hit/len(cols['assembly_accession']):.1%}; "
              f"prokaryotes only by design)", flush=True)

    table = pa.table({c: pa.array(cols[c], type=pa.string()) for c in OUT_COLUMNS})
    del cols

    # sort by lineage: makes identical lineage values contiguous, which both
    # compresses far better and lets Parquet row-group stats skip on rank
    # predicates
    table = table.sort_by([(r, "ascending") for r in RANKS]
                          + [("assembly_accession", "ascending")])
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


def extract_taxdump(tarball_path, dest_dir):
    """Extract only names.dmp and nodes.dmp from taxdump.tar.gz."""
    wanted = {"names.dmp", "nodes.dmp"}
    out = {}
    with tarfile.open(tarball_path, "r:gz") as tar:
        for member in tar.getmembers():
            base = os.path.basename(member.name)
            if base in wanted and member.isfile():
                src = tar.extractfile(member)
                dest = os.path.join(dest_dir, base)
                with src, open(dest, "wb") as fh:
                    fh.write(src.read())
                out[base] = dest
    missing = wanted - set(out)
    if missing:
        raise ValueError(f"taxdump archive missing: {', '.join(sorted(missing))}")
    return out["nodes.dmp"], out["names.dmp"]


def write_date_retrieved(path, when=None):
    from datetime import date as _date
    when = when or _date.today()
    with open(path, "w") as fh:
        fh.write(when.strftime("%Y,%m,%d") + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build the lineage-joined NCBI assembly-info Parquet asset.")
    ap.add_argument("--out-dir", default="build")
    ap.add_argument("--work-dir", default="work")
    ap.add_argument("--genbank", help="local assembly_summary_genbank.txt (skips download)")
    ap.add_argument("--refseq", help="local assembly_summary_refseq.txt (skips download)")
    ap.add_argument("--taxdump", help="local taxdump.tar.gz (skips download)")
    ap.add_argument("--checkm", help="local CheckM_report_prokaryotes.txt (skips download)")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    genbank = args.genbank or download(
        GENBANK_URL, os.path.join(args.work_dir, "assembly_summary_genbank.txt"))
    refseq = args.refseq or download(
        REFSEQ_URL, os.path.join(args.work_dir, "assembly_summary_refseq.txt"))
    taxdump = args.taxdump or download(
        TAXDUMP_URL, os.path.join(args.work_dir, "taxdump.tar.gz"))

    # fail SOFT: a checkm hiccup must not take down the whole weekly asset
    checkm = args.checkm
    if not checkm:
        try:
            checkm = download(CHECKM_URL,
                              os.path.join(args.work_dir, "CheckM_report_prokaryotes.txt"))
        except Exception as e:
            print(f"    WARNING: could not fetch the CheckM report ({e}); "
                  f"checkm columns will be NA", flush=True)
            checkm = None

    nodes_path, names_path = extract_taxdump(taxdump, args.work_dir)

    table = build_table(genbank, refseq, nodes_path, names_path, checkm_path=checkm)

    out_path = os.path.join(args.out_dir, PARQUET_FILENAME)
    size = write_parquet(table, out_path)
    print(f"    wrote {out_path}  ({size/1048576:.1f} MB, {table.num_rows:,} rows)")

    write_date_retrieved(os.path.join(args.out_dir, DATE_FILENAME))
    return 0


if __name__ == "__main__":
    sys.exit(main())
