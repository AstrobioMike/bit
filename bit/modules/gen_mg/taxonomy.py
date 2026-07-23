"""
gen-mg taxonomy-resolution

Fills both taxonomy rank sets (gtdb_* and ncbi_*) on the merged genome table
"""
import os
import pandas as pd # type: ignore
from bit.modules.gen_mg.selection import RANKS, GTDB_RANKS, NCBI_RANKS


def _root_acc(acc):
    """
    strip version suffix for tolerant matching (GCA_000001.2 -> GCA_000001)
    """
    return str(acc).split(".")[0]


def _clean_gtdb_acc(acc):
    """
    strip GTDB's RS_/GB_ prefix from the col-0 accession, leaving e.g.
    'RS_GCF_000005845.2' -> 'GCF_000005845.2'
    """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    return s


def _norm_key(acc):
    """
    the single canonical lookup key for an accession (de-prefix + de-version)
    """
    return _root_acc(_clean_gtdb_acc(acc))


# ---------------------------------------------------------------- GTDB ranks --

def gtdb_taxid_map(gtdb_tab):
    """
    accession(rootless, both GCA & GCF forms) -> ncbi_taxid (str)
    """
    if ("ncbi_genbank_assembly_accession" not in gtdb_tab.columns
            and "accession" not in gtdb_tab.columns):
        return {}
    taxid_col = next((c for c in ("ncbi_taxid", "ncbi_species_taxid") if c in gtdb_tab.columns), None)
    if taxid_col is None:
        return {}

    gca = (gtdb_tab["ncbi_genbank_assembly_accession"].values
           if "ncbi_genbank_assembly_accession" in gtdb_tab.columns else None)
    raw = (gtdb_tab["accession"].values
           if "accession" in gtdb_tab.columns else None)
    tids = gtdb_tab[taxid_col].values

    out = {}
    n = len(tids)
    for i in range(n):
        tid = tids[i]
        if tid is None or (not isinstance(tid, str) and pd.isna(tid)):
            continue
        tid = str(tid).strip()
        if gca is not None:
            k = _norm_key(gca[i])
            if k and k.lower() != "na":
                out[k] = tid
        if raw is not None:
            k = _norm_key(_clean_gtdb_acc(raw[i]))
            if k and k.lower() != "na":
                out[k] = tid
    return out


def gtdb_lineage_map(gtdb_tab):
    """
    accession(rootless, both GCA & GCF forms) -> dict(plain_rank -> gtdb value)
    """
    have_acc = ("ncbi_genbank_assembly_accession" in gtdb_tab.columns
                or "accession" in gtdb_tab.columns)
    if not have_acc or any(c not in gtdb_tab.columns for c in RANKS):
        return {}

    gca = (gtdb_tab["ncbi_genbank_assembly_accession"].values
           if "ncbi_genbank_assembly_accession" in gtdb_tab.columns else None)
    raw = (gtdb_tab["accession"].values
           if "accession" in gtdb_tab.columns else None)
    rank_arrays = {r: gtdb_tab[r].values for r in RANKS}

    out = {}
    n = len(next(iter(rank_arrays.values())))
    for i in range(n):
        lineage = {}
        for r in RANKS:
            v = rank_arrays[r][i]
            lineage[r] = v if (v is not None and not (not isinstance(v, str) and pd.isna(v))) else "NA"
        if gca is not None:
            k = _norm_key(gca[i])
            if k and k.lower() != "na":
                out[k] = lineage
        if raw is not None:
            k = _norm_key(_clean_gtdb_acc(raw[i]))
            if k and k.lower() != "na":
                out[k] = lineage
    return out


def fill_gtdb_taxonomy(merged, gtdb_tab):
    """
    fill gtdb_* columns for any row whose accession is in the GTDB table and
    whose gtdb_* are still NA. Returns the merged frame (modified copy)
    """
    lineage = gtdb_lineage_map(gtdb_tab)
    if not lineage:
        return merged
    out = merged.copy()
    for i, acc in enumerate(out["accession"]):
        # only fill if not already populated (GTDB-selected rows are done)
        if str(out.iloc[i]["gtdb_domain"]) != "NA":
            continue
        hit = lineage.get(_norm_key(acc))
        if hit:
            for gcol, plain in zip(GTDB_RANKS, RANKS):
                out.at[out.index[i], gcol] = hit[plain]
    return out


# ---------------------------------------------------------------- NCBI ranks --

def assembly_info_taxid_map(accessions, assembly_info_path):
    """
    accession(rootless) -> taxid from the NCBI assembly-info Parquet table.
    Only the accession and taxid columns are read, and only rows matching the
    wanted accessions are kept
    """
    wanted = {_norm_key(a) for a in accessions}
    out = {}
    if not assembly_info_path or not os.path.exists(assembly_info_path):
        return out

    import pyarrow.parquet as pq # type: ignore
    tbl = pq.read_table(assembly_info_path, columns=["assembly_accession", "taxid"])
    accs = tbl.column("assembly_accession")
    taxids = tbl.column("taxid")
    for acc, taxid in zip(accs, taxids):
        root = _norm_key(acc.as_py())
        if root in wanted:
            t = (taxid.as_py() or "").strip()
            if t:
                out[root] = t
    return out


_NCBI_CORES_CACHE = {}


def _acc_digits(acc):
    """
    numeric core of an assembly accession, ignoring GCA/GCF prefix and
    version, so GCA_000001.2 and GCF_000001.1 share the key '000001'. GCA/GCF
    twins share their numeric portion
    """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    s = s.split(".")[0]              # drop version
    if "_" in s:
        s = s.split("_", 1)[1]       # drop GCA_/GCF_ prefix -> digits
    return s


def _ncbi_accession_cores(assembly_info_path):
    """
    Arrow array of the numeric cores of every accession in the NCBI assembly-summary
    Parquet. Vectorised equivalent of mapping _acc_digits over the column: strip any
    RS_/GB_ prefix, drop the version suffix, drop the GCA_/GCF_ prefix.
    """
    import pyarrow as pa # type: ignore
    import pyarrow.parquet as pq # type: ignore
    import pyarrow.compute as pc # type: ignore

    key = (assembly_info_path, os.path.getmtime(assembly_info_path))
    hit = _NCBI_CORES_CACHE.get(key)
    if hit is not None:
        return hit

    col = pq.read_table(assembly_info_path, columns=["assembly_accession"]).column(0)
    if isinstance(col, pa.ChunkedArray):
        col = col.combine_chunks()
    cores = pc.replace_substring_regex(col, r'^(RS_|GB_)', '')
    cores = pc.replace_substring_regex(cores, r'\..*$', '')
    cores = pc.replace_substring_regex(cores, r'^[^_]*_', '')

    _NCBI_CORES_CACHE.clear()   # only ever need the current asset
    _NCBI_CORES_CACHE[key] = cores
    return cores


def dead_accession_cores(candidate_accessions, assembly_info_path):
    """
    The numeric cores of `candidate_accessions` that are NOT present in the NCBI
    assembly-summary (i.e., suppressed/removed).

    Intended to be computed ONCE for the whole candidate pool and handed to the
    suppression backfill as `exclude_accessions`, rather than re-screening each round.
    """
    if not assembly_info_path or not os.path.exists(assembly_info_path):
        return set()

    import pyarrow as pa # type: ignore
    import pyarrow.compute as pc # type: ignore

    ncbi_cores = _ncbi_accession_cores(assembly_info_path)
    cand = pa.array([_acc_digits(a) for a in candidate_accessions])
    dead_mask = pc.invert(pc.is_in(cand, value_set=ncbi_cores))
    return set(pc.filter(cand, dead_mask).to_pylist())


def present_accessions(accessions, assembly_info_path):
    """
    subset of `accessions` whose assembly is present in the NCBI
    assembly-summary file. NCBI drops suppressed/removed assemblies from that
    file entirely, so absence == suppressed/removed/version-drifted. Returns a
    set of the ORIGINAL accession strings that are present (live).

    Matching is on the numeric core (via _acc_digits), so a GCF pick counts as
    'present' if its GCA twin is in the table (and vice versa)
    """
    if not assembly_info_path or not os.path.exists(assembly_info_path):
        return set()
    wanted = {}
    for a in accessions:
        wanted.setdefault(_acc_digits(a), a)

    dead = dead_accession_cores(list(wanted.keys()), assembly_info_path)
    return {orig for d, orig in wanted.items() if d not in dead}


def gather_taxids(accessions, gtdb_tab, assembly_info_path):
    """
    accession(rootless) -> taxid, sourced tiered: GTDB table first, then the
    NCBI Parquet (ncbi-data.parquet) for the rest
    """
    acc_roots = [_norm_key(a) for a in accessions]
    taxids = {}

    gmap = gtdb_taxid_map(gtdb_tab) if gtdb_tab is not None else {}
    for a in acc_roots:
        if a in gmap:
            taxids[a] = gmap[a]

    missing = [a for a in acc_roots if a not in taxids]
    if missing:
        amap = assembly_info_taxid_map(missing, assembly_info_path)
        taxids.update(amap)

    return taxids


def parquet_lineage_resolver(assembly_info_path):
    """
    Build a resolver: taxids -> {taxid: 'd;p;c;o;f;g;s'} by reading the NCBI ranks
    straight out of the NCBI Parquet
    """
    def _resolve(taxids):
        uniq = {str(t).strip() for t in taxids if str(t).strip()}
        if not uniq or not assembly_info_path or not os.path.exists(assembly_info_path):
            return {}

        import pyarrow.parquet as pq # type: ignore

        ranks = list(RANKS)
        # push the taxid filter into the read (only rows for wanted taxids), then
        # collapse to one row per taxid IN ARROW before touching Python -- taxid
        # isn't unique per row (many assemblies share a taxid), so a naive
        # row-by-row Python loop over the filtered result is what makes this slow.
        tbl = pq.read_table(assembly_info_path, columns=["taxid"] + ranks,
                            filters=[("taxid", "in", list(uniq))])
        if tbl.num_rows == 0:
            return {}
        uniq_tbl = tbl.group_by("taxid").aggregate([(r, "min") for r in ranks])

        taxid_vals = uniq_tbl.column("taxid").to_pylist()
        rank_vals = {r: uniq_tbl.column(f"{r}_min").to_pylist() for r in ranks}
        lineages = {}
        for i, tid in enumerate(taxid_vals):
            lineages[tid] = ";".join(
                (rank_vals[r][i] if rank_vals[r][i] else "NA") for r in ranks)
        return lineages
    return _resolve


def fill_ncbi_taxonomy(merged, gtdb_tab, assembly_info_path, _resolver=None):
    """
    fill ncbi_* columns for every row, resolving each accession's NCBI tax_id
    (from local tables) to a lineage read out of the NCBI Parquet
    """
    if _resolver is None:
        _resolver = parquet_lineage_resolver(assembly_info_path)

    out = merged.copy()

    # resolve taxids for every accession (not just rows needing lineages) so the
    # taxid column is populated wherever knowable, in both GTDB and NCBI trees.
    all_accs = [out.iloc[i]["accession"] for i in range(len(out))]
    all_taxids = gather_taxids(all_accs, gtdb_tab, assembly_info_path)
    if "taxid" not in out.columns:
        out["taxid"] = "NA"
    for i in range(len(out)):
        tid = all_taxids.get(_norm_key(out.iloc[i]["accession"]))
        out.at[out.index[i], "taxid"] = str(tid) if tid else "NA"

    # which rows still need NCBI ranks
    need_idx = [i for i in range(len(out)) if str(out.iloc[i]["ncbi_domain"]) == "NA"]
    if not need_idx:
        return out

    # reuse already-gathered taxids for the rows needing lineages
    need_accs = [out.iloc[i]["accession"] for i in need_idx]
    taxids = {_norm_key(a): all_taxids[_norm_key(a)]
              for a in need_accs if _norm_key(a) in all_taxids}
    if not taxids:
        return out

    lineages = _resolver(set(taxids.values()))

    for i in need_idx:
        acc_root = _norm_key(out.iloc[i]["accession"])
        tid = taxids.get(acc_root)
        if tid is None:
            continue
        lineage = lineages.get(str(tid))
        if not lineage:
            continue
        parts = lineage.split(";")
        for ncol, val in zip(NCBI_RANKS, parts):
            out.at[out.index[i], ncol] = val if val else "NA"
    return out


# ----------------------------------------------------------------- top level --

def resolve_all(merged, gtdb_tab, assembly_info_path, _resolver=None):
    """
    fill both gtdb_* and ncbi_* on the merged table. GTDB first then NCBI
    (taxid gather + lineage read from the NCBI Parquet)
    """
    merged = fill_gtdb_taxonomy(merged, gtdb_tab)
    merged = fill_ncbi_taxonomy(merged, gtdb_tab, assembly_info_path,
                                _resolver=_resolver)
    return merged
