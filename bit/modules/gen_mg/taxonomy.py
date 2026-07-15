"""
gen-mg taxonomy-resolution layer.

Fills BOTH taxonomy rank sets (gtdb_* and ncbi_*) on the merged genome table,
regardless of how each genome entered the community:

  - GTDB ranks: looked up by accession in the GTDB metadata table. GTDB-selected
    genomes are already filled at normalization; this catches user accessions
    that happen to be in GTDB.

  - NCBI ranks: resolved from each genome's NCBI tax_id via taxonkit. Tax_ids are
    gathered cheaply from local tables first, falling back to the `datasets` CLI
    only for accessions found in neither:
        1. GTDB metadata table   (ncbi_taxid column)        -- covers GTDB picks
        2. ncbi-assembly-info.tsv (taxid, field 5)          -- covers most others
        3. datasets summary CLI                              -- last-resort fallback

A single taxonkit pass then resolves all gathered tax_ids to lineages at once.

Pure-ish logic: no printing. File/subprocess use is confined to the taxid
gathering and the taxonkit call, both injectable/mocked in tests.
"""
import os
import re
import subprocess
import tempfile
from io import StringIO

import pandas as pd # type: ignore

from bit.modules.gen_mg.selection import RANKS, GTDB_RANKS, NCBI_RANKS
from bit.modules.ncbi.get_lineage_from_taxids import get_lineage_from_taxids


def _root_acc(acc):
    """ strip version suffix for tolerant matching (GCA_000001.2 -> GCA_000001). """
    return str(acc).split(".")[0]


def _clean_gtdb_acc(acc):
    """ strip GTDB's RS_/GB_ prefix from the col-0 accession, leaving e.g.
    'RS_GCF_000005845.2' -> 'GCF_000005845.2'. """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    return s


def _norm_key(acc):
    """ the single canonical lookup key for an accession (de-prefix + de-version). """
    return _root_acc(_clean_gtdb_acc(acc))


# ---------------------------------------------------------------- GTDB ranks --

def gtdb_taxid_map(gtdb_tab):
    """ accession(rootless, both GCA & GCF forms) -> ncbi_taxid (str). """
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
    """ accession(rootless, both GCA & GCF forms) -> dict(plain_rank -> gtdb value). """
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
    """ fill gtdb_* columns for any row whose accession is in the GTDB table and
    whose gtdb_* are still NA. Returns the merged frame (modified copy). """
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
    """ accession(rootless) -> taxid from the NCBI assembly-info Parquet table.
    Only the accession and taxid columns are read, and only rows matching the
    wanted accessions are kept.

    NOTE: this is part of the taxid-resolution cascade that Phase 3 of the Parquet
    migration removes entirely (with lineages as columns, no per-accession taxid
    lookup is needed). It's converted to Parquet here only so it keeps working in
    the meantime rather than trying to open a .tsv that no longer exists.
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


def _acc_digits(acc):
    """ numeric core of an assembly accession, ignoring GCA/GCF prefix and
    version, so GCA_000001.2 and GCF_000001.1 share the key '000001'. GCA/GCF
    twins share their numeric portion. """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    s = s.split(".")[0]              # drop version
    if "_" in s:
        s = s.split("_", 1)[1]       # drop GCA_/GCF_ prefix -> digits
    return s


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

    import pyarrow.parquet as pq # type: ignore
    live = set()
    col = pq.read_table(assembly_info_path, columns=["assembly_accession"]).column(0)
    for acc in col:
        d = _acc_digits(acc.as_py())
        if d in wanted:
            live.add(d)
    return {orig for d, orig in wanted.items() if d in live}


def datasets_taxid_map(accessions):
    """ fallback: accession(rootless) -> taxid via the NCBI `datasets` CLI.
    Queries only the given (already-narrowed) accessions. """
    if not accessions:
        return {}
    fields = ["accession", "organism-tax-id"]
    try:
        summary = subprocess.Popen(
            ["datasets", "summary", "genome", "accession", *accessions, "--as-json-lines"],
            stdout=subprocess.PIPE)
        fmt = subprocess.Popen(
            ["dataformat", "tsv", "genome", "--fields", ",".join(fields)],
            stdin=summary.stdout, stdout=subprocess.PIPE)
        summary.stdout.close()
        out, _ = fmt.communicate()
    except (OSError, subprocess.SubprocessError):
        return {}
    tab = pd.read_csv(StringIO(out.decode()), sep="\t", dtype=str)
    if tab.shape[1] < 2:
        return {}
    tab.columns = ["accession", "taxid"][:tab.shape[1]]
    res = {}
    for acc, tid in zip(tab["accession"], tab["taxid"]):
        if pd.notna(acc) and pd.notna(tid) and str(tid).strip():
            res[_norm_key(acc)] = str(tid).strip()
    return res


def gather_taxids(accessions, gtdb_tab, assembly_info_path, use_datasets_fallback=True):
    """ accession(rootless) -> taxid, sourced tiered: GTDB table, then
    assembly-info.tsv, then (optionally) the datasets CLI for leftovers. """
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

    missing = [a for a in acc_roots if a not in taxids]
    if missing and use_datasets_fallback:
        dmap = datasets_taxid_map(missing)
        taxids.update(dmap)

    return taxids


def resolve_lineages(taxids, _runner=get_lineage_from_taxids):
    """ taxid -> 7-rank ';'-joined lineage string, via one taxonkit pass.
    `taxids` is an iterable of taxid strings. `_runner` is injectable for tests. """
    uniq = sorted({str(t).strip() for t in taxids if str(t).strip()})
    if not uniq:
        return {}
    with tempfile.TemporaryDirectory() as td:
        tin = os.path.join(td, "taxids.txt")
        tout = os.path.join(td, "lineages.tsv")
        with open(tin, "w") as fh:
            fh.write("\n".join(uniq) + "\n")
        _runner(tin, tout)
        lineages = {}
        with open(tout) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if not f or not f[0].strip():
                    continue
                taxid = f[0]
                ranks = []
                for r in RANKS:
                    j = idx.get(r)
                    ranks.append(f[j] if (j is not None and j < len(f) and f[j]) else "NA")
                lineages[taxid] = ";".join(ranks)
    return lineages


def fill_ncbi_taxonomy(merged, gtdb_tab, assembly_info_path,
                       use_datasets_fallback=True, _resolver=resolve_lineages):
    """ fill ncbi_* columns for every row, resolving each accession's NCBI tax_id
    (tiered sources) to a lineage. Only fills rows whose ncbi_* are still NA.

    Also persists a `taxid` column for every row whose NCBI tax_id can be resolved
    (the same tax_id is available for GTDB picks via the GTDB metadata table, so it
    is populated in both taxonomy trees), defaulting to "NA" where none resolves.

    Returns the merged frame (modified copy). """
    out = merged.copy()

    # resolve taxids for every accession (not just rows needing lineages) so the
    # taxid column is populated wherever knowable, in both GTDB and NCBI trees.
    all_accs = [out.iloc[i]["accession"] for i in range(len(out))]
    all_taxids = gather_taxids(all_accs, gtdb_tab, assembly_info_path,
                               use_datasets_fallback=use_datasets_fallback)
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

def resolve_all(merged, gtdb_tab, assembly_info_path, use_datasets_fallback=True,
                _resolver=resolve_lineages):
    """ fill both gtdb_* and ncbi_* on the merged table. GTDB first (cheap,
    in-memory), then NCBI (tiered taxid gather + one taxonkit pass). """
    merged = fill_gtdb_taxonomy(merged, gtdb_tab)
    merged = fill_ncbi_taxonomy(merged, gtdb_tab, assembly_info_path,
                                use_datasets_fallback=use_datasets_fallback,
                                _resolver=_resolver)
    return merged
