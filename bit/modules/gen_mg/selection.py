"""
gen-mg selection layer: produce one normalized genome table from
GTDB-generative selection and/or user-supplied accessions (GTDB- or NCBI-sourced)
"""
import pandas as pd # type: ignore

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
GTDB_RANKS = [f"gtdb_{r}" for r in RANKS]
NCBI_RANKS = [f"ncbi_{r}" for r in RANKS]

SCHEMA = ["accession"] + GTDB_RANKS + NCBI_RANKS + [
    "user_supplied",
    "pinned_rel_abundance", "pinned_coverage", "pinned_mutation_rate",
]

_PIN_COLS = ["pinned_rel_abundance", "pinned_coverage", "pinned_mutation_rate"]


# ---------- GTDB selection (wraps the tested one-per-rank selector) ----------

def _resolve_checkm_cols(tab):
    """ return (completeness_col, contamination_col), preferring checkm2_* (GTDB
    R214+) and falling back to checkm_* (<=R207). None if neither present. """
    if "checkm2_completeness" in tab.columns and "checkm2_contamination" in tab.columns:
        return "checkm2_completeness", "checkm2_contamination"
    if "checkm_completeness" in tab.columns and "checkm_contamination" in tab.columns:
        return "checkm_completeness", "checkm_contamination"
    return None, None


def _acc_core(acc):
    """
    Canonical accession key: strip the GTDB RS_/GB_ prefix, the version suffix, AND
    the GCA_/GCF_ prefix, so GCA_000005845.2 and GCF_000005845.1 share the key
    '000005845'
    """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    s = s.split(".")[0]              # drop version
    if "_" in s:
        s = s.split("_", 1)[1]       # drop GCA_/GCF_ prefix -> digits
    return s


# quality floor applied only when derep is off (the representative pool is already
# quality-curated; the full pool is not). Not user-tunable by design.
OFF_MODE_MIN_COMPLETENESS = 70.0
OFF_MODE_MAX_CONTAMINATION = 10.0


def _select_gtdb_one_per_rank(gtdb_tab, derep_rank="species", domains=None,
                              num_genomes=None, seed=None, exclude_groups=None,
                              exclude_accessions=None):
    """
    one GTDB genome per unique value at derep_rank; strict one-per-group.

    When derep_rank == 'off', dereplication is disabled: the pool is the
    FULL GTDB table (all genomes, so multiple strains per species are possible),
    subject only to a fixed quality floor. Otherwise the pool is GTDB species
    representatives, which is the right foundation for one-per-group selection.

    exclude_accessions: accessions (any GCA/GCF/prefixed/versioned form) never to
    select — used by the suppression backfill to skip known-dead genomes while
    still allowing other members of the same derep group
    """
    warnings = []
    if derep_rank != "off" and derep_rank not in RANKS:
        raise ValueError(f"derep_rank must be 'off' or one of {RANKS}, got '{derep_rank}'")

    if derep_rank == "off":
        # full pool (not just representatives) + quality floor
        tab = gtdb_tab
        comp_col, cont_col = _resolve_checkm_cols(tab)
        if comp_col:
            comp = pd.to_numeric(tab[comp_col], errors="coerce")
            cont = pd.to_numeric(tab[cont_col], errors="coerce")
            tab = tab[(comp >= OFF_MODE_MIN_COMPLETENESS) &
                      (cont < OFF_MODE_MAX_CONTAMINATION)]
    else:
        tab = gtdb_tab[gtdb_tab["gtdb_representative"] == "t"]

    if domains:
        wanted = {d.lower(): d for d in domains}
        tab = tab[tab["domain"].str.lower().isin(wanted.keys())]
        present = {d.lower() for d in tab["domain"].unique()}
        for d in sorted(set(wanted) - present):
            warnings.append(f"No GTDB representative genomes for requested domain '{wanted[d]}'.")

    if len(tab) == 0:
        return tab.iloc[0:0].copy(), warnings

    tab = tab.copy()

    # exclude specific known-dead accessions (suppression backfill); version- and
    # prefix-tolerant match on the GCA accession column. Matching is on the NUMERIC
    # CORE (_acc_core), so a GCF_ accession excludes its GCA_ twin and
    # vice versa
    if exclude_accessions:
        dead = {_acc_core(a) for a in exclude_accessions}
        keep = ~tab["ncbi_genbank_assembly_accession"].map(_acc_core).isin(dead)
        tab = tab[keep]
        if len(tab) == 0:
            return tab.iloc[0:0].copy(), warnings

    tab["_is_refseq_ref"] = (tab["ncbi_refseq_category"] == "reference genome").astype(int)
    tab = tab.sample(frac=1.0, random_state=seed)
    tab = tab.sort_values("_is_refseq_ref", ascending=False, kind="stable")

    if derep_rank == "off":
        selected = tab
    else:
        # exclude groups already filled (by user accessions or already-kept picks)
        if exclude_groups:
            tab = tab[~tab[derep_rank].isin(exclude_groups)]
        selected = tab.groupby(derep_rank, sort=False, as_index=False).head(1)

    num_avail = len(selected) if derep_rank == "off" else selected[derep_rank].nunique()

    if num_genomes is not None and num_avail > num_genomes:
        selected = selected.sample(n=num_genomes, random_state=seed)
    elif num_genomes is not None and num_avail < num_genomes and derep_rank != "off":
        warnings.append(
            f"Requested {num_genomes} genomes, but only {num_avail} unique "
            f"'{derep_rank}' group(s) available{' after excluding user-supplied taxa' if exclude_groups else ''}. "
            f"Returning {num_avail}. Use a finer --derep-rank or broaden --domains for more."
        )
    elif num_genomes is not None and num_avail < num_genomes and derep_rank == "off":
        warnings.append(
            f"Requested {num_genomes} genomes, but only {num_avail} genome(s) "
            f"available in this slice. Returning {num_avail}."
        )

    return selected.drop(columns=["_is_refseq_ref"]).reset_index(drop=True), warnings


def _normalize_gtdb_rows(gtdb_selected):
    """
    map a GTDB-selected subframe onto the normalized schema.

    Fills gtdb_* ranks from the GTDB taxonomy columns; ncbi_* ranks are left NA
    here and resolved later (via the accession's NCBI tax_id) by the taxonomy
    layer.
    """
    out = pd.DataFrame()
    out["accession"] = gtdb_selected["ncbi_genbank_assembly_accession"].values
    for gtdb_col, plain in zip(GTDB_RANKS, RANKS):
        out[gtdb_col] = gtdb_selected[plain].values
    for ncbi_col in NCBI_RANKS:
        out[ncbi_col] = "NA"
    out["user_supplied"] = False
    for c in _PIN_COLS:
        out[c] = pd.NA
    return out[SCHEMA].reset_index(drop=True)


# ---------- NCBI normalization (for user accessions / euk path) ----------

def _split_lineage(lineage_str):
    """
    Given a 7-rank ';'-joined lineage (domain..species), return a dict of
    plain-rank -> value. Missing/absent -> 'NA'.
    """
    vals = {r: "NA" for r in RANKS}
    if isinstance(lineage_str, str) and lineage_str.strip():
        parts = [p.strip() for p in lineage_str.split(";")]
        for r, p in zip(RANKS, parts):
            vals[r] = p if p else "NA"
    return vals


# backwards-compatible alias (older callers/tests referenced this name)
_split_ncbi_lineage = _split_lineage


def _normalize_ncbi_rows(ncbi_tab, lineage_lookup=None, user_supplied=True):
    """
    Map an NCBI metadata frame onto the normalized schema.

    Fills ncbi_* ranks from lineage_lookup (accession -> 7-rank ';'-joined
    lineage string, resolved upstream from a tax_id). gtdb_* ranks are left NA
    here and filled later by the taxonomy layer if the accession is present in
    the GTDB metadata table. Ranks are NA when no lineage is available.
    """
    rows = []
    for _, r in ncbi_tab.iterrows():
        acc = r["accession"]
        lineage = (lineage_lookup or {}).get(acc)
        ncbi_ranks = _split_lineage(lineage)
        row = {
            "accession": acc,
            "user_supplied": user_supplied,
        }
        for gr in GTDB_RANKS:
            row[gr] = "NA"
        for nr, plain in zip(NCBI_RANKS, RANKS):
            row[nr] = ncbi_ranks[plain]
        for c in _PIN_COLS:
            row[c] = pd.NA
        rows.append(row)
    return pd.DataFrame(rows, columns=SCHEMA)


# ---------- merge ----------

def merge_sources(generative_df=None, user_df=None):
    """
    Concatenate generative + user-supplied normalized frames, dedup on accession
    (user rows win, since they may carry pins), return (merged, warnings).
    """
    warnings = []
    frames = [f for f in (generative_df, user_df) if f is not None and len(f)]
    if not frames:
        return pd.DataFrame(columns=SCHEMA), ["No genomes selected from any source."]

    merged = pd.concat(frames, ignore_index=True)

    # dedup: prefer user_supplied row when an accession appears twice
    merged = merged.sort_values("user_supplied", ascending=False, kind="stable")
    before = len(merged)
    merged = merged.drop_duplicates(subset="accession", keep="first")
    if len(merged) < before:
        warnings.append(f"Removed {before - len(merged)} duplicate accession(s) across sources "
                        f"(kept user-supplied where overlapping).")

    return merged.reset_index(drop=True), warnings


def user_filled_groups(user_df, derep_rank):
    """
    set of GTDB derep_rank values occupied by user genomes (for generative
    exclusion). Uses the gtdb_<rank> column, since generative selection/derep is
    GTDB-based. Only meaningful after GTDB taxonomy has been resolved for user
    rows; returns empty set if that column is absent or all-NA
    """
    if user_df is None or len(user_df) == 0 or derep_rank == "off":
        return set()
    col = f"gtdb_{derep_rank}"
    if col not in user_df.columns:
        return set()
    return set(user_df[col].dropna().tolist()) - {"NA"}
