"""
gen-metagenome selection layer: produce one normalized genome table from
GTDB-generative selection and/or user-supplied accessions (GTDB- or NCBI-sourced).

Normalized schema (one row per selected genome):
  accession, source_db, domain, phylum, class, order, family, genus, species,
  genome_size, taxonomy_source, user_supplied
plus optional pinned override columns carried through from a user TSV:
  pinned_rel_abundance, pinned_coverage, pinned_mutation_rate

Pure logic: no printing, no sys.exit, no file writes. Returns (df, warnings).
"""
import pandas as pd

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

SCHEMA = ["accession", "source_db"] + RANKS + [
    "metadata_genome_size", "taxonomy_source", "user_supplied",
    "pinned_rel_abundance", "pinned_coverage", "pinned_mutation_rate",
]


# ---------- GTDB selection (wraps the tested one-per-rank selector) ----------

def _select_gtdb_one_per_rank(gtdb_tab, derep_rank="species", domains=None,
                              num_genomes=None, seed=None, exclude_groups=None):
    """ one GTDB genome per unique value at derep_rank; strict one-per-group. """
    warnings = []
    if derep_rank != "off" and derep_rank not in RANKS:
        raise ValueError(f"derep_rank must be 'off' or one of {RANKS}, got '{derep_rank}'")

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
    tab["_is_refseq_ref"] = (tab["ncbi_refseq_category"] == "reference genome").astype(int)
    tab = tab.sample(frac=1.0, random_state=seed)
    tab = tab.sort_values("_is_refseq_ref", ascending=False, kind="stable")

    if derep_rank == "off":
        selected = tab
    else:
        # exclude groups already filled by user-supplied accessions
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
            f"Requested {num_genomes} genomes, but only {num_avail} representative genome(s) "
            f"available in this slice. Returning {num_avail}."
        )

    return selected.drop(columns=["_is_refseq_ref"]).reset_index(drop=True), warnings


def _normalize_gtdb_rows(gtdb_selected):
    """ map a GTDB-selected subframe onto the normalized schema. """
    out = pd.DataFrame()
    out["accession"] = gtdb_selected["ncbi_genbank_assembly_accession"]
    out["source_db"] = "GTDB"
    for r in RANKS:
        out[r] = gtdb_selected[r].values
    # GTDB metadata carries genome_size; fall back to NA if absent
    out["metadata_genome_size"] = gtdb_selected["genome_size"].values if "genome_size" in gtdb_selected.columns else pd.NA
    out["taxonomy_source"] = "GTDB"
    out["user_supplied"] = False
    for c in ("pinned_rel_abundance", "pinned_coverage", "pinned_mutation_rate"):
        out[c] = pd.NA
    return out.reset_index(drop=True)


# ---------- NCBI normalization (for user accessions / euk path) ----------

def _split_ncbi_lineage(lineage_str):
    """
    Given a 7-rank ';'-joined lineage (domain..species) resolved upstream from a
    tax_id, return a dict of rank->value. Missing ranks -> 'NA'. If lineage is
    absent, all ranks 'NA'.
    """
    vals = {r: "NA" for r in RANKS}
    if isinstance(lineage_str, str) and lineage_str.strip():
        parts = [p.strip() for p in lineage_str.split(";")]
        for r, p in zip(RANKS, parts):
            vals[r] = p if p else "NA"
    return vals


def _normalize_ncbi_rows(ncbi_tab, lineage_lookup=None, user_supplied=True):
    """
    Map an NCBI metadata frame (run_ncbi_summary schema) onto the normalized schema.
    lineage_lookup: optional dict accession-> 7-rank ';'-joined lineage string.
    Ranks are NA when no lineage is available (v1 behavior for euk accessions
    until tax_id->lineage resolution is wired in).
    """
    rows = []
    for _, r in ncbi_tab.iterrows():
        acc = r["accession"]
        lineage = (lineage_lookup or {}).get(acc)
        ranks = _split_ncbi_lineage(lineage)
        size = r.get("total_sequence_length", "NA")
        try:
            size = int(size)
        except (TypeError, ValueError):
            size = pd.NA
        rows.append({
            "accession": acc,
            "source_db": r.get("source_database", "NCBI"),
            **ranks,
            "metadata_genome_size": size,
            "taxonomy_source": "GTDB" if lineage and "_NCBI" not in str(lineage) else "NCBI",
            "user_supplied": user_supplied,
            "pinned_rel_abundance": pd.NA,
            "pinned_coverage": pd.NA,
            "pinned_mutation_rate": pd.NA,
        })
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
    """ set of derep_rank values occupied by user genomes (for generative exclusion). """
    if user_df is None or len(user_df) == 0 or derep_rank == "off":
        return set()
    vals = set(user_df[derep_rank].dropna().tolist()) - {"NA"}
    return vals
