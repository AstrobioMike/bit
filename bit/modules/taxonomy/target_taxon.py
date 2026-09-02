"""
`-t/--target-taxon` resolution for the surfaces that DOWNLOAD or CONSUME genomes rather
than just listing accessions.

This is bit's counterpart to GToTree's utils/taxonomy/wanted_ref_tax.py: one place that
turns a taxon name (+ the selection knobs) into a list of assembly accessions
"""

from bit.modules.taxonomy.tax_derep import select_ref_genomes
from bit.modules.taxonomy.tax_targets import expand_all_targets, is_all_target
from bit.modules.gtdb.get_gtdb_data import gtdb_data_table_path, get_gtdb_data
from bit.modules.ncbi.get_ncbi_assembly_data import (ncbi_data_table_path,
                                                     get_ncbi_assembly_data)


# --source value -> (taxonomy-core source key, path fn, fetch fn)
_SOURCE_ASSETS = {
    "gtdb": ("gtdb", gtdb_data_table_path, get_gtdb_data),
    "ncbi": ("ncbi", ncbi_data_table_path, get_ncbi_assembly_data),
}

SOURCE_CHOICES = tuple(_SOURCE_ASSETS)

# --ncbi-section value -> the accession prefixes it admits
NCBI_SECTION_PREFIXES = {
    "refseq": ("GCF_",),
    "genbank": ("GCA_",),
    "both": None,
}

SECTION_CHOICES = tuple(NCBI_SECTION_PREFIXES)


class TargetTaxonError(Exception):
    """A `-t` request that resolved to nothing usable."""


def _assets_for_source(source):
    key = str(source).strip().lower()
    if key not in _SOURCE_ASSETS:
        raise TargetTaxonError(
            f"'{source}' is not a recognized --source "
            f"(expected one of: {', '.join(SOURCE_CHOICES)}).")
    return _SOURCE_ASSETS[key]


def ensure_source_data(source, quiet=True):
    """
    Make sure the asset for `source` is local, and return (core_source, table_path)
    """
    core_source, table_path_fn, fetch_fn = _assets_for_source(source)
    fetch_fn(quiet=quiet)
    if core_source == "gtdb":
        get_ncbi_assembly_data(quiet=quiet)
    return core_source, table_path_fn()


def section_prefixes(source, ncbi_section):
    """
    The accession prefixes for an --ncbi-section value, or None for no restriction.

    Meaningless for GTDB, whose accession column is already GenBank-only, so this
    returns None there rather than filtering out everything.
    """
    if str(source).strip().lower() == "gtdb":
        return None
    return NCBI_SECTION_PREFIXES.get(str(ncbi_section).strip().lower())


def expand_target_taxa(source, taxa):
    """
    Expand any 'all' in `taxa` into the per-domain targets the asset actually holds.

    Returns (expanded_taxa, domains_expanded_from_all). `domains` is empty when no
    'all' was given, which is also the signal that no expansion note is needed.
    """
    taxa = list(taxa or [])
    if not any(is_all_target(t) for t in taxa):
        return taxa, []

    core_source, table_path = ensure_source_data(source)
    return expand_all_targets(table_path, core_source, taxa)


def describe_all_expansion(source, domains):
    """The one-line 'this is what all meant' note, or None."""
    if not domains:
        return None
    return (f"`-t all` was expanded to the domains in the {str(source).upper()} "
            f"table: {', '.join(domains)}.")


def resolve_target_taxon_accessions(source, taxon, target_rank=None,
                                    derep_rank="auto", target_domain=None,
                                    ncbi_section="refseq", assembly_levels=None,
                                    reps_only=None, min_completeness=None,
                                    max_contamination=None, include_rows=True,
                                    exclude_cores=None):
    """
    Resolve one `-t <taxon>` to (accessions, selection).

    `selection` is the RefGenomeSelection it came from, carrying the canonical name,
    resolved rank, effective derep rank, metadata rows, and any warnings -- everything
    the reporting layer needs to explain what a taxon turned into.

    A GTDB-sourced selection is screened against the NCBI table for liveness, so
    suppressed/removed assemblies are dropped BEFORE they reach the downloader (where
    they would otherwise surface as a confusing "not found at NCBI" report).

    `exclude_cores` is an optional set of accession cores from `--exclusion-list`. It
    goes to the selection core rather than being applied to the result, so the
    exclusion happens against the candidate pool before dereplication.

    Raises TargetTaxonError for an unknown source or an empty selection; lets
    TaxonNotFound / AmbiguousTaxon / CrossDomainTaxon / ValueError propagate for the
    CLI layer to translate into its own wording.
    """
    core_source, table_path = ensure_source_data(source)

    screen_against = ncbi_data_table_path() if core_source == "gtdb" else None

    selection = select_ref_genomes(
        table_path, core_source, taxon,
        target_rank=target_rank,
        derep_rank=derep_rank,
        target_domain=target_domain,
        reps_only=reps_only,
        accession_prefixes=section_prefixes(source, ncbi_section),
        assembly_levels=assembly_levels,
        screen_against=screen_against,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        exclude_cores=exclude_cores,
        include_rows=include_rows)

    if not selection.accessions:
        detail = ""
        if min_completeness is not None or max_contamination is not None:
            detail = (" No genomes cleared the requested quality floor, so you may "
                      "want to relax `--min-completeness` / `--max-contamination`.")
        elif exclude_cores:
            detail = (" Every candidate genome for it was named in the "
                      "`--exclusion-list`.")
        raise TargetTaxonError(
            f"No accessions were found for the -t target "
            f"'{selection.canonical}'.{detail}")

    return selection.accessions, selection
