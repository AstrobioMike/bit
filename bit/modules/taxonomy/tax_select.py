"""
Taxon-based genome/accession selection from the GTDB / NCBI Parquet assets

This is the code:
  - `bit get-accs-from-gtdb`  / `bit get-accs-from-ncbi`
  - `bit gen-mg` genome selection (gen_mg/selection.py)
"""

import hashlib

import pyarrow.compute as pc # type: ignore
import pyarrow.parquet as pq # type: ignore

from bit.modules.taxonomy.tax_ranks import (RANKS, NA, REFERENCE_VALUE,
                                            accession_core, rank_index,
                                            validate_derep_rank)


class SourceSpec:
    """Per-source column names. The only thing that differs between GTDB and NCBI."""

    def __init__(self, name, acc_col, rep_filter, quality_cols,
                 ref_col, default_reps_only, level_col=None, contig_col=None,
                 size_col=None, size_fallback_col=None):
        self.name = name
        self.acc_col = acc_col
        self.rep_filter = rep_filter
        self.quality_cols = quality_cols
        self.ref_col = ref_col
        self.default_reps_only = default_reps_only
        self.level_col = level_col
        self.contig_col = contig_col
        self.size_col = size_col
        self.size_fallback_col = size_fallback_col


SOURCES = {
    "gtdb": SourceSpec(
        name="gtdb",
        acc_col="ncbi_genbank_assembly_accession",
        rep_filter=("gtdb_representative", "t"),
        quality_cols=("checkm2_completeness", "checkm2_contamination"),
        ref_col="ncbi_refseq_category",
        default_reps_only=True,
        size_col="genome_size",
        contig_col="contig_count",
    ),
    "ncbi": SourceSpec(
        name="ncbi",
        acc_col="assembly_accession",
        rep_filter=("refseq_category", "reference genome"),
        quality_cols=("checkm_completeness", "checkm_contamination"),
        ref_col="refseq_category",
        default_reps_only=False,
        level_col="assembly_level",
        contig_col="contig_count",
        size_col="genome_size_ungapped",
        size_fallback_col="genome_size",
    ),
}

REFERENCE_VALUE = "reference genome"

ASSEMBLY_LEVEL_ORDER = {
    "complete genome": 4,
    "chromosome": 3,
    "scaffold": 2,
    "contig": 1,
}

# A RefSeq "reference genome" is preferred over a merely higher-checkm genome,
# but only if it is not actually bad. These gates are deliberately loose, here to
# catch junk, not to replace all "reference" genomes
REF_MIN_COMPLETENESS = 85.0
REF_MAX_CONTAMINATION = 10.0

MISSING_CONTAMINATION = float("inf")


# How many ranks FINER than the target the default derep-rank
DEREP_STEPS = 2

# The default derep rank is DERIVED, not hand-written, from one rule:
#
#     "two ranks finer than the target, clamped at species, and off if that lands
#      on species"
#
# The second half of that rule is the important half. Dereplicating at the SAME rank as
# the target always yields exactly ONE genome (there is only one group: the target
# itself). That is true at EVERY rank, not just species:
#
#     --wanted-ref-tax Bacteria            --derep-rank domain  -> 1 genome
#     --wanted-ref-tax Alphaproteobacteria --derep-rank class   -> 1 genome
#     --wanted-ref-tax Escherichia         --derep-rank genus   -> 1 genome
#
# I went with this instead of a fixed derep-rank default because, e.g., 'species' is
# good for a genus but bad for a domain, e.g.:
#
#   target                 full pool   fixed 'species'   2-lower
#   Bacteria (domain)        878,998         189,801        574   <- class
#   Pseudomonadota (phylum)  324,613          46,828        413   <- order
#   Alphaproteobacteria (cl)  51,346          21,473        428   <- family
#   Escherichia (genus)       37,464              13         13   <- species


def _derive_default_derep(target_rank):
    """Two ranks finer, clamped at species; None if that would be the target's own rank."""
    i = rank_index(target_rank)
    j = min(i + DEREP_STEPS, len(RANKS) - 1)
    return None if j == i else RANKS[j]


DEFAULT_DEREP_BY_TARGET_RANK = {r: _derive_default_derep(r) for r in RANKS}
#   domain->class  phylum->order  class->family  order->genus
#   family->species  genus->species (clamped)  species->None (degenerate -> off)

# counts outside this band get a nudge toward overriding --derep-rank
SANE_LOW, SANE_HIGH = 10, 5000


class TaxonNotFound(Exception):
    pass


class AmbiguousTaxon(Exception):
    """The taxon name exists at more than one rank; caller must disambiguate."""

    def __init__(self, taxon, ranks_found):
        self.taxon = taxon
        self.ranks_found = ranks_found
        super().__init__(
            f"'{taxon}' occurs at more than one rank ({', '.join(ranks_found)}); "
            f"specify which rank is wanted")


# ---------------------------------------------------------------------------
# resolving a user-supplied taxon
# ---------------------------------------------------------------------------

def find_ranks_for_taxon(path, taxon):
    """
    Which of the 7 ranks contain `taxon`? Case-insensitive.
    Returns (canonical_name, [ranks]). Reads one column at a time, so this stays
    cheap even on the 4M-row NCBI table.
    """
    target = str(taxon).strip().lower()
    canonical = None
    found = []
    for rank in RANKS:
        col = pq.read_table(path, columns=[rank]).column(rank)
        uniq = pc.unique(col).to_pylist()
        for name in uniq:
            if name is not None and name != NA and name.lower() == target:
                canonical = name
                found.append(rank)
                break
    if canonical is None:
        raise TaxonNotFound(f"'{taxon}' doesn't exist at any rank in this source")
    return canonical, found


def resolve_taxon(path, taxon, rank=None):
    """
    Resolve a user-supplied taxon (+ optional explicit rank) to (canonical, rank).

    Raises AmbiguousTaxon if the name lives at multiple ranks and no rank was
    given. E.g., a name used as both an order and a family. On the NCBI side a user can
    sidestep this entirely by passing a taxid for now (though i might remove this; revisit Mike)
    """
    canonical, found = find_ranks_for_taxon(path, taxon)

    if rank:
        r = str(rank).strip().lower()
        rank_index(r)
        if r not in found:
            raise TaxonNotFound(
                f"'{canonical}' exists at rank(s) {', '.join(found)}, not '{r}'")
        return canonical, r

    if len(found) > 1:
        raise AmbiguousTaxon(canonical, found)

    return canonical, found[0]


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------

def select(path, source, rank, taxon, reps_only=False, columns=None):
    """
    All rows under `taxon` at `rank`. Returns a pyarrow Table.

    The rank predicate is pushed down to Parquet, so on the lineage-sorted assets
    this skips whole row groups rather than scanning.
    """
    spec = SOURCES[source]
    filters = [(rank, "=", taxon)]
    if reps_only:
        filters.append((spec.rep_filter[0], "=", spec.rep_filter[1]))

    cols = columns or [spec.acc_col]
    # always need the accession back
    if spec.acc_col not in cols:
        cols = [spec.acc_col] + list(cols)

    return pq.read_table(path, columns=cols, filters=filters)


def select_accessions(path, source, rank, taxon, reps_only=False):
    """The common case: just the accession list."""
    spec = SOURCES[source]
    tab = select(path, source, rank, taxon, reps_only=reps_only)
    return [a for a in tab.column(spec.acc_col).to_pylist() if a and a != NA]


def select_by_taxid(path, rank, taxid):
    """
    NCBI only: select by a lineage TAXID rather than a name
    """
    col = f"{rank}_taxid"
    tab = pq.read_table(path, columns=["assembly_accession"],
                        filters=[(col, "=", str(taxid))])
    return tab.column("assembly_accession").to_pylist()


# ---------------------------------------------------------------------------
# liveness screening (suppressed / removed / version-drifted assemblies)
# ---------------------------------------------------------------------------

def live_accession_cores(ncbi_table_path):
    """
    The set of assembly core accs currently present in the NCBI assembly summary

    NCBI drops suppressed/removed assemblies from the assembly summary files, so
    ABSENCE == suppressed / removed. GTDB is a snapshot and can still have them
    """
    tab = pq.read_table(ncbi_table_path, columns=["assembly_accession"])
    return {accession_core(a) for a in tab.column("assembly_accession").to_pylist() if a}


# ---------------------------------------------------------------------------
# derep rank resolution
# ---------------------------------------------------------------------------

def default_derep_rank(wanted_rank):
    """The default derep rank for a target at `wanted_rank`. None == off."""
    return DEFAULT_DEREP_BY_TARGET_RANK[RANKS[rank_index(wanted_rank)]]


def resolve_derep_rank(wanted_rank, derep_rank="auto"):
    """
    Work out the effective derep rank. Returns (rank_or_None, warnings);
    None means no dereplication.

    derep_rank:
      "auto"            -> DEFAULT_DEREP_BY_TARGET_RANK (two ranks finer, clamped)
      "off"/None/"none" -> no dereplication
      an explicit rank  -> honoured, but must not be COARSER than the target
    """
    warnings = []
    w_rank = RANKS[rank_index(wanted_rank)]

    if derep_rank in (None, "off", "none", "None"):
        return None, warnings

    if derep_rank == "auto":
        eff = DEFAULT_DEREP_BY_TARGET_RANK[w_rank]
        if eff is None:
            warnings.append(
                f"Dereplication is off by default for a target at '{w_rank}' rank: the "
                f"default derep rank would be '{w_rank}' itself, which returns a single "
                f"genome. For this run, all genomes under the taxon will be used. If you DID want a "
                f"single best representative (e.g., as an outgroup), pass "
                f"--derep-rank {w_rank} explicitly.")
        return eff, warnings

    problem = validate_derep_rank(wanted_rank, derep_rank)
    if problem:
        raise ValueError(problem)

    return derep_rank, warnings


def size_advice(n_selected, wanted_rank, derep_rank):
    """
    Nudge at the tails.

    The 2-rank-lower default derep-rank keeps most cases in a sane band
      - small phyla may have very few orders
      - an explicit fine derep on a broad taxon can explode (e.g., Bacteria + species -> ~190k)
      - and derep OFF can explode, e.g., `--wanted-ref-tax "Escherichia coli"`
        defaults to derep off and can return tens of thousands of genomes
    """
    if derep_rank is None:
        if n_selected > SANE_HIGH:
            w = rank_index(wanted_rank)
            if w == len(RANKS) - 1:
                return [f"{n_selected:,} reference genomes selected with dereplication "
                        f"off. That is a very large tree. There is no rank finer than "
                        f"'{wanted_rank}' to group by, so the options are all of them like this, or "
                        f"a single best representative (--derep-rank {wanted_rank})."]
            suggestion = DEFAULT_DEREP_BY_TARGET_RANK[RANKS[w]]
            return [f"{n_selected:,} reference genomes selected with dereplication off. "
                    f"That is a very large tree. Consider --derep-rank "
                    f"'{suggestion}'."]
        return []

    d = rank_index(derep_rank)
    if n_selected > SANE_HIGH and d > 0:
        return [f"{n_selected:,} reference genomes selected. That is a very large "
                f"tree. Consider a coarser --derep-rank (e.g., '{RANKS[d - 1]}')."]
    if n_selected < SANE_LOW and d < len(RANKS) - 1:
        return [f"Only {n_selected:,} reference genome(s) selected. Consider a finer "
                f"--derep-rank (e.g., '{RANKS[d + 1]}') for more genomes to be included."]
    return []


# ---------------------------------------------------------------------------
# one genome per rank
# ---------------------------------------------------------------------------

def _num(value, default=None):
    if value in (None, "", NA, "na"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def mean_contig_length(row, spec):
    """
    genome_size_ungapped / contig_count for a proxy of quality for those with no checkm values
    """
    if not spec.contig_col or not spec.size_col:
        return 0.0

    contigs = _num(row.get(spec.contig_col))
    if not contigs or contigs <= 0:
        return 0.0

    size = _num(row.get(spec.size_col))
    if size is None and spec.size_fallback_col:
        size = _num(row.get(spec.size_fallback_col))
    if size is None or size <= 0:
        return 0.0

    return size / contigs


def first_key(row, spec):
    """Cheapest possible picker: first row in lineage-sorted order. Deterministic."""
    return (row.get(spec.acc_col) or "",)


def seeded_random_picker(seed, prefer_references=True):
    """
    gen_mg's pick policy, as a picker you can hand to derep(pick=...).

    gen_mg does NOT want the best genome per group, we want an
    unbiased draw within the group
    """
    def pick(row, spec):
        acc = row.get(spec.acc_col) or ""
        h = hashlib.blake2b(f"{seed}:{acc}".encode(), digest_size=8).digest()
        r = int.from_bytes(h, "big")

        if prefer_references:
            is_ref = str(row.get(spec.ref_col) or "").strip().lower() == REFERENCE_VALUE
            return (0 if is_ref else 1, r, acc)
        return (r, acc)

    return pick


#: pickers addressable by name from a CLI flag if i want
PICKERS = {
    "quality": None,
    "first": first_key,
}


def resolve_picker(pick):
    """
    `pick` may be a NAME ("quality" / "first") or a CALLABLE (row, spec) -> sort key
    """
    if callable(pick):
        return pick
    if isinstance(pick, str) and pick in PICKERS:
        return PICKERS[pick]
    raise ValueError(
        f"unknown pick policy {pick!r}. Pass one of {sorted(PICKERS)}, or a callable "
        f"(row, spec) -> sort key (see seeded_random_picker).")


def quality_key(row, spec):
    """
    Deterministic "best genome in this group" ordering. Lower tuple sorts first.

    Ordering rationale:

    1. RefSeq REFERENCE genomes first (gated on not being too crappy)

    2. Genomes WITH checkm before genomes without. Prokaryotes have checkm and euks don't.
       In practice a group is shouldn't span those, so this should never matter, but a mixed group
       shouldn't rank a euk above a prokaryote just because its NA sorted low

    3. checkm: highest completeness, then lowest contamination

    4. NO checkm (euks): best assembly_level (Complete > Chromosome >
       Scaffold > Contig), then LONGEST MEAN CONTIG

    5. accession, so the result is stable/reproducible
    """
    comp_col, cont_col = spec.quality_cols
    comp = _num(row.get(comp_col))
    cont = _num(row.get(cont_col))
    has_checkm = comp is not None

    is_ref = str(row.get(spec.ref_col) or "").strip().lower() == REFERENCE_VALUE
    if is_ref and has_checkm:
        is_ref = (comp >= REF_MIN_COMPLETENESS
                  and (cont is None or cont <= REF_MAX_CONTAMINATION))

    lvl = ASSEMBLY_LEVEL_ORDER.get(
        str(row.get(spec.level_col) or "").strip().lower(), 0) if spec.level_col else 0

    return (
        0 if is_ref else 1,                                      # 1. reference genomes first
        0 if has_checkm else 1,                                  # 2. checkm'd before not
        -(comp if has_checkm else 0.0),                          # 3. completeness (higher bette)
        (cont if cont is not None                                #    contamination (lower better)
         else (MISSING_CONTAMINATION if has_checkm else 0.0)),   #    unknown sorts worst)
        -lvl,                                                    # 4. euk fallback, assembly level
        -mean_contig_length(row, spec),                          #    then longest mean contig
        row.get(spec.acc_col) or "",                             # 5. accession tiebreak
    )


PICKERS["quality"] = quality_key


def derep(path, source, wanted_rank, wanted_taxon, derep_rank,
          reps_only=None, pick="quality", screen_against=None):
    """
    One genome per unique value of `derep_rank`, within `wanted_taxon`

    E.g. --wanted-ref-tax bacteria --derep-rank class
         -> one genome for each bacterial class

    reps_only:
        None (default) -> use the SOURCE's default: True for GTDB, False for NCBI.

        This must not be a plain True. GTDB representatives are comprehensive, but
        RefSeq "reference genomes" are sparse

    pick:
        A picker NAME or a CALLABLE (row, spec) -> sort key (lowest wins).

        "quality" (default) -- best-per-group
        "first"             -- first in lineage-sorted order
        callable            -- gen_mg uses `seeded_random_picker(seed)`: references first, then
                               uniformly RANDOM within a tier, because for metagenome
                               simulation we want an unbiased draw rather than
                               systematically the most-complete things


    screen_against:
        Path to the NCBI Parquet asset. When given, the candidate pool is
        pre-filtered to assemblies that still exist at NCBI BEFORE grouping.

    Returns (accessions, groups_seen, warnings).
    """
    warnings = []

    problem = validate_derep_rank(wanted_rank, derep_rank)
    if problem:
        raise ValueError(problem)

    if rank_index(derep_rank) == rank_index(wanted_rank):
        # Sometimes may be wanted, as in "the single best genome for taxon X" is
        # how you could pick an OUTGROUP. But it is also an easy mis-type for
        # someone who wanted a tree SPANNING the taxon, and the difference is 1 genome
        # vs hundreds. So mentioning either way
        finer = RANKS[rank_index(derep_rank) + 1] \
            if rank_index(derep_rank) < len(RANKS) - 1 else None
        msg = (f"--derep-rank '{derep_rank}' is the SAME rank as the target taxon, so "
               f"exactly ONE genome will be selected (the single best genome for "
               f"'{wanted_taxon}').")
        if finer:
            msg += (f" If you wanted a tree spanning '{wanted_taxon}', use a finer "
                    f"--derep-rank (e.g. '{finer}').")
        warnings.append(msg)

    spec = SOURCES[source]
    if reps_only is None:
        reps_only = spec.default_reps_only

    cols = [spec.acc_col, derep_rank, spec.ref_col] + list(spec.quality_cols)
    for extra in (spec.level_col, spec.contig_col, spec.size_col,
                  spec.size_fallback_col):
        if extra and extra not in cols:
            cols.append(extra)

    tab = select(path, source, wanted_rank, wanted_taxon,
                 reps_only=reps_only, columns=cols)

    if screen_against:
        live = live_accession_cores(screen_against)
        before = tab.num_rows
        rows = [r for r in tab.to_pylist()
                if r.get(spec.acc_col) and accession_core(r[spec.acc_col]) in live]
        n_dead = before - len(rows)
        if n_dead:
            warnings.append(
                f"{n_dead:,} candidate genome(s) are no longer available at NCBI "
                f"(suppressed/removed) and were excluded before selection.")
        tab = None
    else:
        rows = tab.to_pylist()

    if not rows:
        warnings.append(f"No genomes found under {wanted_rank} '{wanted_taxon}'"
                        + (" in the representatives pool." if reps_only else "."))
        return [], 0, warnings

    keyfn = resolve_picker(pick)

    best = {}
    n_na_group = 0
    for row in rows:
        group = row.get(derep_rank)
        if not group or group == NA:
            n_na_group += 1
            continue
        k = keyfn(row, spec)
        if group not in best or k < keyfn(best[group], spec):
            best[group] = row

    if n_na_group:
        warnings.append(
            f"{n_na_group:,} genome(s) have no assigned '{derep_rank}' and were "
            f"skipped (unnamed/unclassified at that rank).")

    accs = [best[g][spec.acc_col] for g in sorted(best)]
    return accs, len(best), warnings
