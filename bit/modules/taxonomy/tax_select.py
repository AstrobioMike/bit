"""
Taxon-based genome selection from the GTDB / NCBI Parquet assets

This is the code path behind:
  - `bit get-accs-from-gtdb` / `bit get-accs-from-ncbi`
  - `bit gen-mg` genome selection (gen_mg/selection.py)
"""

import pyarrow.compute as pc # type: ignore
import pyarrow.dataset as ds # type: ignore
import pyarrow.parquet as pq # type: ignore

from bit.modules.taxonomy.tax_ranks import RANKS, NA, rank_index, validate_derep_rank


class SourceSpec:
    """Per-source column names. The only thing that differs between GTDB and NCBI."""

    def __init__(self, name, acc_col, rep_filter, quality_cols,
                 ref_col, level_col=None, contig_col=None, size_col=None,
                 size_fallback_col=None):
        self.name = name
        self.acc_col = acc_col
        self.rep_filter = rep_filter
        self.quality_cols = quality_cols
        self.ref_col = ref_col
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
        size_col="genome_size",
        contig_col="contig_count",
    ),
    "ncbi": SourceSpec(
        name="ncbi",
        acc_col="assembly_accession",
        rep_filter=("refseq_category", "reference genome"),
        quality_cols=("checkm_completeness", "checkm_contamination"),
        ref_col="refseq_category",
        level_col="assembly_level",
        contig_col="contig_count",
        size_col="genome_size_ungapped",
        size_fallback_col="genome_size",
    ),
}

REFERENCE_VALUE = "reference genome"

# ranks genomes with NO checkm
ASSEMBLY_LEVEL_ORDER = {
    "complete genome": 4,
    "chromosome": 3,
    "scaffold": 2,
    "contig": 1,
}

# A RefSeq "reference genome" is preferred over a merely higher-checkm genome --
# but only if it is not actually bad. These gates are deliberately loose: they exist to
# catch junk, not to bump a reference genome that is slightly lower than another
REF_MIN_COMPLETENESS = 85.0
REF_MAX_CONTAMINATION = 10.0

MISSING_CONTAMINATION = float("inf")


# Default dereplication rank, as a function of the TARGET taxon's rank: two ranks
# finer than the target, clamped at species.
#
# I did this instead of a fixed 'species' default because derep at species level is fine for a genus wanted taxon
# but not great for a domain. E.g.:
#
#   target                 full pool   fixed 'species'   2-lower
#   Bacteria (domain)        878,998         189,801        574   <- class
#   Pseudomonadota (phylum)  324,613          46,828        413   <- order
#   Alphaproteobacteria (cl)  51,346          21,473        428   <- family
#   Escherichia (genus)       37,464              13         13   <- species
#
DEFAULT_DEREP_BY_TARGET_RANK = {
    "domain":  "class",
    "phylum":  "order",
    "class":   "family",
    "order":   "genus",
    "family":  "species",
    "genus":   "species",
    "species": None,
}

# counts outside this band get a nudge toward modifying --derep-rank
SANE_LOW, SANE_HIGH = 20, 5000


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
    given -- e.g. a name used as both an order and a family. This is the homonym
    case; on the NCBI side a caller can sidestep it entirely by passing a taxid
    (see select_by_taxid; though i might drop this, revist, Mike...)
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
    NCBI only: select by a lineage TAXID rather than a name. Unambiguous by
    construction -- this is why the NCBI
    asset carries domain_taxid..species_taxid (for now)
    """
    col = f"{rank}_taxid"
    tab = pq.read_table(path, columns=["assembly_accession"],
                        filters=[(col, "=", str(taxid))])
    return tab.column("assembly_accession").to_pylist()


# ---------------------------------------------------------------------------
# liveness screening (suppressed / removed assemblies)
# ---------------------------------------------------------------------------

def _acc_core(acc):
    """
    Canonical key for an assembly, agnostic to GTDB's RS_/GB_ prefix, the
    GCA_/GCF_ pairing, and the version suffix. 'RS_GCF_000005845.2' -> '000005845'
    """
    s = str(acc)
    if s.startswith(("RS_", "GB_")):
        s = s[3:]
    return s.split(".")[0].split("_")[-1]


def live_accession_cores(ncbi_table_path):
    """
    The set of assembly cores currently present in the NCBI assembly summary.

    NCBI drops suppressed/removed assemblies from that file entirely, so
    ABSENCE == suppressed / removed. GTDB is a snapshot so still has some of them
    """
    tab = pq.read_table(ncbi_table_path, columns=["assembly_accession"])
    return {_acc_core(a) for a in tab.column("assembly_accession").to_pylist() if a}


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
                "Target taxon is at species rank, so dereplication was turned off "
                "(one genome per species within a species returns a single genome).")
        return eff, warnings

    problem = validate_derep_rank(wanted_rank, derep_rank)
    if problem:
        raise ValueError(problem)

    return derep_rank, warnings


def size_advice(n_selected, wanted_rank, derep_rank):
    """
    Nudge at the tails. 2-lower keeps most cases in a sane band, but small phyla
    can under-sample (Cyanobacteriota -> 37 orders) and an explicit fine derep on a
    broad taxon can explode (Bacteria + species -> ~190k).
    """
    if derep_rank is None:
        return []
    d = rank_index(derep_rank)
    if n_selected > SANE_HIGH and d > 0:
        return [f"{n_selected:,} reference genomes selected -- that is a very large "
                f"tree. Consider a coarser --derep-rank (e.g. '{RANKS[d - 1]}')."]
    if n_selected < SANE_LOW and d < len(RANKS) - 1:
        return [f"Only {n_selected:,} reference genome(s) selected. Consider a finer "
                f"--derep-rank (e.g. '{RANKS[d + 1]}') for more resolution."]
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
    genome_size_ungapped / contig_count -- an N50 proxy that is normalised for genome
    size, unlike a raw contig count.

    Returns 0.0 (i.e. "worst") when it can't be computed, so a row with missing size
    or count never outranks one we can actually assess.
    """
    if not spec.contig_col or not spec.size_col:
        return 0.0

    contigs = _num(row.get(spec.contig_col))
    if not contigs or contigs <= 0:          # missing, zero, or junk
        return 0.0

    size = _num(row.get(spec.size_col))
    if size is None and spec.size_fallback_col:
        size = _num(row.get(spec.size_fallback_col))   # genome_size includes N-gaps
    if size is None or size <= 0:
        return 0.0

    return size / contigs


def quality_key(row, spec):
    """
    Deterministic "best genome in this group" ordering. Lower tuple sorts first.

    Ordering rationale:

    1. RefSeq REFERENCE genomes first (gated on not being junk). These are curated
       and recognisable -- having E. coli K-12 represent its group makes a tree far
       easier to read than an arbitrary accession. At the top end, a 99.5% vs 98.8%
       checkm difference is noise for SCG recovery, so interpretability wins. The
       gate exists so that a pathological reference cannot beat a good non-reference.

    2. Genomes WITH checkm before genomes without. Prokaryotes have checkm; euks and
       viruses do not (NCBI's report is prokaryotes-only). In practice a group is
       taxonomically homogeneous so this rarely bites, but a mixed group must not
       rank a euk above a prokaryote just because its NA sorted low.

    3. checkm: highest completeness, then lowest contamination.

    4. NO checkm (euks/viruses): best assembly_level (Complete > Chromosome >
       Scaffold > Contig), then LONGEST MEAN CONTIG. assembly_level alone is far too
       coarse for euks -- most sit at "Scaffold".

       Mean contig length, NOT raw contig count: a 12 Mb fungal genome in 200 contigs
       is better assembled than a 2 Mb genome in 150, so a bare count is confounded
       by genome size. NCBI's assembly_summary carries no N50 (the metric you'd
       actually want), so mean contig length is the best available proxy.

       It is also the mechanistically right quantity HERE. GToTree recovers SCGs by
       calling genes and HMM-searching them; a gene straddling a contig break gets
       truncated and missed, and the genome then loses hits and may be dropped by the
       genome-hits cutoff. The chance any given gene is broken scales with break
       DENSITY (contigs per base) -- whose reciprocal is exactly mean contig length.

       Numerator is genome_size_ungapped: contigs contain no N-gaps by definition,
       whereas genome_size includes scaffolding gaps. Falls back to genome_size.

    5. accession, purely so the result is stable/reproducible.
    """
    comp_col, cont_col = spec.quality_cols
    comp = _num(row.get(comp_col))
    cont = _num(row.get(cont_col))
    has_checkm = comp is not None

    is_ref = str(row.get(spec.ref_col) or "").strip().lower() == REFERENCE_VALUE
    if is_ref and has_checkm:
        # only trust the reference flag if the genome isn't junk
        is_ref = (comp >= REF_MIN_COMPLETENESS
                  and (cont is None or cont <= REF_MAX_CONTAMINATION))

    lvl = ASSEMBLY_LEVEL_ORDER.get(
        str(row.get(spec.level_col) or "").strip().lower(), 0) if spec.level_col else 0

    return (
        0 if is_ref else 1,           # 1. reference genomes first
        0 if has_checkm else 1,       # 2. checkm'd (prokaryotes) before not
        -(comp if has_checkm else 0.0),       # 3. completeness (negated: higher better)
        (cont if cont is not None                #    contamination (lower better;
         else (MISSING_CONTAMINATION if has_checkm else 0.0)),
                                              #    unknown sorts WORST, not best)
        -lvl,                         # 4. euk fallback: assembly level
        -mean_contig_length(row, spec),  #    ...then longest mean contig
        row.get(spec.acc_col) or "",  # 5. stable tiebreak
    )


def derep_rank(path, source, wanted_rank, wanted_taxon, per_rank,
                 reps_only=True, pick="quality", screen_against=None):
    """
    One genome per unique value of `per_rank`, within `wanted_taxon`.

    E.g. --wanted-ref-taxon bacteria --one-genome-per-rank class
         -> one genome for each bacterial class.

    reps_only:
        GTDB -- keep True. GTDB representatives are comprehensive (~one per species
        cluster), so one-per-class over the rep pool gives excellent coverage.

        NCBI -- consider False. RefSeq "reference genomes" are SPARSE (a few
        thousand genome-wide), so many classes have zero and the result silently
        under-covers. With reps_only=False the pool is all genomes, ranked by
        assembly_level instead.

    pick:
        "quality" -- deterministic best-per-group (GToTree: a fragmentary genome
                     just loses SCG hits and gets filtered out anyway).
        "first"   -- first in lineage-sorted order; cheap and deterministic.

        gen_mg wants a *seeded random* pick preferring RefSeq reference genomes
        instead; that's a different policy, so it passes its own picker rather
        than reusing these.

    screen_against:
        Path to the NCBI Parquet asset. When given, the candidate pool is
        pre-filtered to assemblies that still exist at NCBI BEFORE grouping.

        This matters MUCH more here than in the un-dereplicated case: losing 1 of
        8,830 Cyanobacteriota genomes is noise, but if the single genome chosen to
        represent an entire class has been suppressed, that whole class silently
        vanishes from the tree.

        Because the pick is DETERMINISTIC, pre-filtering the pool is equivalent to
        backfilling a dead pick -- no iterative re-draw loop is needed (gen_mg
        needs one only because it picks randomly and must return exactly N).

    Returns (accessions, groups_seen, warnings).
    """
    warnings = []

    problem = validate_one_per_rank(wanted_rank, per_rank)
    if problem:
        raise ValueError(problem)

    if rank_index(per_rank) == rank_index(wanted_rank):
        warnings.append(
            f"--one-genome-per-rank '{per_rank}' is the same rank as the target "
            f"taxon, so this will return a single genome.")

    spec = SOURCES[source]

    cols = [spec.acc_col, per_rank, spec.ref_col] + list(spec.quality_cols)
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
                if r.get(spec.acc_col) and _acc_core(r[spec.acc_col]) in live]
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

    best = {}
    n_na_group = 0
    for row in rows:
        group = row.get(per_rank)
        # a genome whose per_rank is unnamed must NOT become its own "NA" group --
        # that would silently allocate a genome to a group that isn't real.
        if not group or group == NA:
            n_na_group += 1
            continue
        if pick == "first":
            best.setdefault(group, row)
            continue
        k = quality_key(row, spec)
        if group not in best or k < quality_key(best[group], spec):
            best[group] = row

    if n_na_group:
        warnings.append(
            f"{n_na_group:,} genome(s) have no assigned '{per_rank}' and were "
            f"skipped (unnamed/unclassified at that rank).")

    accs = [best[g][spec.acc_col] for g in sorted(best)]
    return accs, len(best), warnings
