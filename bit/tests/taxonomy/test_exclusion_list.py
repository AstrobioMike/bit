"""
`--exclusion-list` keeps named accessions out of what `-t` pulls.

The exclusion is applied INSIDE the selection core, against the candidate pool, before
dereplication and before any best-per-group pick. That ordering is the point: excluding
a species' chosen representative promotes the next-best genome in that species rather
than dropping the species. Post-filtering a dereplicated set would do the latter.

This mirrors GToTree's gtotree/tests/utils/misc/test_exclusion_list.py -- the two
implementations are meant to behave identically, so the coverage is kept in step.
"""

import types

import pyarrow as pa  # type: ignore
import pyarrow.parquet as pq  # type: ignore
import pytest  # type: ignore

from bit.modules.taxonomy.tax_ranks import RANKS
from bit.modules.taxonomy.tax_derep import select_ref_genomes
from bit.modules.taxonomy.tax_counts import count_genomes, derep_size
from bit.modules.taxonomy.exclusion_list import (load_exclusion_cores,
                                                 read_exclusion_list,
                                                 filter_rows_by_exclusion,
                                                 filter_accessions_by_exclusion,
                                                 filter_table_by_exclusion,
                                                 exclusion_warning,
                                                 exclusion_list_help)


def _write(path, entries):
    path.write_text("\n".join(entries) + "\n")
    return str(path)


# ---------------------------------------------------------------------------
# a tiny GTDB-shaped asset: one genus, two species, two genomes each
# ---------------------------------------------------------------------------

_EXTRA_COLS = ["gtdb_representative", "ncbi_refseq_category",
               "checkm2_completeness", "checkm2_contamination",
               "genome_size", "contig_count"]


def _rec(acc, species, comp):
    d = {"ncbi_genbank_assembly_accession": acc,
         "gtdb_representative": "t", "ncbi_refseq_category": "",
         "checkm2_completeness": comp, "checkm2_contamination": "0.5",
         "genome_size": "4000000", "contig_count": "1"}
    lineage = ("Bacteria", "Testophyla", "ClassA", "OrdA", "FamA", "GenA", species)
    for i, r in enumerate(RANKS):
        d[r] = lineage[i]
    return d


# GenA sp1's best genome is ...01 (99.0); its runner-up is ...02 (80.0).
# GenA sp2's best genome is ...03 (99.0); its runner-up is ...04 (70.0).
_RECORDS = [
    _rec("GCA_000000001.1", "GenA sp1", "99.0"),
    _rec("GCA_000000002.1", "GenA sp1", "80.0"),
    _rec("GCA_000000003.1", "GenA sp2", "99.0"),
    _rec("GCA_000000004.1", "GenA sp2", "70.0"),
]


@pytest.fixture
def asset(tmp_path):
    keys = ["ncbi_genbank_assembly_accession"] + list(RANKS) + _EXTRA_COLS
    cols = {k: [rec[k] for rec in _RECORDS] for k in keys}
    path = tmp_path / "gtdb.parquet"
    pq.write_table(
        pa.table({k: pa.array(v, type=pa.string()) for k, v in cols.items()}),
        str(path))
    return str(path)


def _select(asset, exclude=None, derep_rank="species"):
    return select_ref_genomes(asset, "gtdb", "GenA", derep_rank=derep_rank,
                              reps_only=False, exclude_cores=exclude)


class TestExclusionHappensBeforeDereplication:
    """
    The behaviour this feature exists for: an excluded genome leaves the candidate
    pool, so its dereplication group still contributes a genome.
    """

    def test_baseline_picks_the_best_genome_in_each_group(self, asset):
        selection = _select(asset)

        assert selection.accessions == ["GCA_000000001.1", "GCA_000000003.1"]
        assert selection.num_excluded == 0

    def test_excluding_a_groups_best_promotes_its_runner_up(self, asset):
        selection = _select(asset, exclude={"000000001"})

        # sp1 is still represented, by its second-best genome rather than not at all
        assert selection.accessions == ["GCA_000000002.1", "GCA_000000003.1"]

    def test_the_group_is_only_lost_when_every_member_is_excluded(self, asset):
        selection = _select(asset, exclude={"000000001", "000000002"})

        assert selection.accessions == ["GCA_000000003.1"]

    def test_excluded_count_is_candidates_removed_not_genomes_lost(self, asset):
        # two candidates removed, but the selection only shrinks by one group
        selection = _select(asset, exclude={"000000001", "000000002"})

        assert selection.num_excluded == 2

    def test_the_advisory_is_raised_on_the_selection(self, asset):
        selection = _select(asset, exclude={"000000001"})

        assert exclusion_warning(1) in selection.warnings

    def test_derep_off_simply_drops_the_listed_genomes(self, asset):
        selection = _select(asset, exclude={"000000001"}, derep_rank="off")

        assert selection.accessions == ["GCA_000000002.1", "GCA_000000003.1",
                                        "GCA_000000004.1"]
        assert selection.num_excluded == 1

    def test_no_exclusion_keeps_everything(self, asset):
        selection = _select(asset, exclude=None)

        assert len(selection.accessions) == 2
        assert selection.num_excluded == 0
        assert not any("--exclusion-list" in w for w in selection.warnings)


class TestCountsAgreeWithWhatAPullReturns:
    """
    The counts flags preview a pull, so they honor the same list. A count you have to
    mentally correct would be worse than no count.
    """

    def test_total_genome_count_drops_by_the_excluded_candidates(self, asset):
        assert count_genomes(asset, "gtdb") == 4
        assert count_genomes(asset, "gtdb", exclude_cores={"000000001"}) == 3

    def test_derep_size_is_unchanged_when_a_runner_up_can_step_in(self, asset):
        # exactly mirrors the selection: sp1 still contributes a genome
        base = derep_size(asset, "gtdb", "genus", "GenA", "species")
        excluded = derep_size(asset, "gtdb", "genus", "GenA", "species",
                              exclude_cores={"000000001"})

        assert base == 2
        assert excluded == 2

    def test_derep_size_drops_when_a_whole_group_is_excluded(self, asset):
        assert derep_size(asset, "gtdb", "genus", "GenA", "species",
                          exclude_cores={"000000001", "000000002"}) == 1

    def test_derep_size_matches_the_selection_it_previews(self, asset):
        for cores in (None, {"000000001"}, {"000000001", "000000002"}):
            previewed = derep_size(asset, "gtdb", "genus", "GenA", "species",
                                   exclude_cores=cores)
            actual = len(_select(asset, exclude=cores).accessions)

            assert previewed == actual, f"mismatch for {cores}"


class TestLoadingTheList:

    def test_entries_become_accession_cores(self, tmp_path):
        excl = _write(tmp_path / "excl.txt", ["GCF_000005845.2"])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_no_path_is_an_empty_set(self):
        assert load_exclusion_cores(None) == set()

    def test_blanks_and_comments_are_ignored(self, tmp_path):
        excl = _write(tmp_path / "excl.txt",
                      ["# a comment", "", "GCF_000005845.2", "   "])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_entries_naming_the_same_assembly_collapse(self, tmp_path):
        excl = _write(tmp_path / "excl.txt",
                      ["GCF_000005845.2", "GCA_000005845.1", "RS_GCF_000005845.3"])

        assert load_exclusion_cores(excl) == {"000005845"}

    def test_junk_entries_are_dropped_rather_than_kept_as_empty_cores(self, tmp_path):
        # otherwise a junk line would match every candidate that also fails to parse
        excl = _write(tmp_path / "excl.txt", ["na", "not-an-accession"])

        assert load_exclusion_cores(excl) == set()

    def test_read_preserves_first_seen_order_and_dedupes(self, tmp_path):
        excl = _write(tmp_path / "excl.txt",
                      ["GCF_2", "GCF_1", "GCF_2", "GCF_3"])

        assert read_exclusion_list(excl) == ["GCF_2", "GCF_1", "GCF_3"]


class TestMatchingIsOnAccessionCore:
    """
    GCA_/GCF_ pairing, GTDB's RS_/GB_ prefixes, and version suffixes all normalise
    away, so either member of a pair at any version excludes the assembly however the
    asset happens to spell it.
    """

    def _rows(self):
        return [{"acc": "GCF_000005845.1"}, {"acc": "GCA_111111111.1"}]

    def test_a_listed_gca_excludes_the_gcf_in_the_pool(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["GCA_000005845.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert [r["acc"] for r in kept] == ["GCA_111111111.1"]
        assert n == 1

    def test_a_gtdb_prefixed_entry_normalises_too(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["RS_GCF_000005845.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert [r["acc"] for r in kept] == ["GCA_111111111.1"]
        assert n == 1

    def test_a_different_assembly_is_left_alone(self, tmp_path):
        cores = load_exclusion_cores(
            _write(tmp_path / "excl.txt", ["GCA_000999999.2"]))

        kept, n = filter_rows_by_exclusion(self._rows(), "acc", cores)

        assert len(kept) == 2
        assert n == 0

    def test_junk_rows_are_not_matched_by_junk_entries(self, tmp_path):
        cores = load_exclusion_cores(_write(tmp_path / "excl.txt", ["na"]))
        rows = [{"acc": "some-weird-id"}, {"acc": "another-weird-id"}]

        kept, n = filter_rows_by_exclusion(rows, "acc", cores)

        assert len(kept) == 2
        assert n == 0


class TestTheThreeFilterFormsAgree:
    """
    Rows, plain accessions, and Arrow tables are three views of the same pool, so all
    three must normalise accessions identically.
    """

    ACCS = ["GCF_000005845.1", "GCA_111111111.1", "GCA_000005845.3"]

    def test_row_and_accession_forms_agree(self):
        cores = {"000005845"}
        rows = [{"acc": a} for a in self.ACCS]

        kept_rows, n_rows = filter_rows_by_exclusion(rows, "acc", cores)
        kept_accs, n_accs = filter_accessions_by_exclusion(self.ACCS, cores)

        assert [r["acc"] for r in kept_rows] == kept_accs
        assert n_rows == n_accs == 2

    def test_table_form_agrees(self):
        cores = {"000005845"}
        table = pa.table({"acc": pa.array(self.ACCS, type=pa.string())})

        kept_table, n_table = filter_table_by_exclusion(table, "acc", cores)
        kept_accs, n_accs = filter_accessions_by_exclusion(self.ACCS, cores)

        assert kept_table.column("acc").to_pylist() == kept_accs
        assert n_table == n_accs

    def test_an_empty_list_is_a_no_op_in_every_form(self):
        table = pa.table({"acc": pa.array(self.ACCS, type=pa.string())})

        assert filter_accessions_by_exclusion(self.ACCS, set()) == (self.ACCS, 0)
        assert filter_rows_by_exclusion([{"acc": a} for a in self.ACCS],
                                        "acc", None)[1] == 0
        assert filter_table_by_exclusion(table, "acc", None)[1] == 0


class TestExclusionWarningWording:

    def test_none_excluded_produces_no_line(self):
        assert exclusion_warning(0) is None

    def test_the_line_says_candidate_and_before_selection(self):
        line = exclusion_warning(3)

        assert "3 candidate genomes" in line
        assert "--exclusion-list" in line
        assert "before selection" in line

    def test_singular_and_plural_agree_with_their_verb(self):
        assert "1 candidate genome was removed" in exclusion_warning(1)
        assert "2 candidate genomes were removed" in exclusion_warning(2)


class TestHelpText:

    def test_help_names_the_taxon_flag_it_constrains(self):
        assert "`-t`" in exclusion_list_help("-t")


class TestDownloaderPreflightValidation:
    """
    `bit dl-ncbi-assemblies` rejects a list it could not act on, rather than finishing
    a run having silently ignored it. Mirrors GToTree's equivalent checks.
    """

    def _args(self, tmp_path, **kw):
        base = dict(wanted_accessions=None, target_taxon=None, exclusion_list=None,
                    dry_run=True, output_dir=str(tmp_path / "out"))
        base.update(kw)
        return types.SimpleNamespace(**base)

    def test_exclusion_list_without_a_target_taxon_exits(self, tmp_path):
        from bit.modules.ncbi.dl_ncbi_assemblies import preflight_checks

        accs = _write(tmp_path / "accs.txt", ["GCF_000005845.2"])
        excl = _write(tmp_path / "excl.txt", ["GCF_000005845.2"])
        # `-w` only: the list has nothing to act on, since named accessions are
        # downloaded as provided
        args = self._args(tmp_path, wanted_accessions=accs, exclusion_list=excl)

        with pytest.raises(SystemExit):
            preflight_checks(args)

    def test_a_missing_exclusion_file_exits(self, tmp_path):
        from bit.modules.ncbi.dl_ncbi_assemblies import preflight_checks

        args = self._args(tmp_path, target_taxon=["Alteromonas"],
                          exclusion_list=str(tmp_path / "nope.txt"))

        with pytest.raises(SystemExit):
            preflight_checks(args)
