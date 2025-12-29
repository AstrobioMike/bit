# Change Log

<!-- 
## v (dd-mmm-yyyy)

### Added

### Changed

### Fixed

---

-->


## v1.13.5 (NOT RELEASED YET)

### Added

### Changed

### Fixed
- `bit-assemble` now properly filters spades-assembled contigs based on user-specific min-contig length

---



## v1.13.4 (20-Nov-2025)

### Changed
- updates to `bit-get-cov-stats`
  - can start from bam file now in addition to mosdepth per-base.bed.gz (will generate the mosdepth output if starting from bam)
  - modularized, integration test added
- updates to `bit-check-for-fastq-dup-headers`
  - autodetect gzipped or not
  - modularized, test added

---

## v1.13.3 (29-Sep-2025)

### Added
- more test coverage of `bit-ez-screen`
- unit tests for `bit-gen-kraken2-tax-plots`, `bit-kraken2-to-taxon-summaries`, and `bit-calc-variation-in-msa`
- integration test for `bit-cov-analyzer`

### Changed
- modularized `bit-calc-variation-in-msa`
- updates to `bit-gen-kraken2-tax-plots`
  - modularized
  - appropriately adds domain letter to plots from GTDB tax kraken2 reports now

---

## v1.13.2 (22-Aug-2025)

### Changed
- updates to `bit-kraken2-to-taxon-summaries`
  - modularized
  - no longer takes the larger kraken.out file, it now works off of the kraken.report
  - no longer works based on taxid lookup, it now works based on the taxonomy in the kraken report
    - this means it will exactly match the taxonomy used in the kraken2 db, and not swap anything if taxids or rank names changed
    - this also means it now works with GTDB-kraken2-db produced reports (which previously would not work with the standard taxid lookup method)

---

## v1.13.1 (06-Aug-2025)

### Changed
- modularized `bit-filter-seqs-by-length` and added a test for it
  - also changed the name to `bit-filter-fasta-by-length`, though i'm retaining a stub for the old name so it still works when called that way too
- modularized `bit-summarize-column` and added tests
- added a temp fix for the latest megahit osx-64 build not working (see smk/envs/assemble-osx-64.yaml; if that's gone it was no longer needed and removed in the future)

### Fixed 
- fix to `bit-summarize-column` when standard input is only one column (was erroring out before)

---

## v1.13.0 (22-Jul-2025)

### Added
- added `bit-assemble`
  - command-line wrapper for an assembly workflow with optional qc and digital normalization
- added more integration tests

### Changed
- modularized `bit-gen-reads`

---

## v1.12.3 (26-Jun-2025)

### Changed
- modifications to `bit-cov-analyzer`
  - default --min-region-length set to 500 to help reduce overwhelming output and focus on larger regions
  - added a column for "zero_cov_bases" to output low- and high-coverage region tsvs

---

## v1.12.2 (25-Jun-2025)

### Changed
- modifications to `bit-cov-analyzer`
  - added N filtering
    - if an an identified low-coverage region is more than 50% Ns, it won't be reported in the output low-coverage table
  - added `--min-region-length` parameter, though the default is 0 (so smallest window size)
  - modularized this script

---

## v1.12.1 (25-Jun-2025)

### Changed
- adding seed option to `bit-gen-reads`

---

## v1.12.0 (24-Jun-2025)

### Added
- added a reads mode to `bit-ez-screen` 

### Changed
- `bit-ez-screen` now has subcommands for assembly vs reads modes

---

## v1.11.1 (20-Jun-2025)

### Fixed 
- fixed read-header formatting for paired read from `bit-gen-reads`
- fixed conda recipe meta.yaml and update-conda-package.sh

---

## v1.11.0 (18-Jun-2025)

### Changed
- started making bit a python package (to facilitate sharing functions, formatting, etc. overtime and moving forward)
- added conda recipe here instead of being maintained elsewhere

---

## v1.10.9 (04-Jun-2025)

### Fixed 
- fixed `bit-get-workflow` to be able to pull from all prior released workflows instead of just recent ones

---

> Previous version changes are only tracked on the [releases page](https://github.com/AstrobioMike/bit/releases).

---