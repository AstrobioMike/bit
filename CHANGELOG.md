# Change Log

<!-- 
## v (dd-mmm-yyyy)

### Added

### Changed

### Fixed

---

-->

## v1.16.0 (NOT YET RELEASED)
### Added
- setup for argcomplete during conda install for commands with subcommands
- added flag `--skip-read-pids` to `bit-get-cov-stats` so user can save time if they don't want that info


### Changed
- removed `bit-gff-to-anvio`, that is now just stored as a gist here: https://gist.github.com/AstrobioMike/45e10adb3eaceb338a7eb10f49038355
- removed `bit-split-multifasta`, that is now just stored as a gist here: https://gist.github.com/AstrobioMike/28cae086241bc0d68e05c215acd91a69
- `bit-get-cov-stats`
  - now just called `bit-cov-stats`
  - output *-per-ref.tsv now includes "length" and "num_contigs" columns


---


## v1.15.0 (14-Apr-2026)

### Added
- `bit-mutate-seqs` now has a tunable parameter for transitions/transversions
- `bit-extract-seqs`
  - added `by-headers` subcommand (**this replaces `bit-parse-fasta-by-headers`**)
  - `by-coords` subcommand **replaces `bit-extract-seqs-by-coords`**
- `bit-genbank` which has subcommmands for all genbank helpers
  - `to-fasta` replaces `bit-genbank-to-fasta`
  - `to-AA-seqs` replaces `bit-genbank-to-AA-seqs`
  - `to-cds-seqs` 
  - `to-cds-tsv` replaces `bit-genbank-to-cds-table`


### Changed
- removed `bit-figshare-upload`, that is now just stored as a gist here: https://gist.github.com/AstrobioMike/9f86931747357ae8949f715145c5eec4
- removed `bit-locus-clean-slate`, that is now just stored as a gist here: https://gist.github.com/AstrobioMike/17324583568b051cd3b190cf4767b524
- moved many imports from top of cli modules until they are needed to increase snappiness when just doing things like trying to view the help menu
- `bit-summarize-column` no longer has an `-i` flag, it expects the input table/file to be given as a positional argument
- `bit-count-bases` can now take STDIN
- `bit-normalize-table` modularized and tests added
- `bit-parse-fasta-by-headers` has been removed and replaced with `bit-extract-seqs by-headers`
- `bit-extract-seqs-by-coords` has been removed and replaced with `bit-extract-seqs by-coords`
- `bit-filter-seqs-by-length` has been renamed to `bit-filter-fasta-by-length`
- `bit-genbank-to-fasta` has been removed and replaced with `bit-genbank to-fasta`
- `bit-genbank-to-AA-seqs` has been removed and replaced with `bit-genbank to-AA-seqs`
- `bit-genbank-to-cds-table` has been removed and replaced with `bit-genbank to-cds-tsv`
- `bit-genbank-to-clean-slate` has been removed and replaced with `bit-genbank clean-slate`
- moved `bit-GL-combine-contig-tax-tables` and `bit-GL-combine-KO-and-tax-tables` from the primary bit package into the metagenomics-wf scripts
- `bit-get-test-data metagenome` pulls a smaller, super simple, single-sample dataset I'm hosting on [github](https://github.com/AstrobioMike/test-data/releases/tag/test-metagenomics-reads-v1) now rather than figshare

---

## v1.14.0 (07-Apr-2026)

### Added
- to `bit-cov-analyzer`
  - progress updates while running
  - zero-coverage region outputs now also generated
  - output region tsvs now have a "low_complexity" column that holds True or False
    - this is based on:
      -  low_complexity = True if: unique 3-mers / all-possible-3mers <= 0.4
- `bit-extract-seqs`
  - enables pulling out target seqs from a fasta by bed file or specified primers via subcommands
- to `bit-gen-reads`
  - `--fragment-size-range` option added, defaults to 10% of fragment size

### Changed
- `bit-cov-analyzer`
  - `-s | --sliding-window-size` changed to `-w | --window-size`, and `-S | --step-size` changed to `-s | --step-size` (lower-case)
  - default window size change from 50 to 100, and default step size changed from 10 to 20
  - drastic improvements to efficiency when working with large genomes (e.g., 3GB)
  - histogram of coverages no longer plotted by default, only done now when adding the `--write-window-stats` flag
  - no longer produces window-coverage-overview.txt as all of that info is captured within window-coverage-overview.tsv
- `bit-get-mapped-reads-pid`
  - minor improvements to efficiency
- `bit-get-cov-stats`
  - improvements to efficiency
  - now also reports median percent id of mapped reads per ref and per contig (when provided an input bam file)
- `bit-summarize-assembly`
  - adds commas when printing stats to terminal for readability
- `bit-extract-seqs-by-coords` is now combined into `bit-extract-seqs`
- `bit-gen-reads`
  - now has a `--fragment-size-range` that defaults to 10% of fragment size
  - by default will not include regions with Ns in generated reads, add `--include-Ns` to allow that

---

## v1.13.15 (13-Mar-2026)

### Changed
- `bit-assemble`
  - the threads parameter is now passed to bbnorm and fastp (if run) in addition to the assemblers
- `bit-gen-reads`
  - `--type long` will no longer preferentially start reads at position 0 if the requested read size is larger than the contig; now it will start randomly and just produce a read that ends where the contig ends (unless `--circularize` is added)
- `bit-cov-analyzer`
  - no longer writes out individual window stats by default (to save spacetime), it needs to be turned on with the `--write-window-stats` now if wanted

---

## v1.13.14 (13-Mar-2026)

### Fixed
- `bit-gen-reads` previously may have by chance created reads with identical headers (since only coordinates were being added), now there is also a counter to prevent this

---

## v1.13.13 (12-Mar-2026)

### Added
- to `bit-gen-reads`
  - added single-end and long-read capabilities (through `--type` argument now, see *Changed* below)
  - single can be used up to any size, but if specifying `--long`, it will also generate reads with lengths spanning a range around the specified read size (50% by default)
- to `bit-calc-variation-in-msa`
  - 3Di as an option for `--type`

### Changed
- `bit-gen-reads`
  - now has `--type` flag for paired-end, single-end, or long (paired-end still by default)
  - did more work than it's worth to ensure the *exact* number of requested reads are always returned
- `bit-calc-variation-in-msa`
  - `--gaps-treatment` changed to "include" by default

---

## v1.13.12 (06-Mar-2026)

### Added
- `bit-add-insertion`

---


## v1.13.11 (03-Mar-2026)

### Changed
- `bit-update-ncbi-taxonomy` replaced with `get-ncbi-tax-data` (prior still retained for now)
- dropped `bit-calc`
  - if you are the one other person that ever used this and you want it back, you can add this to your ~/.bashrc: `bit-calc () { awk "BEGIN { print $1 }"; }` :)
- modified `bit-colnames` to try to autodetect delimeter

### Fixed
- added back in setup.py glob portion needed for scripts not fully integrated into python-packaging yet

---

## v1.13.10 (27-Feb-2026)

### Changed
- `bit-get-cov-stats`
  - the `--include-non-primary` flag now in addition to calculating percent ID including supplemental and secondary alignments also runs mosdepth with `--flag 1540`
- `bit-dl-ncbi-assemblies` 
  - in python now instead of bash (i hope this doesn't hinder performance too much...)
  - default concurrent downloads is 10 now instead of 1
  - default format is fasta now instead of gbk
  - downloads only happen in http now, no more ftp, so the -P flag to specify http has been removed
  - added optional output dir
- no longer keeping stubs in scripts/, instead keeping a ton of entry points in pyproject.toml
- `bit-filter-seqs-by-length` renamed to `bit-filter-fasta-by-length` to be more specific (prior retained for now)


---

## v1.13.9 (25-Feb-2026)

### Added
- `bit-get-cov-stats` by default now produces per-contig level info also (can be shut off with `--skip-per-contig`)

### Changed
- `bit-get-cov-stats`
  - the original ref-based output file is now called \<output-prefix\>-per-ref.tsv (changed from \<output-prefix\>.tsv)
  - outputs include median coverage in addition to mean
  - for speed (and consistency with expectations of known most-frequent users), when `bit-get-cov-stats` runs mosdepth, it uses the `-x | --fast-mode` flag now
  - added progress bar when parsing coverage info
- `bit-assemble`
  - re-arranging of help menu
  - memory setting now passable to spades too
- `report_message` function from modules.general slightly altered
  - this is more a note to myself for if/when i see weird things in terminal-printing format show up later
- general help-menu formatting

---

## v1.13.8 (05-Feb-2026)

### Added
- `--circularize` option added to `bit-gen-reads`

---

## v1.13.7 (20-Jan-2026)

### Changed
- `bit-get-cov-stats` now also reports mean percent ID of mapped reads for each input reference when the input includes a bam file (leveraging `bit-get-mapped-reads-pid`)

---

## v1.13.6 (16-Jan-2026)

### Added
- `bit-get-mapped-reads-pid` to pull out percent-identity information of mapped reads from an input bam
  - for each mapped read:
    - calculated percent ID = (full_aligned_length - NM) / full_aligned_length * 100
      - where full_aligned_length = Matches + Mismatches + Insertions + Deletions
- added a 'genome' option to `bit-get-test-data` that pulls an E. coli genome

### Changed
- modularized and added tests for `bit-get-workflow`, `bit-get-test-data`, `bit-dedupe-fasta-headers`, `bit-fasta-to-genbank`, and `bit-fasta-to-bed`
- added more tests to `bit-gen-reads`
- improved coverage on some other modules with more unit tests
- moved more setup info into pyproject.toml, but retained minimal setup.py to be able to glob because bit has a lot of separate scripts/entry points

---


## v1.13.5 (31-Dec-2025)

### Added

### Changed
- `bit-count-bases-per-seq` has been removed with its function combined into `bit-count-bases`
  - if input fasta has one sequence, it prints the length to the terminal; if it has 2 or more, it will print out summary stats; in either of the two prior cases, if an output file is specified, the program will additionally write lengths of all sequences to that specified file
  - modularized, test added
- modularized and added unit tests for `bit-lineage-to-tsv`
- `bit-mutate-seqs`
  - `--seed-for-randomization` long-parameter shortened to just `--seed`
  - modularized and test added

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
