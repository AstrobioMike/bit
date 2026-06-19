# Metagenomics workflow change log

## 1.0.4
- removed f-strings and made shell blocks raw/r so Snakefile was compatiable with bit env python 3.12

## 1.0.3
- CAT was having trouble downloading its DB, so added `-k` to curl command
- moved `bit-GL-combine-contig-tax-tables` and `bit-GL-combine-KO-and-tax-tables` from bit package into workflow scripts

## 1.0.2
- pinned specific version of diamond (2.0.6) to the CAT environment

## 1.0.1
- can optionally skip binning and MAG recovery and characterization with new option in config.yaml, "perform_binning_and_MAG_recovery"

## 1.0.0
- initial workflow release
