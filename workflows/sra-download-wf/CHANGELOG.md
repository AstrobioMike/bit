# SRA-download workflow change log

## 1.1.0
- workflow can now also handle sra objects that hold single-end data

## 1.0.1
- updates to `scripts/combine-sra-accessions.sh`
  - more efficient now by not cat'ing if there is only one SRR for a sample
  - default is to remove original files now, and `-k` needs to be added in order to keep them

## 1.0.0
- initial workflow release
