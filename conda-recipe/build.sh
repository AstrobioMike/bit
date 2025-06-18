#!/usr/bin/env bash

# setting up directories for stored stuff
mkdir -p ${PREFIX}/share/bit/gtdb-tax-info/
mkdir -p ${PREFIX}/share/bit/ncbi-tax-info/
mkdir -p ${PREFIX}/share/bit/go-dbs/
mkdir -p ${PREFIX}/share/bit/ncbi-assembly-summaries/


# adding placeholder files so the directories are not deleted on creation
touch ${PREFIX}/share/bit/gtdb-tax-info/GTDB-arc-and-bac-metadata.tsv
touch ${PREFIX}/share/bit/ncbi-tax-info/nodes.dmp
touch ${PREFIX}/share/bit/go-dbs/conda-placeholder
touch ${PREFIX}/share/bit/ncbi-assembly-summaries/date-retrieved.txt


# adding database paths to conda environment activation script
mkdir -p ${PREFIX}/etc/conda/activate.d/
echo 'export GTDB_DIR=${CONDA_PREFIX}/share/bit/gtdb-tax-info/' >> ${PREFIX}/etc/conda/activate.d/bit.sh
echo 'export TAXONKIT_DB=${CONDA_PREFIX}/share/bit/ncbi-tax-info/' >> ${PREFIX}/etc/conda/activate.d/bit.sh
echo 'export GO_DB_DIR=${CONDA_PREFIX}/share/bit/go-dbs/' >> ${PREFIX}/etc/conda/activate.d/bit.sh
echo 'export NCBI_assembly_data_dir=${CONDA_PREFIX}/share/bit/ncbi-assembly-summaries/' >> ${PREFIX}/etc/conda/activate.d/bit.sh

# running pip install of bit
"${PYTHON}" -m pip install . --no-deps --ignore-installed -vv
