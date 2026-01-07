#!/usr/bin/env bash

# this is meant to be done after a full conda install (or setting up all dependencies)
## making sure we are in the bit-dev conda environment
if [ "${CONDA_DEFAULT_ENV}" != "bit-dev" ]; then
    printf "\n    This should be run in the 'bit-dev' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    exit 1
fi

pip install -e .

BIN_DIR=$(dirname $(which python))

printf "\n\n  Linking scripts to ${BIN_DIR}/ for dev...\n\n"

for script in bit/scripts/*; do
    ln -sf "$(realpath "${script}")" "${BIN_DIR}/$(basename "${script}")"
done
