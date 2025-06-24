#!/usr/bin/env bash

# this is meant to be done after a full conda install (or setting up all dependencies)
# and after a `pip install -e .`, as that doesn't properly link scripts

BIN_DIR=$(dirname $(which python))

printf "\n  Linking scripts to ${BIN_DIR}/ for dev...\n\n"

for script in bit/scripts/*; do
    ln -sf "$(realpath "${script}")" "${BIN_DIR}/$(basename "${script}")"
done
