#!/usr/bin/env bash

# this is meant to be done after a full conda install (or setting up all dependencies)
## making sure we are in the bit-dev conda environment
if [ "${CONDA_DEFAULT_ENV}" != "bit-dev" ]; then
    printf "\n    This should be run in the 'bit-dev' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    exit 1
fi

rm -rf build/ bit.egg-info/

BIN_DIR=$(dirname $(which python))

# removing stale symlinks that point into bit/scripts/
# (these conflict with pyproject.toml entry points)
for f in ${BIN_DIR}/*; do
    target=$(readlink "$f" 2>/dev/null)
    if [[ "$target" == *"bit/scripts"* ]]; then
        rm -f "$f"
    fi
done

pip install --no-build-isolation -e .

printf "\n\n  Linking unported scripts to ${BIN_DIR}/ for dev...\n\n\n"

for script in bit/scripts/*; do
    name=$(basename "${script}")
    # skip if pip already created an entry point for this name
    if [ -f "${BIN_DIR}/${name}" ]; then
        continue
    fi
    ln -sf "$(realpath "${script}")" "${BIN_DIR}/${name}"
done
