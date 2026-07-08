#!/usr/bin/env bash

# this is meant to be done after a full conda install (or setting up all dependencies)
## making sure we are in the bit-dev conda environment
if [ "${CONDA_DEFAULT_ENV}" != "bit-dev" ]; then
    printf "\n    This should be run in the 'bit-dev' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    return 1 2>/dev/null || exit 1
fi

rm -rf build/ bit.egg-info/

pip install --no-build-isolation -e .

if command -v register-python-argcomplete >/dev/null 2>&1; then
    eval "$(register-python-argcomplete bit)"
fi

## if changing conda versions and wanting to install locally entirely (rather than using a prior official conda install of bit)
# conda build -c conda-forge -c bioconda conda-recipe/
# conda create -n bit-dev -c conda-forge -c bioconda --use-local bit
