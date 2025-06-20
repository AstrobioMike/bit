#!/bin/bash
set -e

# this script as written expects the base conda directory to be located at ~/miniconda3/

# auto-works on conda-recipe/ dir
# example usage: bash update-conda-package.sh

##########################################################################
### if setting up on a new system for the first time, these are needed ###
##########################################################################

## needed for build
# conda install conda-build

## need to go through this process for uploads to anaconda to work as performed below
# conda install anaconda-client
# anaconda login

##########################################################################

### some pre-flight checks ###

recipe_dir="conda-recipe"

## making sure it is a directory here
if [ ! -d "${recipe_dir}" ]; then
    printf "\n    'conda-recipe/' needs to be a directory here.\n"
    printf "    You know this, Mike...\n\n"
    exit 1
fi

## making sure we are in the base conda environment
if [ "${CONDA_DEFAULT_ENV}" != "base" ]; then
    printf "\n    This should be run in the 'base' conda environment..\n"
    printf "    You know this, Mike...\n\n"
    exit 1
fi

### getting started ###

## getting current version from setup.py
version=$(grep "version" setup.py | cut -f 2 -d '=' | tr -d '",')
## getting current build from meta.yaml
build=$(grep -A 1 ^build ${recipe_dir}/meta.yaml | grep number | cut -f 2 -d ":" | tr -s " " "\t" | cut -f 2)


## conda way
conda-build -c conda-forge -c bioconda -c defaults ${recipe_dir}/

program=$(grep "name" setup.py | cut -f 2 -d '=' | tr -d '",')

## converting to other platforms
conda convert --platform linux-64 ~/miniconda3/conda-bld/osx-64/${program}-${version}-*_${build}.tar.bz2 -o ~/miniconda3/conda-bld/
# conda convert --platform osx-arm64 ~/miniconda3/conda-bld/osx-64/${program}-${version}-*_${build}.tar.bz2 -o ~/miniconda3/conda-bld/
    # i am leaving off the new mac chip way because most of the dependencies still don't have arm64 versions


## uploading to anaconda
anaconda upload ~/miniconda3/conda-bld/osx-64/${program}-${version}-*_${build}.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-64/${program}-${version}-*_${build}.tar.bz2
# anaconda upload ~/miniconda3/conda-bld/osx-arm64/${program}-${version}-*_${build}.tar.bz2
    # i am leaving off the new mac chip way because most of the dependencies still don't have arm64 versions
