#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "This script removes hardwraps from a fasta file.\n\n"
  printf "Usage:\n\t remove-wrap <input.fasta> [-k]\n\n"
  printf 'optional arguments:\n  -k, KEEP\t adding this flag will keep the original file\n\n'
  exit 1
fi

if [ -f $1 ]; then
  awk '!/^>/ { printf "%s", $0; n="\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' "$1" > temp

  if [ $2 ]; then

    if [ $2 == "-k" ]; then
      mv "$1" "$1".tmp
      mv temp "$1"
      echo -e "     ${GREEN}Fasta replaced with non-wrapping version, original stored with \".wrap\" appended. Cheers!${NC}"
    fi

  else
    mv temp "$1"
    echo -e "     ${GREEN}Original fasta replaced with non-wrapping version. Cheers!${NC}"

  fi
else
  echo -e "     ${RED}Input file not found :/${NC}"
  exit 1
fi
