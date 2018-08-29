#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "This script returns the total number of bases (or amino acids) in a fasta file.\n\n"
  printf "Usage:\n\t count-bases <input.fasta>\n\n"
  exit 1
fi

if [ -f $1 ]; then
  bases=$(grep -v ">" $1 | wc | awk '{print $3-$1}')
  echo -e "     ${GREEN}$bases bases in $1${NC}"

else
  echo -e "     ${RED}Input file not found :/${NC}"
  exit 1

fi
