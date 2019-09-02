#!/usr/bin/env bash

# setting colors to use
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

my_ext=$2
ext=$3

assembly=$(echo "$1" | cut -f 1)
downloaded_accession=$(echo "$1" | cut -f 2)

# storing and building links
base_link=$(echo "$1" | cut -f 9)


# checking link was actually present (sometimes, very rarely, it is not there)
# if not there, attempting to build ourselves
if [ $base_link == "na" ] || [ -z $base_link ]; then

    p1=$(printf "ftp://ftp.ncbi.nlm.nih.gov/genomes/all")

    # checking if GCF or GCA
    if [[ $assembly == "GCF"* ]]; then 
        p2="GCF"
    else
        p2="GCA"
    fi
    
    p3=$(echo $assembly | cut -f 2 -d "_" | cut -c 1-3)
    p4=$(echo $assembly | cut -f 2 -d "_" | cut -c 4-6)
    p5=$(echo $assembly | cut -f 2 -d "_" | cut -c 7-9)

    ass_name=$(echo "$1" | cut -f 3)
    end_path=$(paste -d "_" <(echo "$assembly") <(echo "$ass_name"))

    base_link=$(paste -d "/" <(echo "$p1") <(echo "$p2") <(echo "$p3") <(echo "$p4") <(echo "$p5") <(echo "$end_path"))

else

    end_path=$(basename $base_link)

fi

curl --silent --retry 10 -o ${assembly}${my_ext}.gz "${base_link}/${end_path}${ext}.gz"

if [ -s ${assembly}${my_ext}.gz ]; then

    printf "\r\t  Successfully downloaded: $assembly"

else

    printf "\n     ${RED}******************************* ${NC}NOTICE ${RED}*******************************${NC}  \n"
    printf "\t\t  $assembly's $format file didn't download successfully.\n\n"
    printf "\t\t  Written to \"NCBI-accessions-not-downloaded.txt\".\n"
    printf "     ${RED}********************************************************************** ${NC}\n\n"

    echo ${assembly} >> NCBI-accessions-not-downloaded.txt

fi