#!/usr/bin/env bash

RED='\033[0;31m'
YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

# a function for providing help
print_help() {


    printf "\n                              ${YELLOW}HELP MENU${NC}"
    printf "\n  ${YELLOW}********************************************************************${NC}\n"

    printf "\n  Sometimes multiple SRA accessions comprise 1 sample. This script\n"
    printf "  is a helper to combine multiple SRA fastq files together.\n\n"
    printf "  It expects as input a tsv, with no header, where the first column\n"
    printf "  is the sample name, and the second column is the SRA accession, e.g.:\n\n"

    printf "              Sample-1\tSRR123456\n"
    printf "              Sample-1\tSRR123457\n"
    printf "              Sample-2\tSRR123458\n"
    printf "              Sample-3\tSRR123459\n"
    printf "              Sample-3\tSRR123460\n\n"

    printf "  It takes two positional arguments, the first being the tsv mapping file,\n"
    printf "  and the second being the path to the directory holding all the fastq files.\n\n"

    printf "    Ex. Usage:\n\t bash scripts/combine-sra-accessions.sh -i map.tsv -d fastq-files/ \n\n"

    printf "  Note that this is a simple bash script, it will work on one sample-set at a time,\n"
    printf "  and there is not much checked to catch human error on the input table.\n\n"

    printf "  By default it will remove the initial fastq files. Provide the '-k' flag if you want to\n"
    printf "  keep them.\n"

    printf "\n  ${YELLOW}********************************************************************${NC}\n\n"

    exit

}

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
    print_help
fi


########################################
### Setting up and parsing arguments ###
########################################
remove_original_fastqs="true"

while getopts ":i:d:k" args; do
    case "${args}"
    in
        i) map_file=$OPTARG;;
        d) fastq_dir=$OPTARG;;
        k) remove_original_fastqs="false";;
        \?) printf "\n  ${RED}Invalid argument: -${OPTARG}${NC}\n" 1>&2
            print_help
            ;;
        :)
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            print_help
            ;;
    esac
done

##################################################
## Making sure required arguments were provided ##
##################################################
if [ ! -n "${map_file}" ] || [ ! -n "${fastq_dir}" ]; then

    printf "\n    ${RED}ERROR${NC}: The required arguments were not provided. See help below.\n"
    print_help

fi


########################################
########### Pre-flight checks ##########
########################################


# check that the first positional argument is a file
if [ ! -f ${map_file} ]; then

    printf "\n    ${RED}ERROR${NC}: The file '${map_file}' does not exist.\n"
    print_help

fi

# check that the second positional argument is a directory
if [ ! -d ${fastq_dir} ]; then

    printf "\n    ${RED}ERROR${NC}: The directory '${fastq_dir}' does not exist.\n"
    print_help

fi

# checking input table has 2 columns
if [ $(head -n 1 ${map_file} | awk '{print NF}') -ne 2 ]; then

    printf "\n    ${RED}ERROR${NC}: The input table must have 2 columns. See help below.\n"
    print_help

fi


########################################
########## Getting to work #############
########################################

starting_dir=$(pwd)
path_to_map=$(realpath ${map_file})

# moving into directory to make it easier to run cat
cd ${fastq_dir}

printf "\n"

for sample in $(cut -f 1 ${path_to_map} | sort -u)
do

    printf "  Currently working on: ${sample} ...\r"

    target_R1s=$(grep ${sample} ${path_to_map} | cut -f 2 | sed 's/$/_R1.fastq.gz/' | tr '\n' ' ' | sed 's/ $//')
    target_R2s=$(grep ${sample} ${path_to_map} | cut -f 2 | sed 's/$/_R2.fastq.gz/' | tr '\n' ' ' | sed 's/ $//')


    if [ ${remove_original_fastqs} == "true" ]; then

        # checking if there are multiple, if so we cat them; if just one, we just mv/rename it
        if printf "${target_R1s}" | grep -q " "; then

            cat ${target_R1s} > ${sample}_R1.fastq.gz
            cat ${target_R2s} > ${sample}_R2.fastq.gz

            rm ${target_R1s}
            rm ${target_R2s}

        else

            mv ${target_R1s} ${sample}_R1.fastq.gz
            mv ${target_R2s} ${sample}_R2.fastq.gz

        fi

    else

        cat ${target_R1s} > ${sample}_R1.fastq.gz
        cat ${target_R2s} > ${sample}_R2.fastq.gz

    fi

done

# moving back to initial dir
cd ${starting_dir}

printf "\n\n                     ${GREEN}DONE!${NC}\n\n"
printf "  ${YELLOW}The combined fastq files are in the directory: ${fastq_dir}${NC}\n\n"

if [ ${remove_original_fastqs} == "true" ]; then

    printf "  ${YELLOW}The original fastq files were removed because the '-k' flag was not provided.${NC}\n\n"

else

    printf "  ${YELLOW}Note that the original fastq files were left because the '-k' flag was provided.${NC}\n\n"

fi
