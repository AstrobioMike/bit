#!/usr/bin/env python

import sys
import argparse
import os

parser = argparse.ArgumentParser(description = 'This script is for parsing NCBI\'s assembly summary file down\
                                              to the provided accessions. It is used by the `bit-dl-ncbi-assemblies`\
                                              script. For version info, run `bit-version`.')

required = parser.add_argument_group('required arguments')

required.add_argument("-a", "--assembly-summary", metavar = "<FILE>", help = "NCBI's assembly summary file", action = "store", dest = "all_assemblies", required = True)
required.add_argument("-w", "--wanted-accessions", metavar = "<FILE>", help = "Single-column file with wanted accessions", action = "store", dest = "wanted_accs", required = True)
parser.add_argument("-o", "--output-file", help = 'Output file of wanted summary info only (default: "wanted.tsv")', action = "store", default = "wanted.tsv")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

wanted_dict = {}

with open(args.wanted_accs, "r") as wanted_accs:

    for line in wanted_accs:
        root_acc = line.strip().split(".")[0]
        wanted_dict[str(root_acc)] = line.strip()

out_file = open(args.output_file, "w")

with open(args.all_assemblies) as assemblies:

    for line in assemblies:
        line = line.split("\t")

        if line[0].split(".")[0] in wanted_dict:
            dl_acc = str(line[0])

            if not dl_acc:
                dl_acc = "NA"

            ass_name = str(line[15])
            if not ass_name:
                ass_name = "NA"

            taxid = str(line[5])
            if not taxid:
                taxid = "NA"

            org_name = str(line[7])
            if not org_name:
                org_name = "NA"

            infra_name = str(line[8])
            if not infra_name:
                infra_name = "NA"

            version_status = str(line[10])
            if not version_status:
                version_status = "NA"

            ass_level = str(line[11])
            if not ass_level:
                ass_level = "NA"

            ftp_path = str(line[19])
            if not ftp_path:
                ftp_path = "NA"

            out_file.write(str(wanted_dict[str(line[0].split(".")[0])]) + "\t" + str(dl_acc) + "\t"  + str(ass_name) + "\t" + str(taxid) + "\t" + str(org_name) + "\t" + str(infra_name) + "\t" + str(version_status) + "\t" + str(ass_level) + "\t" + str(ftp_path) + "\n")
