#!/usr/bin/env python

import os
import argparse
import sys

parser = argparse.ArgumentParser(description = 'This script is for combining bracken output tables. It was modified\
                                              from the `combine_bracken_outputs.py` script provided by Jennifer\
                                              Lu (jlu26@jhmi.edu) that comes with bracken for use with the\
                                              `bit-combine-bracken-and-add-lineage` script. For version info,\
                                              run `bit-version`.')


required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input-files", metavar = "<FILE(s)>", nargs = "+", type = str, help = "space-delimited list of bracken output files", action = "store", required = True)
parser.add_argument("-n", "--sample-names", metavar = "<NAME(s)>", help = 'Sample names provided as a comma-delimited list (by default will use basename of input files)', action = "store", default = '')
parser.add_argument("-o", "--output-file", metavar = "<FILE>", help='Output file of combined tables (default: "combined-bracken.tsv")', action = "store", default = "combined-bracken.tsv")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

# setting up variables
sample_counts = {}
total_counts = {}
all_samples = []

# setting sample names and intializing counts
if len(args.sample_names) == 0:
    for file in args.input_files:
        curr_sample = os.path.basename(file)
        total_counts[curr_sample] = 0
        all_samples.append(curr_sample)

else:
    for curr_sample in args.sample_names.split(","):
        total_counts[curr_sample] = 0
        all_samples.append(curr_sample)


# working on each file
# initialize level variable
level = ''
# initializiing iterator for grabbing sample names
i = 0

for file in args.input_files:

    # storing current sample name
    curr_name = all_samples[i]

    # incrementing iterator
    i += 1

    with open(file) as f:
        # skipping header
        next(f)
        for line in f:
            [name, taxid, taxlvl, kreads, areads, estreads, frac] = line.strip().split("\t")
            estreads = int(estreads)

            # error checks
            if name not in sample_counts:
                sample_counts[name] = {}
                sample_counts[name][taxid] = {}
            elif taxid != list(sample_counts[name].keys())[0]:
                sys.exit("Taxonomy IDs not matching for species %s: (%s\t%s)" % (name, taxid, list(sample_counts[name].keys())[0]))
            if len(level) == 0:
                level = taxlvl
            elif level != taxlvl:
                sys.exit("Taxonomy level not matching between samples :(")

            # summing counts for current sample
            total_counts[curr_name] += estreads
            # adding read counts for that taxa for this sample to the dict holding all samples
            sample_counts[name][taxid][curr_name] = estreads


# opening output file
output_file = open(args.output_file, "w")

# writing header
output_file.write("name\ttax_id\ttax_level")
for name in all_samples:
    output_file.write("\t%s_num\t%s_frac" % (name, name))
output_file.write("\n")

# writing out each sample
for name in sample_counts:
    taxid = list(sample_counts[name].keys())[0]
    output_file.write("%s\t%s\t%s" % (name, taxid, level)) # seems like "level" variable is trusting the last thing was the same for all as it was for the last file's last line, probably true, but then not sure why we check. might return to this

    #Calculate and print information per sample
    for sample in all_samples:
        if sample in sample_counts[name][taxid]:
            num = sample_counts[name][taxid][sample]
            perc = float(num)/float(total_counts[sample])
            output_file.write("\t%i\t%0.5f" % (num, perc))

        # if sample doesn't have counts for this taxa, adding zeroes
        else:
            output_file.write("\t0\t0.00000")

    output_file.write("\n")

output_file.close()
