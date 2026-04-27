#!/usr/bin/env python

import pandas as pd # type: ignore


def combine_kraken2_tax_summaries(all_samples, output_file):

    taxid_dict = {}

    # building final table
    building_tab = pd.DataFrame(columns=["taxid"])

    ## working on each file
    for sample_key in all_samples:

        # reading current file into pandas dataframe
        curr_tab = pd.read_csv(sample_key, sep="\t")

        # adding to building taxid dictionary
        for row in curr_tab.itertuples():
            if row[1] not in taxid_dict:
                taxid_dict[row[1]] = {'domain': row[2], 'phylum': row[3], 'class': row[4], 'order': row[5], 'family': row[6], 'genus': row[7], 'species': row[8]}

        # trimming down current table
        curr_sub_tab = curr_tab[["taxid", "read_counts", "percent_of_reads"]]

        # and changing count and percent column names to match sample ID
        curr_sub_tab.columns = ['taxid', str(all_samples[sample_key]) + "_read_counts", str(all_samples[sample_key]) + "_perc_of_reads"]

        # merging with master tab on GO_term
        building_tab = building_tab.merge(curr_sub_tab, on="taxid", how="outer")

    ## replacing NAs with 0s
    building_tab = building_tab.fillna(0)

    ## making taxid dictionary into dataframe and merging into final table
    taxid_df = pd.DataFrame.from_dict(taxid_dict, orient="index")

    taxid_df.reset_index(inplace=True)
    taxid_df.rename(columns = {'index': 'taxid'}, inplace=True)

    final_tab = taxid_df.merge(building_tab, on="taxid", how="outer")

    final_tab.sort_values(by=["taxid"], inplace=True)

    final_tab = final_tab.fillna("NA")

    with open(output_file, "w") as out:
        out.write(final_tab.to_csv(index=False, sep="\t"))
