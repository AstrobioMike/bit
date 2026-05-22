import pandas as pd # type: ignore


def filter_kofamscan_results(input_file, output_file):

    annot_dict = {}
    e_value_dict = {}

    with open(input_file, "r") as annots:

        for line in annots:

            if line.startswith("#"):
                continue

            line = line.lstrip("*").strip().split("\t")

            # adding gene ID if not present, to ensure all end up in final table
            if line[0] not in annot_dict:
                annot_dict[line[0]] = {"KO_ID": "NA", "KO_function": "NA"}

            # nothing there if no annotations for current gene, skipping
            if len(line) == 1:
                continue

            else:

                # only considering if its score is above the threshold
                # some, though very few, like K15869, don't have a threshold score due to having too few representatives, so if no threshold, just taking
                if line[2] == "" or float(line[3]) > float(line[2]):

                    # adding to e_value_dict if not represented already, adding annotation to annot_dict, and moving on
                    if line[0] not in e_value_dict:

                        annot_dict[line[0]] = {"KO_ID": line[1], "KO_function": line[5].strip('"')}
                        e_value_dict[line[0]] = line[4]
                        continue

                    else:

                        # replacing current annotation only if e-value is lower than current
                        if float(line[4]) < float(e_value_dict[line[0]]):

                            annot_dict[line[0]] = {"KO_ID": line[1], "KO_function": line[5].strip('"')}
                            e_value_dict[line[0]] = line[4]

    annot_tab = pd.DataFrame.from_dict(annot_dict, orient="index")
    annot_tab.reset_index(inplace=True)
    annot_tab.rename(columns={"index": "gene_ID"}, inplace=True)

    with open(output_file, "w") as out:
        out.write(annot_tab.to_csv(index=False, sep="\t"))
