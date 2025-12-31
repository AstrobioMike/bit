import re

def get_rank(lineage, prefix):

    if lineage.startswith(prefix):

        curr_rank = lineage.split(";")[0].replace(prefix, "", 1)

        lineage = re.sub(f"^{prefix}{curr_rank};", "", lineage)

    else:

        curr_rank = "NA"

    return(lineage, curr_rank)


def convert_lineage_to_tsv(input_tsv, output_tsv, make_taxid):

    with open(output_tsv, "w") as out_tab:

        if make_taxid:
            out_tab.write("seq_ID\ttaxid\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
        else:
            out_tab.write("seq_ID\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")

        for line in open(input_tsv):

            line = line.strip().split("\t")
            ID = line[0]

            # this is if there is no second column (no lineage)
            if len(line) == 1:

                out_line = f"{ID}\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

            else:

                lineage = line[1]

                # removing "root" if that's at the start
                if lineage.startswith("root;"):
                    lineage = re.sub("^root;", "", lineage)

                # getting all ranks present, setting to NA if not
                lineage, t_domain = get_rank(lineage, "d__")
                lineage, t_phylum = get_rank(lineage, "p__")
                lineage, t_class = get_rank(lineage, "c__")
                lineage, t_order = get_rank(lineage, "o__")
                lineage, t_family = get_rank(lineage, "f__")
                lineage, t_genus = get_rank(lineage, "g__")
                lineage, t_species = get_rank(lineage, "s__")

                if make_taxid:

                    taxid_string = "_".join([t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species]).replace(" ", "_")
                    out_line = "\t".join([ID, taxid_string, t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species])

                else:

                    out_line = "\t".join([ID, t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species])

            out_tab.write(out_line + "\n")
