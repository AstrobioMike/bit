from Bio import SeqIO # type: ignore


def genbank_to_fasta(input_genbank, output_fasta):
    with open(output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(input_genbank, "genbank"):
            out_fasta.write(f">{record.name}\n{record.seq}\n")


def genbank_to_AA_seqs(input_genbank, output_fasta):

    note_terms_to_exclude = ["frameshifted", "internal stop", "incomplete"] # dumping gene if noted as these in the "note" section of the call to keep only complete genes
    location_terms_to_exclude = ["join", "<", ">"] # dumping gene if "location" section contains any of these: "join" means the gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig

    with open(output_fasta, "w") as out_fasta:

        for record in SeqIO.parse(input_genbank, "genbank"):

            for feature in record.features:

                if feature.type == "CDS" and "translation" in feature.qualifiers:

                    location = str(feature.location)

                    if any(exclusion_term in location for exclusion_term in location_terms_to_exclude):
                        continue

                    if "note" in feature.qualifiers:
                        note = str(feature.qualifiers["note"][0])
                    else:
                        note = ""

                    if any(exclusion_term in note for exclusion_term in note_terms_to_exclude):
                        continue

                    if "transl_except" in feature.qualifiers:
                        continue

                    if "pseudo" in feature.qualifiers:
                        continue

                    if "locus_tag" in feature.qualifiers:
                        locus_tag = str(feature.qualifiers["locus_tag"][0])
                    else:
                        locus_tag = "No_locus_tag"

                    if "protein_id" in feature.qualifiers:
                        protein_id = str(feature.qualifiers["protein_id"][0])
                    else:
                        protein_id = "No_protein_id"

                    if "product" in feature.qualifiers:
                        product = str(feature.qualifiers["product"][0]).replace(" ", "-")
                    else:
                        product = "No_product"

                    out_fasta.write(f">{product}|{locus_tag}|{protein_id}\n{feature.qualifiers['translation'][0]}\n")
