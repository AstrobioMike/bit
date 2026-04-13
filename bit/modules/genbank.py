from Bio import SeqIO # type: ignore
import pandas as pd # type: ignore


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

                    if "transl_except" in feature.qualifiers or "pseudo" in feature.qualifiers:
                        continue


                    locus_tag = feature.qualifiers.get("locus_tag", ["No_locus_tag"])[0]
                    gene = feature.qualifiers.get("gene", ["No_gene"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["No_protein_id"])[0]
                    product = feature.qualifiers.get("product", ["No_product"])[0].replace(" ", "-")

                    out_fasta.write(f">{locus_tag}|{protein_id}|{gene}|{product}\n{feature.qualifiers['translation'][0]}\n")


def genbank_to_cds_seqs(input_genbank, output_fasta):

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

                    if "transl_except" in feature.qualifiers or "pseudo" in feature.qualifiers:
                        continue

                    locus_tag = feature.qualifiers.get("locus_tag", ["No_locus_tag"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["No_protein_id"])[0]
                    gene = feature.qualifiers.get("gene", ["No_gene"])[0]
                    product = feature.qualifiers.get("product", ["No_product"])[0].replace(" ", "-")
                    seq = feature.extract(record.seq)

                    out_fasta.write(f">{locus_tag}|{protein_id}|{gene}|{product}\n{seq}\n")


def genbank_to_cds_tsv(input_genbank, output_tsv):

    cds_dataframe = parse_genbank_cds_to_dataframe(input_genbank)

    cds_dataframe.to_csv(output_tsv, sep="\t", index=False)


def parse_genbank_cds_to_dataframe(genbank_path):

    cds_entries = []

    with open(genbank_path, "r") as handle:

        for record in SeqIO.parse(handle, "genbank"):

            for feature in record.features:
                if feature.type == "CDS":

                    gene = feature.qualifiers.get("gene", ["NA"])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", ["NA"])[0]
                    product = feature.qualifiers.get("product", ["NA"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["NA"])[0]

                    cds_entries.append({
                        "gene": gene,
                        "protein_id": protein_id,
                        "locus_tag": locus_tag,
                        "product": product,
                    })

    cds_df = pd.DataFrame(cds_entries)

    return cds_df
