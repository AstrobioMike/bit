def parse_ncbi_assembly_summary(assembly_summary_file, run_data):

    wanted_dict = {}

    for acc in run_data.wanted_accs:
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    FORMAT_EXTENSIONS = {
        "genbank":     ("_genomic.gbff.gz",            ".gb.gz"),
        "fasta":       ("_genomic.fna.gz",             ".fa.gz"),
        "protein":     ("_protein.faa.gz",             ".faa.gz"),
        "gff":         ("_genomic.gff.gz",             ".gff.gz"),
        "nt_cds":      ("_cds_from_genomic.fna.gz",    "_cds_from_genomic.fna.gz"),
        "feature_tab": ("_feature_table.txt.gz",       ".tsv.gz"),
        "report":      ("_assembly_report.txt",         "_assembly_report.txt"),
        "stats":       ("_assembly_stats.txt",          "_assembly_stats.txt"),
    }

    with open(run_data.ncbi_sub_table_path, "w") as out_file:

        cols = ["input_accession", "found_accession", "assembly_name", "taxid", "organism_name", "infraspecific_name", "version_status", "assembly_level", "http_base_link"]
        if run_data.wanted_format:
            cols.extend(["target_link", "local_destination"])
        out_file.write("\t".join(cols) + "\n")

        with open(assembly_summary_file, "r") as assemblies:
            for line in assemblies:
                fields = line.strip().split("\t")
                if not fields:
                    continue
                root = fields[0].split(".")[0]
                if root in wanted_dict:
                    found.add(wanted_dict[root])

                    dl_acc = fields[0].strip() if fields[0].strip() else "NA"
                    assembly_name = fields[15].strip() if len(fields) > 15 and fields[15].strip() else "NA"
                    taxid = fields[5].strip() if len(fields) > 5 and fields[5].strip() else "NA"
                    org_name = fields[7].strip() if len(fields) > 7 and fields[7].strip() else "NA"
                    infra_name = fields[8].strip() if len(fields) > 8 and fields[8].strip() else "NA"
                    version_status = fields[10].strip() if len(fields) > 10 and fields[10].strip() else "NA"
                    assembly_level = fields[11].strip() if len(fields) > 11 and fields[11].strip() else "NA"

                    http_path = build_base_link(dl_acc, assembly_name)

                    out_line = "\t".join([
                        wanted_dict[root],
                        dl_acc,
                        assembly_name,
                        taxid,
                        org_name,
                        infra_name,
                        version_status,
                        assembly_level,
                        http_path
                    ])
                    if run_data.wanted_format:
                        ncbi_ext, local_ext = FORMAT_EXTENSIONS[run_data.wanted_format]
                        target_link = f"{http_path}{dl_acc}_{assembly_name}{ncbi_ext}"
                        local_path = f"{run_data.output_dir}/{dl_acc}{local_ext}"
                        out_line += "\t" + target_link + "\t" + local_path

                    out_line += "\n"

                    out_file.write(out_line)

    not_found = set(run_data.wanted_accs) - found

    if len(not_found) > 0:
        with open(run_data.not_found_path, "w") as not_found_file:
            for acc in not_found:
                not_found_file.write(acc + "\n")

    run_data.num_found = len(found)
    run_data.num_not_found = len(not_found)

    return run_data


def build_base_link(dl_acc, assembly_name):
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    prefix, number = dl_acc.split("_")

    p1 = number[0:3]
    p2 = number[3:6]
    p3 = number[6:9]

    path = f"{prefix}/{p1}/{p2}/{p3}/{dl_acc}_{assembly_name}/"

    return base_url + path
