import re
import itertools

def sanitize_assembly_name(name):
    sanitized = re.sub(r"[\s/,#()\[\]]", "_", name)
    sanitized = re.sub(r"_+", "_", sanitized)
    return sanitized.strip("_")


def build_base_link(dl_acc, assembly_name):
    """
    fallback directory URL builder, used only when the summary has no ftp_path
    """
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    prefix, rest = dl_acc.split("_", 1)
    number = rest.split(".")[0]
    p1, p2, p3 = number[0:3], number[3:6], number[6:9]
    dir_basename = f"{dl_acc}_{sanitize_assembly_name(assembly_name)}"
    path = f"{prefix}/{p1}/{p2}/{p3}/{dir_basename}/"
    return base_url + path, dir_basename


def parse_ncbi_assembly_summary(assembly_summary_file, run_data):

    wanted_dict = {}

    for acc in run_data.wanted_accs:
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    FORMAT_EXTENSIONS = {
        "genbank":     ("_genomic.gbff.gz",            ".gb.gz"),
        "fasta":       ("_genomic.fna.gz",             ".fasta.gz"),
        "protein":     ("_protein.faa.gz",             ".faa.gz"),
        "gff":         ("_genomic.gff.gz",             ".gff.gz"),
        "nt_cds":      ("_cds_from_genomic.fna.gz",    "_cds_from_genomic.fna.gz"),
        "feature_tab": ("_feature_table.txt.gz",       ".tsv.gz"),
        "report":      ("_assembly_report.txt",         "_assembly_report.txt"),
        "stats":       ("_assembly_stats.txt",          "_assembly_stats.txt"),
    }

    # column positions are resolved from the file's header row (the slim table
    # written by slim_ncbi_assembly_summary carries a clean header). If the file
    # has no header (a legacy positional NCBI summary, or an old cached copy),
    # fall back to the fixed NCBI assembly_summary positions so existing data
    # still parses. NAMES maps the logical fields this reader needs to their
    # NCBI column names.
    NAMES = {
        "accession": "assembly_accession",
        "assembly_name": "asm_name",
        "taxid": "taxid",
        "org_name": "organism_name",
        "infra_name": "infraspecific_name",
        "version_status": "version_status",
        "assembly_level": "assembly_level",
        "ftp_path": "ftp_path",
    }
    LEGACY_POS = {
        "accession": 0, "assembly_name": 15, "taxid": 5, "org_name": 7,
        "infra_name": 8, "version_status": 10, "assembly_level": 11,
        "ftp_path": 19,
    }

    def _resolve_positions(header_line):
        """header line (clean or '#'-prefixed) -> {logical: index}."""
        names = header_line.lstrip("#").rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(names)}
        return {logical: idx[col] for logical, col in NAMES.items()}

    with open(run_data.ncbi_sub_table_path, "w") as out_file:

        cols = ["input_accession", "found_accession", "assembly_name", "taxid", "organism_name", "infraspecific_name", "version_status", "assembly_level", "http_base_link"]
        if run_data.wanted_format:
            cols.extend(["target_link", "local_destination"])
        out_file.write("\t".join(cols) + "\n")

        with open(assembly_summary_file, "r") as assemblies:
            first = assemblies.readline()
            first_stripped = first.lstrip("#").rstrip("\n")
            if first_stripped.split("\t", 1)[0] == "assembly_accession":
                pos = _resolve_positions(first)          # header present
                remaining = assemblies                   # data starts next line
            else:
                pos = dict(LEGACY_POS)                    # headerless -> legacy
                # the line we already read is data; chain it back on
                remaining = itertools.chain([first], assemblies)

            p_acc = pos["accession"]; p_name = pos["assembly_name"]
            p_tax = pos["taxid"]; p_org = pos["org_name"]; p_infra = pos["infra_name"]
            p_ver = pos["version_status"]; p_lvl = pos["assembly_level"]
            p_ftp = pos["ftp_path"]
            max_pos = max(pos.values())

            for line in remaining:
                fields = line.rstrip("\n").split("\t")
                if not fields or len(fields) <= max_pos:
                    continue
                root = fields[p_acc].split(".")[0]
                if root in wanted_dict:
                    found.add(wanted_dict[root])

                    dl_acc = fields[p_acc].strip() if fields[p_acc].strip() else "NA"
                    assembly_name = fields[p_name].strip() if fields[p_name].strip() else "NA"
                    taxid = fields[p_tax].strip() if fields[p_tax].strip() else "NA"
                    org_name = fields[p_org].strip() if fields[p_org].strip() else "NA"
                    infra_name = fields[p_infra].strip() if fields[p_infra].strip() else "NA"
                    version_status = fields[p_ver].strip() if fields[p_ver].strip() else "NA"
                    assembly_level = fields[p_lvl].strip() if fields[p_lvl].strip() else "NA"
                    ftp_path = fields[p_ftp].strip()
                    if ftp_path and ftp_path.lower() != "na":
                        http_path = ftp_path.replace("ftp://", "https://").rstrip("/") + "/"
                        dir_basename = http_path.rstrip("/").split("/")[-1]
                    elif assembly_name != "NA" and dl_acc != "NA":
                        http_path, dir_basename = build_base_link(dl_acc, assembly_name)
                    else:
                        http_path, dir_basename = "NA", "NA"

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
                        if http_path != "NA":
                            target_link = f"{http_path}{dir_basename}{ncbi_ext}"
                        else:
                            target_link = "NA"
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
