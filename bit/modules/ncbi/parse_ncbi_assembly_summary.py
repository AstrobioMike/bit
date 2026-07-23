import re
import pyarrow as pa # type: ignore
import pyarrow.compute as pc # type: ignore
import pyarrow.dataset as ds # type: ignore


def sanitize_assembly_name(name):
    sanitized = re.sub(r"[\s/,#()\[\]]", "_", name)
    sanitized = re.sub(r"_+", "_", sanitized)
    return sanitized.strip("_")


def build_base_link(dl_acc, assembly_name):
    # only used if the ftp_path is missing or NA; otherwise we use that directly
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    prefix, rest = dl_acc.split("_", 1)
    number = rest.split(".")[0]
    p1, p2, p3 = number[0:3], number[3:6], number[6:9]
    dir_basename = f"{dl_acc}_{sanitize_assembly_name(assembly_name)}"
    path = f"{prefix}/{p1}/{p2}/{p3}/{dir_basename}/"
    return base_url + path, dir_basename


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


_NEEDED_COLUMNS = [
    "assembly_accession", "asm_name", "taxid", "organism_name",
    "infraspecific_name", "version_status", "assembly_level", "ftp_path",
]


def _clean(value):
    if value is None:
        return "NA"
    v = str(value).strip()
    return v if v else "NA"


def _resolve_links(dl_acc, assembly_name, ftp_path):
    """
      - prefer ftp_path (ftp:// -> https://, single trailing slash)
      - else rebuild from accession + assembly name (build_base_link)
      - else NA
    """
    if ftp_path and ftp_path.lower() != "na":
        http_path = ftp_path.replace("ftp://", "https://").rstrip("/") + "/"
        dir_basename = http_path.rstrip("/").split("/")[-1]
    elif assembly_name != "NA" and dl_acc != "NA":
        http_path, dir_basename = build_base_link(dl_acc, assembly_name)
    else:
        http_path, dir_basename = "NA", "NA"
    return http_path, dir_basename


def parse_ncbi_assembly_summary(assembly_summary_file, run_data):
    """
    Look up run_data.wanted_accs in the hosted NCBI Parquet table and write the
    per-accession info table (run_data.ncbi_sub_table_path), plus the not-found list
    """
    # root (version-stripped) accession -> the exact string the user asked for
    wanted_dict = {}
    for acc in run_data.wanted_accs:
        root_acc = acc.strip().split(".")[0]
        wanted_dict[root_acc] = acc.strip()

    found = set()

    dataset = ds.dataset(str(assembly_summary_file), format="parquet")

    roots = list(wanted_dict)
    root_field = pc.replace_substring_regex(ds.field("assembly_accession"),
                                            r"\..*$", "")
    predicate = pc.is_in(root_field,
                         value_set=pa.array(sorted(set(roots)), type=pa.string()))

    with open(run_data.ncbi_sub_table_path, "w") as out_file:

        cols = ["input_accession", "found_accession", "assembly_name", "taxid",
                "organism_name", "infraspecific_name", "version_status",
                "assembly_level", "http_base_link"]
        if run_data.wanted_format:
            cols.extend(["target_link", "local_destination"])
        out_file.write("\t".join(cols) + "\n")

        scanner = dataset.scanner(columns=_NEEDED_COLUMNS, filter=predicate)

        for batch in scanner.to_batches():
            rows = batch.to_pylist()
            for row in rows:
                dl_acc_raw = _clean(row.get("assembly_accession"))
                root = dl_acc_raw.split(".")[0]
                if root not in wanted_dict:
                    continue

                found.add(wanted_dict[root])

                dl_acc = dl_acc_raw
                assembly_name = _clean(row.get("asm_name"))
                taxid = _clean(row.get("taxid"))
                org_name = _clean(row.get("organism_name"))
                infra_name = _clean(row.get("infraspecific_name"))
                version_status = _clean(row.get("version_status"))
                assembly_level = _clean(row.get("assembly_level"))
                ftp_path = (row.get("ftp_path") or "").strip()

                http_path, dir_basename = _resolve_links(dl_acc, assembly_name, ftp_path)

                out_fields = [
                    wanted_dict[root], dl_acc, assembly_name, taxid, org_name,
                    infra_name, version_status, assembly_level, http_path,
                ]
                out_line = "\t".join(out_fields)

                if run_data.wanted_format:
                    ncbi_ext, local_ext = FORMAT_EXTENSIONS[run_data.wanted_format]
                    if http_path != "NA":
                        target_link = f"{http_path}{dir_basename}{ncbi_ext}"
                    else:
                        target_link = "NA"
                    local_path = f"{run_data.output_dir}/{dl_acc}{local_ext}"
                    out_line += "\t" + target_link + "\t" + local_path

                out_file.write(out_line + "\n")

    not_found = set(run_data.wanted_accs) - found

    if len(not_found) > 0:
        with open(run_data.not_found_path, "w") as not_found_file:
            for acc in not_found:
                not_found_file.write(acc + "\n")

    run_data.num_found = len(found)
    run_data.num_not_found = len(not_found)

    return run_data
