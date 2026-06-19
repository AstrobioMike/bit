import sys
import subprocess
import pandas as pd # type: ignore
from io import StringIO
from bit.modules.general import color_text, wprint, report_message


FIELD_MAP = {
    "accession": "accession",
    "organism-name": "organism_name",
    "organism-tax-id": "tax_id",
    "assminfo-name": "assembly_name",
    "assminfo-level": "assembly_level",
    "assminfo-refseq-category": "refseq_category",
    "assmstats-total-sequence-len": "total_sequence_length",
    "checkm-completeness": "checkm_completeness",
    "checkm-contamination": "checkm_contamination",
    "source_database": "source_database",
}

def get_accessions_from_ncbi(args):

    if args.get_taxon_counts:
        get_taxon_count(args)
        sys.exit(0)

    get_accessions(args)
    sys.exit(0)


def build_summary_command(args):

    cmd = ["datasets", "summary", "genome", "taxon", str(args.target_taxon),
           "--as-json-lines"]

    if args.source == "refseq":
        cmd += ["--assembly-source", "RefSeq"]
    elif args.source == "genbank":
        cmd += ["--assembly-source", "GenBank"]

    if args.reference_genomes_only:
        cmd += ["--reference"]

    if args.assembly_level:
        cmd += ["--assembly-level", ",".join(args.assembly_level)]

    if args.annotated_only:
        cmd += ["--annotated"]

    return cmd


def run_ncbi_summary(args):
    """ 
    runs `datasets summary ... | dataformat tsv genome ...` and returns a parsed DataFrame
    """

    summary_cmd = build_summary_command(args)
    dataformat_cmd = ["dataformat", "tsv", "genome", "--fields", ",".join(FIELD_MAP.keys())]

    report_message(f"Querying NCBI for genomes under target '{args.target_taxon}'...", "yellow")

    summary_proc = subprocess.Popen(summary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dataformat_proc = subprocess.Popen(dataformat_cmd, stdin=summary_proc.stdout,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    summary_proc.stdout.close()

    out, dataformat_err = dataformat_proc.communicate()
    summary_err = summary_proc.stderr.read()
    summary_proc.stderr.close()
    summary_proc.wait()

    out_text = out.decode("utf-8", errors="replace")
    summary_err_text = summary_err.decode("utf-8", errors="replace") if summary_err else ""
    dataformat_err_text = dataformat_err.decode("utf-8", errors="replace") if dataformat_err else ""

    # a field-name mismatch (or other dataformat problem) shows up as a non-zero exit from
    # dataformat with a message on its stderr. printing that error out if it happens
    if dataformat_proc.returncode != 0:
        report_message("The NCBI `dataformat` step failed:", "red")
        print("")
        wprint(dataformat_err_text.strip() if dataformat_err_text.strip() else
               "(no error message was returned)")
        print("")
        sys.exit(1)

    if summary_proc.returncode != 0 and not out_text.strip():

        lowered = summary_err_text.lower()
        if "is not recognized" in lowered or "is not exact" in lowered:
            return pd.DataFrame(columns=list(FIELD_MAP.values()))

        report_message("The NCBI `datasets` query failed:", "red")
        print("")
        wprint(summary_err_text.strip() if summary_err_text.strip() else
               "(no error message was returned)")
        print("")
        sys.exit(1)

    if not out_text.strip():
        return pd.DataFrame(columns=list(FIELD_MAP.values()))

    tab = pd.read_csv(StringIO(out_text), sep="\t", dtype=str, low_memory=False)

    if len(tab.columns) == len(FIELD_MAP):
        tab.columns = list(FIELD_MAP.values())
    else:
        report_message("The columns returned by `dataformat` didn't match what was "
                       "expected. The installed NCBI Datasets version may use different "
                       "field names. Got these columns:", "red")
        print("")
        wprint(", ".join(str(c) for c in tab.columns))
        print("")
        sys.exit(1)

    def normalize_cell(value):
        if pd.isna(value):
            return "NA"
        stripped = str(value).strip()
        if stripped == "" or stripped.lower() == "na":
            return "NA"
        return stripped

    tab = tab.apply(lambda col: col.map(normalize_cell))

    if "source_database" in tab.columns:
        source_db_map = {
            "SOURCE_DATABASE_GENBANK": "genbank",
            "SOURCE_DATABASE_REFSEQ": "refseq",
        }
        tab["source_database"] = tab["source_database"].map(
            lambda value: source_db_map.get(value, value)
        )

    return tab


def get_accessions(args):
    """ 
    pulls accessions + metadata for the target taxon and writes the two output files
    """

    tab = run_ncbi_summary(args)

    if tab.empty:
        report_message(f"No genomes were found under taxon '{args.target_taxon}' with the "
                       "specified filters.", "none", initial_indent="    ", subsequent_indent="    ")
        print("")
        sys.exit(0)

    # build output filenames from the (filesystem-safe) taxon
    taxon_for_filename = str(args.target_taxon).replace(" ", "-").replace("/", "-")

    suffix_bits = []
    if args.reference_genomes_only:
        suffix_bits.append("refseq-ref")
    elif args.source != "both":
        suffix_bits.append(args.source.lower())
    suffix = ("-" + "-".join(suffix_bits)) if suffix_bits else ""

    acc_out_filename = f"NCBI-{taxon_for_filename}{suffix}-accessions.txt"
    tab_out_filename = f"NCBI-{taxon_for_filename}{suffix}-metadata.tsv"

    tab.to_csv(tab_out_filename, sep="\t", index=False)

    target_accs = tab["accession"].tolist()
    with open(acc_out_filename, "w") as out:
        for acc in target_accs:
            out.write(acc + "\n")

    print("")
    wprint(f"Wrote {len(target_accs):,} accession(s) to:")
    wprint("  " + color_text(acc_out_filename))
    print("")
    wprint("Associated metadata of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def get_taxon_count(args):

    tab = run_ncbi_summary(args)

    count = len(tab.index)

    report_message(f"There are {count:,} genome(s) under target '{args.target_taxon}' with the "
                    "specified filters.", "none", initial_indent="    ", subsequent_indent="    ",
                    trailing_newline=True)
