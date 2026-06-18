import sys
import os
import shutil
import subprocess
import pandas as pd
from io import StringIO
from bit.modules.general import color_text, wprint, report_message


# the dataformat field names we pull and the friendlier column names we rename them to;
# order here is also the output column order.
# field names verified against `dataformat tsv genome --help` for datasets v18.30.1;
# they have shifted across major versions (e.g. `source-database` existed in older docs
# but is not a valid field in v18), so re-check this list if bumping the pinned version.
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

    check_datasets_available()

    # mutually-exclusive / co-required arg checks (mirrors the GTDB tool's up-front validation)
    if args.get_taxon_counts and not args.target_taxon:
        report_message("A specific taxon needs to also be provided to the `-t` flag in "
                       "order to use `--get-taxon-counts`.", "yellow")
        print("")
        wprint("  E.g.: bit get-accs-from-ncbi --get-taxon-counts -t Alteromonas")
        print("")
        sys.exit(0)

    if not args.target_taxon:
        report_message("A target taxon (name or NCBI taxID) needs to be provided to the "
                       "`-t` flag.", "yellow")
        print("")
        wprint("  E.g.: bit get-accs-from-ncbi -t Alteromonas")
        print("")
        sys.exit(0)

    if args.get_taxon_counts:
        get_taxon_count(args)
        sys.exit(0)

    get_accessions(args)
    sys.exit(0)


def check_datasets_available():
    """ ensures the NCBI `datasets` and `dataformat` CLIs are on PATH

    these are hard conda deps (ncbi-datasets-cli) in meta.yaml, but unlike a python
    import they won't fail loudly on their own when called via subprocess, so we check
    here to give a clean message rather than a cryptic FileNotFoundError
    """

    missing = [tool for tool in ("datasets", "dataformat") if shutil.which(tool) is None]

    if missing:
        report_message("This program requires the NCBI Datasets command-line tools, but "
                       f"the following were not found on your PATH: {', '.join(missing)}.",
                       "red")
        print("")
        wprint("These are normally installed as a dependency (conda package "
               "`ncbi-datasets-cli`). You can install them with:")
        print("")
        print(color_text("    conda install -c conda-forge ncbi-datasets-cli\n"))
        sys.exit(1)


def build_summary_command(args):
    """ assembles the `datasets summary genome taxon ...` argv list

    returned as a list (not a shell string) so it's passed to subprocess without a shell,
    which avoids any quoting issues with taxon names that contain spaces
    """

    cmd = ["datasets", "summary", "genome", "taxon", str(args.target_taxon),
           "--as-json-lines"]

    # --source: genbank vs refseq vs both
    if args.source == "refseq":
        cmd += ["--assembly-source", "RefSeq"]
    elif args.source == "genbank":
        cmd += ["--assembly-source", "GenBank"]

    # only NCBI RefSeq "reference" genomes (the broad-diversity subset);
    # corollary to the GTDB tool's --refseq-reference-genomes-only
    if args.reference_only:
        cmd += ["--reference"]

    # assembly level(s), comma-separated, e.g. "complete,chromosome"
    if args.assembly_level:
        cmd += ["--assembly-level", args.assembly_level]

    # only annotated genomes
    if args.annotated_only:
        cmd += ["--annotated"]

    return cmd


def run_ncbi_summary(args):
    """ runs `datasets summary ... | dataformat tsv genome ...` and returns a parsed DataFrame

    we run the two processes explicitly and pipe between them rather than using shell=True,
    so we keep control of argument quoting and can surface stderr from each stage clearly
    """

    summary_cmd = build_summary_command(args)
    dataformat_cmd = ["dataformat", "tsv", "genome", "--fields", ",".join(FIELD_MAP.keys())]

    report_message(f"Querying NCBI for genomes under taxon '{args.target_taxon}'...", "yellow")
    print("")

    summary_proc = subprocess.Popen(summary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dataformat_proc = subprocess.Popen(dataformat_cmd, stdin=summary_proc.stdout,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # allow summary_proc to receive SIGPIPE if dataformat exits early
    summary_proc.stdout.close()

    # collect dataformat's stdout+stderr first (it's the downstream process), then drain
    # summary's stderr. we read summary_proc.stderr directly rather than calling
    # communicate() a second time, since stdout was already closed above and a second
    # communicate() after the pipe is torn down doesn't reliably return stderr.
    out, dataformat_err = dataformat_proc.communicate()
    summary_err = summary_proc.stderr.read()
    summary_proc.stderr.close()
    summary_proc.wait()

    out_text = out.decode("utf-8", errors="replace")
    summary_err_text = summary_err.decode("utf-8", errors="replace") if summary_err else ""
    dataformat_err_text = dataformat_err.decode("utf-8", errors="replace") if dataformat_err else ""

    # a field-name mismatch (or other dataformat problem) shows up as a non-zero exit from
    # dataformat with a message on ITS stderr, not datasets'. surface that explicitly,
    # since it's the most likely failure when the installed CLI version drifts.
    if dataformat_proc.returncode != 0:
        report_message("The NCBI `dataformat` step failed:", "red")
        print("")
        wprint(dataformat_err_text.strip() if dataformat_err_text.strip() else
               "(no error message was returned)")
        print("")
        wprint("This often means a requested field name isn't valid for your installed "
               "NCBI Datasets version. Check the valid names with "
               + color_text("dataformat tsv genome --help", "none") + ".")
        print("")
        sys.exit(1)

    # `datasets` exits non-zero (and prints to stderr) when a taxon yields no genomes;
    # we treat "no records" as an empty result rather than a hard failure
    if summary_proc.returncode != 0 and not out_text.strip():
        # distinguish "no genomes found" from a real error by inspecting stderr
        lowered = summary_err_text.lower()
        if "does not match any" in lowered or \
           "no genomes" in lowered or \
           "no assemblies" in lowered or \
           "no records" in lowered:
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

    # dataformat emits human-readable headers (e.g. "Assembly Accession"); we re-key by
    # position to our snake_case names since we control the --fields order above
    if len(tab.columns) == len(FIELD_MAP):
        tab.columns = list(FIELD_MAP.values())
    else:
        # column count drift means our assumed field set no longer matches the installed
        # dataformat; surface it rather than silently misaligning columns
        report_message("The columns returned by `dataformat` didn't match what was "
                       "expected. The installed NCBI Datasets version may use different "
                       "field names. Got these columns:", "red")
        print("")
        wprint(", ".join(str(c) for c in tab.columns))
        print("")
        sys.exit(1)

    # normalize missing values to a consistent "NA" in the output table. this covers:
    #   - true NaN (empty fields under dtype=str)
    #   - empty / whitespace-only strings from dataformat
    #   - NCBI's own literal "na" placeholders in any casing (na, NA, Na, ...)
    # every column is read as string (dtype=str) or NaN, so we map all of them rather than
    # guarding on dtype == "object", which doesn't reliably hold across pandas versions.
    def normalize_cell(value):
        if pd.isna(value):
            return "NA"
        stripped = str(value).strip()
        if stripped == "" or stripped.lower() == "na":
            return "NA"
        return stripped

    tab = tab.apply(lambda col: col.map(normalize_cell))

    return tab


def get_accessions(args):
    """ pulls accessions + metadata for the target taxon and writes the two output files """

    tab = run_ncbi_summary(args)

    if tab.empty:
        report_message(f"No genomes were found under taxon '{args.target_taxon}' with the "
                       "specified filters.", "yellow")
        print("")
        sys.exit(0)

    # build output filenames from the (filesystem-safe) taxon
    taxon_for_filename = str(args.target_taxon).replace(" ", "-").replace("/", "-")

    suffix_bits = []
    if args.reference_only:
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
    wprint("Associated taxonomy and assembly info of these targets written to:")
    wprint("  " + color_text(tab_out_filename))
    print("")


def get_taxon_count(args):
    """ reports how many genomes match the target taxon under the current filters

    implemented via the same summary query (counting returned rows) rather than a
    separate count endpoint, so the count always reflects the same filters that
    `get_accessions` would apply
    """

    tab = run_ncbi_summary(args)

    count = len(tab.index)

    print("")
    wprint(f"  There are {count:,} genome(s) under taxon '{args.target_taxon}' with the "
           "specified filters.")
    print("")