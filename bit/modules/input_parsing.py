from pathlib import Path
from bit.modules.general import (report_message,
                                 notify_premature_exit)


accepted_R1_designations = ["_R1_", "_R1.", "-R1.", "-R1-", ".R1.", "_1."]
accepted_R2_designations = ["_R2_", "_R2.", "-R2.", "-R2-", ".R2.", "_2."]
accepted_read_extensions  = [".fq.gz", ".fastq.gz"]


def get_input_reads_dict_from_dir(reads_dir, SE = False):
    """
    This scans the input directory for read files and returns a dictionary of prefix (hopefully sample names), and read paths
    Place-holder for if i want to implement single-end in the future
    """

    p = Path(reads_dir)
    fastqs = sorted(
        [f for f in p.iterdir() if f.is_file() and any(str(f).endswith(ext) for ext in accepted_read_extensions)],
        key=lambda x: x.name
    )
    reads_dict = {}

    for fq in fastqs:
        name = fq.name
        samp, which = parse_read_filename(name)
        if which:
            reads_dict.setdefault(samp, {})[which] = str(fq.resolve())

    bad = []
    for s, pair in reads_dict.items():
        if "R1" not in pair or "R2" not in pair:
            missing = "R1" if "R1" not in pair else "R2"
            bad.append(f"    sample {s!r} didn't have an {missing} detected")
    if bad:
        print("")
        print("\n".join(bad))
        message = f"Check the input reads dir ('{reads_dir}') for the above issue(s)."
        report_message(message, initial_indent="    ", subsequent_indent="    ")
        notify_premature_exit()

    return reads_dict


def parse_read_filename(file_path):
    """
    Parses a read filename (not full path) to extract the sample name and read designation (R1 or R2).
    Returns a tuple of (sample_name, read_designation).
    """
    filename = Path(file_path).name
    for tag in accepted_R1_designations:
        if tag in filename:
            return filename.split(tag)[0], "R1"
    for tag in accepted_R2_designations:
        if tag in filename:
            return filename.split(tag)[0], "R2"
    # If no valid designation found
    return None, None


def get_input_reads_dict_from_paths(r1_path, r2_path = None):
    """
    Takes explicint input paths for reads and returns a dictionary of prefix (hopefully sample names), and read paths
    """
    reads_dict = {}

    for designation, path in (("R1", r1_path), ("R2", r2_path)):
        if path is None:
            continue
        p = Path(path)
        validate_extension(p)
        samp, which = parse_read_filename(p.name)
        if which is None:
            message = f"We couldn't detect an R1/R2 designation in {p.name!r}."
            report_message(message, initial_indent="    ", subsequent_indent="    ")
            joined = "\n      ".join(accepted_R1_designations)
            print(f"\n    Supported designations are:\n      {joined}")
            print("\n     Or their R2 equivalents.")
            notify_premature_exit()

        if which != designation:
            # user passed as -1 but filename says R2, for example
            message = f"Expected designation {designation!r} but file {p.name!r} parsed as {which}."
            report_message(message, initial_indent="    ", subsequent_indent="    ")
            notify_premature_exit()
        reads_dict.setdefault(samp, {})[which] = str(p.resolve())

    return reads_dict


def validate_extension(file_path):
    if not any(str(file_path).endswith(ext) for ext in accepted_read_extensions):
        message = f"\n{file_path.name!r} doesn't have a supported extension\n"
        report_message(message, initial_indent="   ", subsequent_indent="   ")
        joined = "\n      ".join(accepted_read_extensions)
        print(f"\n    Supported extensions are:\n      {joined}")
        notify_premature_exit()
