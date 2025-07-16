import os
import re
import sys
import textwrap
from importlib.resources import files
from pathlib import Path
from importlib.metadata import version


def get_package_path(rel_path = ""):
    return os.path.join(files("bit") / rel_path)


def color_text(text, color = 'green'):
    tty_colors = {
        'green' : '\033[0;32m%s\033[0m',
        'yellow' : '\033[0;33m%s\033[0m',
        'red' : '\033[0;31m%s\033[0m'
    }
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text, width = 80, initial_indent = "  ",
           subsequent_indent = "  "):
    print(textwrap.fill(text, width = width, initial_indent = initial_indent,
          subsequent_indent = subsequent_indent, break_on_hyphens = False))


def report_message(message, color = "yellow", width = 80, initial_indent = "  ",
                   subsequent_indent = "  "):
    print("")
    wprint(color_text(message, color), width = width,
           initial_indent = initial_indent,
           subsequent_indent = subsequent_indent)


def report_failure(message, color = "red"):
    print("")
    wprint(color_text(message, color))
    print("\nExiting for now :(\n")
    sys.exit(1)


def get_input_reads_dict(reads_dir, SE = False):
    """
    This scans the input directory for read files and returns a dictionary of prefix (hopefully sample names), and read paths
    Place-holder for if i want to implement single-end in the future
    """

    accepted_R1_designations = ["_R1_", "_R1.", "-R1.", "-R1-", ".R1.", "_1."]
    accepted_R2_designations = ["_R2_", "_R2.", "-R2.", "-R2-", ".R2.", "_2."]
    accepted_read_extensions  = [".fq.gz", ".fastq.gz"]

    p = Path(reads_dir)
    fastqs = sorted(
        [f for f in p.iterdir() if f.is_file() and any(str(f).endswith(ext) for ext in accepted_read_extensions)],
        key=lambda x: x.name
    )
    reads_dict = {}

    for fq in fastqs:
        name = fq.name
        # find which read (R1 or R2) it is
        which = None
        for tag in accepted_R1_designations:
            if tag in name:
                which = "R1"
                samp = name.split(tag)[0]
                break
        for tag in accepted_R2_designations:
            if tag in name:
                which = "R2"
                samp = name.split(tag)[0]
                break
        if which is None:
            continue

        # record it
        reads_dict.setdefault(samp, {})[which] = str(fq.resolve())

    bad = []
    for s, pair in reads_dict.items():
        if "R1" not in pair or "R2" not in pair:
            missing = "R1" if "R1" not in pair else "R2"
            bad.append(f"  ❌ sample {s!r} didn't have an {missing} detected")
    if bad:
        print("")
        print("\n".join(bad))
        reads_dir = reads_dir if reads_dir.endswith("/") else reads_dir + "/"
        report_failure(f"Check the input reads dir ('{reads_dir}') for the above issues.")

    return reads_dict


def log_command_run(full_cmd_executed, log_dir, log_file = None):
    if log_file is None:
        log_file = Path(log_dir) / "command-execution-info.txt"
    else:
        log_file = Path(log_file)
    with log_file.open("w") as f:
        f.write(f"bit version: v{version('bit')}\n\n")
        f.write(f"Rendered command:\n{full_cmd_executed}\n\n")


def check_files_are_found(paths_list):
    for path in paths_list:
        if not Path(path).is_file():
            print(f"\n    We were not able to find the input file: {path}")
            notify_premature_exit()


def notify_premature_exit():
    print("\n    Exiting for now :(\n")
    sys.exit(1)


def tee(msg, log_path, end="\n"):
    ANSI_ESCAPE = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    print(msg, end=end, file=sys.stdout)
    clean = ANSI_ESCAPE.sub('', msg)
    with open(log_path, "a") as f:
        f.write(clean + end)
