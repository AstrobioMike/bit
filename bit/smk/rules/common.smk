# helpers shared by the bit workflows (pulled in with `include:`)

import os
import sys

_tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}

def env(name):
    """Path of one of the shared conda environment files (SMK is set by the workflow's Snakefile)."""
    return os.path.join(SMK, "envs", name)


def color_text(text, color = 'green'):
    if sys.stdout.isatty():
        return _tty_colors[color] % text
    return text


def exit_with_message(message):
    print(color_text(f"\n    {message}\n", "yellow"))
    print("    Exiting for now.\n")
    sys.exit(1)


def read_id_list(path):
    """Reads a single-column file of IDs (blank lines ignored); they must be unique."""
    if not os.path.isfile(path):
        exit_with_message(f"The file '{path}' was not found :(")
    with open(path) as f:
        ids = [line.strip() for line in f if line.strip()]
    if not ids:
        exit_with_message(f"The file '{path}' holds no IDs :(")
    if len(set(ids)) != len(ids):
        exit_with_message(f"Not all IDs in the '{path}' file are unique :(")
    return ids


def combine_benchmarks(benchmarks_dir = "benchmarks", out_name = "ALL-benchmarks.tsv"):
    """Combines the per-job snakemake benchmark files into one table (first column: the process)."""
    if not os.path.isdir(benchmarks_dir):
        return
    files = sorted(f for f in os.listdir(benchmarks_dir) if f.endswith("-benchmarks.tsv") and f != out_name)
    if not files:
        return
    lines, header = [], None
    for name in files:
        with open(os.path.join(benchmarks_dir, name)) as f:
            rows = [l.rstrip("\n") for l in f if l.strip()]
        if not rows:
            continue
        header = header or rows[0]
        process = name[:-len("-benchmarks.tsv")]
        lines += [f"{process}\t{row}" for row in rows[1:]]
    if header is None:
        return
    with open(os.path.join(benchmarks_dir, out_name), "w") as out:
        out.write(f"process\t{header}\n")
        out.write("".join(f"{l}\n" for l in lines))
