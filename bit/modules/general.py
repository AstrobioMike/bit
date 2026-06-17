import os
import re
import shutil
import sys
import textwrap
import csv
import subprocess
from importlib.resources import files
from contextlib import contextmanager
import itertools
import threading
import time
from pathlib import Path
from importlib.metadata import version
import urllib.request
from tqdm import tqdm # type: ignore


def get_package_path(rel_path = ""):
    return os.path.join(files("bit") / rel_path)


def color_text(text, color = 'green'):
    tty_colors = {
        'green' : '\033[0;32m%s\033[0m',
        'orange': '\033[38;5;208m%s\033[0m',
        'red' : '\033[0;31m%s\033[0m',
        'teal' : '\033[0;36m%s\033[0m',
        'yellow' : '\033[0;33m%s\033[0m'
    }

    if sys.stdout.isatty() and color != "none":
        return tty_colors[color] % text
    else:
        return text


def wprint(text, width = 80, initial_indent = "  ",
           subsequent_indent = "  "):
    print(textwrap.fill(text, width = width, initial_indent = initial_indent,
          subsequent_indent = subsequent_indent, break_on_hyphens = False))


def report_message(message, color = "yellow", width = 80,
                   initial_indent = "  ", subsequent_indent = "  ",
                   join = True, leading_newline = True, trailing_newline = False):

    if leading_newline:
        print("")

    if join:
        message = " ".join(message.split())

    wrapped = textwrap.fill(
        message,
        width = width,
        initial_indent = initial_indent,
        subsequent_indent = subsequent_indent,
        break_on_hyphens = False,
    )

    print(color_text(wrapped, color))

    if trailing_newline:
        print("")


def report_failure(message="", color = "red"):
    if message:
        print("")
        wprint(color_text(message, color))
    else:
        print("")
    print("\nExiting for now :(\n")
    sys.exit(1)


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


def check_if_output_dir_exists(output_dir, force_overwrite=False):
    if Path(output_dir).is_dir():
        if not force_overwrite:
            print(f"\n    {color_text(f'Output directory already exists: {output_dir}', 'yellow')}")
            print("\n    Please specify a different output directory or add the `-F/--force-overwrite` flag.")
            notify_premature_exit()
        else:
            shutil.rmtree(output_dir)


def notify_premature_exit(exit_code=1):
    print("\n  Exiting for now :(\n", file=sys.stderr)
    sys.exit(exit_code)


def tee(msg, log_path, end="\n"):
    ANSI_ESCAPE = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    print(msg, end=end, file=sys.stdout)
    clean = ANSI_ESCAPE.sub('', msg)
    with open(log_path, "a") as f:
        f.write(clean + end)


def is_gzipped(file_path):
    with open(file_path, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def download_with_tqdm(url, target, filename=None, urlopen=False):
    opener = urllib.request.build_opener()
    opener.addheaders = [('User-Agent', 'curl/8.0')]
    urllib.request.install_opener(opener)
    with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=target, ncols = 90) as t:
        def reporthook(block_num, block_size, total_size):
            if total_size > 0:
                t.total = total_size
            t.update(block_size)
        if not urlopen:
            downloaded_path, _ = urllib.request.urlretrieve(url, filename, reporthook=reporthook)
            sys.stdout.write("")
            return downloaded_path
        else:
            dl = urllib.request.urlopen(url, reporthook=reporthook)
            sys.stdout.write("")
            return dl


def attempt_to_make_dir(dir_path):
    try:
        os.makedirs(dir_path, exist_ok=True)
    except Exception as e:
        report_failure(f"Failed to create directory '{dir_path}' with the following error:\n{e}")


def sniff_delimiter(line):

    try:
        dialect = csv.Sniffer().sniff(line, delimiters='\t,|; ')
        return dialect.delimiter
    except csv.Error:
        return None


def colnames(args):

    input_source = args.input_file
    if hasattr(input_source, "readline"):
        header_line = input_source.readline().rstrip('\n\r')
    else:
        with open(input_source) as f:
            header_line = f.readline().rstrip('\n\r')

    if not header_line:
        print("  File appears to be empty.", file=sys.stderr)
        sys.exit(1)

    delimiter = sniff_delimiter(header_line)
    if delimiter is None:
        print(f"\n  {color_text('Could not detect a delimiter :(', 'yellow')}")
        print("\n  But here's the first line:")
        print(f"    {header_line}")
        print()
        sys.exit(1)

    columns = header_line.split(delimiter)

    width = len(str(len(columns)))

    for i, name in enumerate(columns, 1):
        print(f"  {i:>{width}}  {name}")


def report_version():

    from datetime import datetime
    from importlib.metadata import version

    ver = f"v{version('bit')}"
    print()
    print(f"{' ' * 33} bit {color_text(ver, 'green')}")
    print(f"{' ' * 25} github.com/AstrobioMike/bit\n")

    print("    If you happen to find this toolset useful in your work, please be sure to")
    print("    cite it :)\n")

    print("  Lee M. bit: a multipurpose collection of bioinformatics tools. F1000Research 2022, 11:122")
    print("  https://doi.org/10.12688/f1000research.79530.1\n")

    today = datetime.today().strftime('%A')
    signoff = f"Happy {today} :)"
    print(f"{'':>51}{color_text(signoff,'green')}\n")


@contextmanager
def spinner(in_progress_msg, complete_msg):
    """Show a spinner while a block runs; report elapsed time only if >= 60 s."""
    done = threading.Event()
    elapsed = [0.0]

    def spin():
        for char in itertools.cycle("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"):
            if done.is_set():
                break
            sys.stderr.write(f"\r    {char} {in_progress_msg}")
            sys.stderr.flush()
            time.sleep(0.1)
        if elapsed[0] >= 60:
            mins, secs = divmod(int(elapsed[0]), 60)
            time_str = f"(took ~{mins} min and {secs} sec)"
        else:
            time_str = ""
        sys.stderr.write(f"\r    ✔ {complete_msg}{time_str}          \n")
        sys.stderr.flush()
        time.sleep(0.1)

    t = threading.Thread(target=spin)
    t.start()
    start_time = time.monotonic()
    try:
        yield
    finally:
        elapsed[0] = time.monotonic() - start_time
        done.set()
        t.join()


def check_bam_file_is_indexed(bam_file):

    header_result = subprocess.run(
        ["samtools", "view", "-H", bam_file],
        capture_output=True, text=True
    )
    is_sorted = any(
        "SO:coordinate" in line
        for line in header_result.stdout.splitlines()
        if line.startswith("@HD")
    )

    sorted_for_you = False
    if not is_sorted:

        message = """
                    We're sorting and indexing the BAM for you. Why? Because it's common courtesy.
                    All the programs that don't do this for us when it's needed are just big jerks!
                    """

        report_message(message, color="orange", initial_indent="    ", subsequent_indent="    ")

        stem = bam_file[:-4] if bam_file.endswith(".bam") else bam_file
        sorted_bam = stem + ".sorted.bam"
        report_message("Sorting BAM file...")
        with spinner("", ""):
            subprocess.run(["samtools", "sort", "-@", "4", "-o", sorted_bam, bam_file], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        bam_file = sorted_bam
        sorted_for_you = True

    if not Path(bam_file + ".bai").is_file():
        if not sorted_for_you:
            message = """
                      We're indexing the BAM for you. Why? Because it's common courtesy.
                      All the programs that don't do this for us when it's needed are just big jerks!
                      """

            report_message(message, color="orange", initial_indent="    ", subsequent_indent="    ")

        report_message("Indexing BAM file...")
        with spinner("", ""):
            subprocess.run(["samtools", "index", "-@", "4", bam_file], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return bam_file
