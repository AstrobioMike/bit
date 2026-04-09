import os
import re
import sys
import textwrap
import csv
from importlib.resources import files
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


def report_failure(message, color = "red"):
    print("")
    wprint(color_text(message, color))
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


def notify_premature_exit():
    print("\n  Exiting for now :(\n")
    sys.exit(1)


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
    with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=target, ncols = 90) as t:
        def reporthook(block_num, block_size, total_size):
            if total_size > 0:
                t.total = total_size
            t.update(block_size)
        if not urlopen:
            urllib.request.urlretrieve(url, filename, reporthook=reporthook)
            sys.stdout.write("")
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

    with open(args.input_file) as f:
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
        exit()

    columns = header_line.split(delimiter)

    width = len(str(len(columns)))

    for i, name in enumerate(columns, 1):
        print(f"  {i:>{width}}  {name}")


def report_version():

    from datetime import datetime
    from importlib.metadata import version

    ver = f"v{version('bit')}"
    print()
    print(f"{' ' * 22} Bioinformatics Tools (bit) {color_text(ver, 'green')}")
    print(f"{' ' * 25} github.com/AstrobioMike/bit\n")

    print("    If you happen to find this toolset useful in your work, please be sure to")
    print("    cite it :)\n")

    print("  Lee M. bit: a multipurpose collection of bioinformatics tools. F1000Research 2022, 11:122")
    print("  https://doi.org/10.12688/f1000research.79530.1\n")

    today = datetime.today().strftime('%A')
    signoff = f"Happy {today} :)"
    print(f"                                                   {color_text(signoff,'green')}\n")
