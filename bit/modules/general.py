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
