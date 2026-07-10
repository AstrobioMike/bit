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
import socket
from pathlib import Path
from importlib.metadata import version
import urllib.request
import urllib.error
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


class _TooSlow(Exception):
    def __init__(self, mbps):
        self.mbps = mbps
        super().__init__(f"too slow: {mbps:.2f} MB/s")


def download_with_tqdm(url, target, filename=None, urlopen=False, leave=True,
                       retries=True, attempts=6, retry_wait=3,
                       speed_gate=False, min_mbps=2.0, probe_seconds=5.0):
    """
    Download `url` to `filename`, showing a tqdm progress bar, with optional
    transient-error retries and optional speed-gated route rerolling

    retries=True (default): transient failures (timeouts, connection resets,
    transient HTTP/URL errors) are retried up to `attempts` times, waiting
    `retry_wait` seconds between tries. A 404 is never retried (raised at once).
    Set retries=False for a single-shot download that raises on any error.

    speed_gate=False (default): don't judge throughput. speed_gate=True: GitHub
    release assets are served from a CDN whose per-connection throughput varies
    by edge node; a slow connection tends to stay slow, and reconnecting often
    lands on a faster edge. With the gate on, non-final attempts measure
    throughput over an initial `probe_seconds` window and, if it's below
    `min_mbps`, abort and reconnect to reroll the route. The final attempt
    accepts whatever speed it gets and runs to completion, so a persistently
    slow network still succeeds.

    Raises the last underlying error if every attempt fails.
    """
    if urlopen:
        opener = urllib.request.build_opener()
        opener.addheaders = [('User-Agent', 'curl/8.0')]
        urllib.request.install_opener(opener)
        return urllib.request.urlopen(url)

    # a single attempt if retries are off
    if not retries:
        attempts = 1

    # resolve total size up front so the bar is bounded/persistent (a plain GET
    # carries Content-Length through GitHub's redirect more reliably than HEAD)
    total = None
    try:
        with urllib.request.urlopen(url) as r:
            cl = r.headers.get("Content-Length")
            total = int(cl) if cl else None
    except Exception:
        total = None

    floor_bytes_per_s = (min_mbps * 1024 * 1024) if speed_gate else 0.0
    last_err = None

    for attempt in range(1, attempts + 1):
        is_final = (attempt == attempts)
        # enforce the speed floor only on non-final attempts (and only if gated)
        floor = 0.0 if is_final else floor_bytes_per_s

        try:
            _stream_once(url, filename, target, total, leave, floor, probe_seconds)
            if leave:
                sys.stderr.write("\n")
            return filename

        except _TooSlow as e:
            # performance failure: reroll immediately, no wait
            last_err = e
            report_message(
                f"that was a slow route (try {attempt}/{attempts - 1}), trying to get a faster one...",
                "yellow", initial_indent="          ", subsequent_indent="          ",
                width=90)
            print()
            continue

        except urllib.error.HTTPError as err:
            # a definitive 404 is never worth retrying
            if err.code == 404:
                raise
            last_err = err
            if is_final:
                raise
            wprint(color_text(
                f"    download failed (attempt {attempt}/{attempts}); retrying...",
                "yellow"))
            time.sleep(retry_wait)
            continue

        except (urllib.error.URLError, socket.timeout, TimeoutError,
                ConnectionError, OSError) as err:
            # transient network failure: wait, then retry
            last_err = err
            if is_final:
                raise
            wprint(color_text(
                f"    download failed (attempt {attempt}/{attempts}); retrying...",
                "yellow"))
            time.sleep(retry_wait)
            continue

    if last_err:
        raise last_err


def _stream_once(url, filename, desc, total, leave, floor_bytes_per_s, probe_seconds):
    """
    Stream url->filename with a tqdm bar. If floor_bytes_per_s > 0, measure
    throughput over the first `probe_seconds` and raise _TooSlow if under floor.
    Network/HTTP errors propagate to the caller's attempt loop.
    """
    chunk = 1024 * 256  # 256 KB reads
    start = time.monotonic()
    probed = False
    req = urllib.request.Request(url, headers={'User-Agent': 'curl/8.0'})
    with urllib.request.urlopen(req) as resp, open(filename, "wb") as out, \
         tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
              desc=desc, ncols=90, leave=leave, total=total) as bar:
        downloaded = 0
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            out.write(buf)
            downloaded += len(buf)
            if total is not None:
                bar.update(min(downloaded, total) - bar.n)
            else:
                bar.update(len(buf))
            if floor_bytes_per_s > 0 and not probed:
                elapsed = time.monotonic() - start
                if elapsed >= probe_seconds:
                    probed = True
                    rate = downloaded / elapsed if elapsed > 0 else 0.0
                    if rate < floor_bytes_per_s:
                        bar.close()
                        raise _TooSlow(rate / (1024 * 1024))
        if total is not None and bar.n < total:
            bar.update(total - bar.n)


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
def spinner(in_progress_msg, complete_msg, indent="    "):
    """Show a spinner while a block runs; report elapsed time only if >= 60 s."""
    done = threading.Event()
    elapsed = [0.0]

    def spin():
        for char in itertools.cycle("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"):
            if done.is_set():
                break
            sys.stderr.write(f"\r{indent}{char} {in_progress_msg}")
            sys.stderr.flush()
            time.sleep(0.1)
        if elapsed[0] >= 60:
            mins, secs = divmod(int(elapsed[0]), 60)
            time_str = f"(took ~{mins} min and {secs} sec)"
        else:
            time_str = ""
        sys.stderr.write(f"\r{indent}✔ {complete_msg}{time_str}          \n")
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
