import os
import sys
import time
import gzip
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests # type: ignore
from dataclasses import dataclass
from tqdm import tqdm # type: ignore
from bit.modules.general import (color_text, check_files_are_found,
                                 attempt_to_make_dir)
from bit.modules.ncbi.parse_ncbi_assembly_summary import parse_ncbi_assembly_summary
from bit.modules.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data


TRANSIENT_STATUS = {429, 500, 502, 503, 504}

max_threads=20
max_retries=20

def dl_ncbi_assemblies(args):

    preflight_checks(args)

    run_data = setup(args)

    run_data = parse_main_assembly_table(run_data)

    summarize_search(run_data)

    run_data = download_assemblies(run_data)

    report_finish(run_data)


def preflight_checks(args):

    check_files_are_found([args.wanted_accessions])

    if args.output_dir and not os.path.exists(args.output_dir):
        attempt_to_make_dir(args.output_dir)

    get_ncbi_assembly_data(quiet=True)


def setup(args):

    with open(args.wanted_accessions, "r") as f:
        wanted_accs = [line.strip() for line in f if line.strip()]

    run_data = RunData(wanted_format=args.format,
        num_jobs=args.jobs,
        output_dir=args.output_dir,
        wanted_accs=wanted_accs,
        num_wanted=len(wanted_accs),
        ncbi_sub_table_path=Path(args.output_dir) / "wanted-ncbi-accessions-info.tsv",
        not_found_path=Path(args.output_dir) / "ncbi-accessions-not-found.txt",
        not_downloaded_path=Path(args.output_dir) / "ncbi-accessions-not-downloaded.tsv",
        quiet=getattr(args, "quiet", False)
        )

    return run_data


def parse_main_assembly_table(run_data):

    if not run_data.quiet:
        print(color_text(f"\n    Targeting {run_data.num_wanted} accession(s) in {run_data.wanted_format} format...\n", "yellow"))

    assembly_summary_file = Path(os.environ['NCBI_assembly_data_dir']) / "ncbi-assembly-info.tsv"

    run_data = parse_ncbi_assembly_summary(assembly_summary_file, run_data)

    return run_data


@dataclass
class RunData:
    wanted_format: str = None
    num_jobs: int = 10
    output_dir: str = None
    wanted_accs: list = None
    num_wanted: int = 0
    num_found: int = 0
    num_not_found: int = 0
    num_downloaded: int = 0
    num_skipped: int = 0
    num_not_downloaded: int = 0
    ncbi_sub_table_path: str = None
    not_found_path: str = None
    not_downloaded_path: str = None
    quiet: bool = False


def summarize_search(summary):

    if summary.num_found != summary.num_wanted:
        if summary.num_found > 0:
            print(color_text(f"{' ' * 34}NOTICE", "orange"))
            print(f"      {summary.num_not_found} accession(s) not found at NCBI. They may be invalid or suppressed.")
            print(f"          See '{summary.not_found_path}'.\n")
        else:
            print(color_text(f"{' ' * 34}NOTICE", "orange"))
            print(f"      None of the {summary.num_wanted} target accession(s) were found at NCBI...")
            print(f"      This is kinda weird. Are the inputs assembly accessions?\n")
            os.remove(summary.not_found_path)
            sys.exit(1)

        if not summary.quiet:
            print(f"    Remaining total targets: {summary.num_found}\n")


def valid_gzip(path):
    """
    integrity check for a gzipped file (only used for retried runs)
    """
    if str(path).endswith(".gz"):
        try:
            with gzip.open(path, "rb") as fh:
                while fh.read(1024 * 1024):
                    pass
            return True
        except (OSError, EOFError):
            return False
    return True


def sleep_backoff(attempt, resp=None):
    """
    sleep before a retry. Honors an NCBI-provided Retry-After header when present,
    otherwise falls back to exponential backoff with jitter
    """
    if resp is not None:
        retry_after = resp.headers.get("Retry-After")
        if retry_after:
            try:
                time.sleep(float(retry_after))
                return
            except ValueError:
                pass
    time.sleep((2 ** (attempt - 1)) + random.uniform(0, 1))


def download_one(target_link, local_dest, retries=max_retries):

    local_path = Path(local_dest)

    if local_path.exists() and local_path.stat().st_size > 0 and valid_gzip(local_path):
        return (local_dest, None, "skipped")

    local_path.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(target_link, stream=True, timeout=60)

            if resp.status_code == 404:
                return (local_dest, "Not available in requested format (404)", "failed")

            if resp.status_code in TRANSIENT_STATUS:
                if attempt == retries:
                    return (local_dest, f"HTTP {resp.status_code} after {retries} attempts", "failed_transient")
                sleep_backoff(attempt, resp)
                continue

            resp.raise_for_status()

            content_type = resp.headers.get("Content-Type", "")
            if "xml" in content_type.lower() or "html" in content_type.lower():
                if attempt == retries:
                    return (local_dest, f"NCBI returned an error page after {retries} attempts", "failed_transient")
                sleep_backoff(attempt)
                continue

            with open(local_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1024 * 64):
                    fh.write(chunk)

            if local_path.stat().st_size == 0:
                local_path.unlink(missing_ok=True)
                if attempt == retries:
                    return (local_dest, "Downloaded file was empty", "failed_transient")
                sleep_backoff(attempt)
                continue

            return (local_dest, None, "downloaded")

        except (requests.RequestException, OSError) as e:
            if attempt == retries:
                local_path.unlink(missing_ok=True)
                return (local_dest, str(e), "failed_transient")
            sleep_backoff(attempt)


def run_download_pass(targets, run_data, desc="Progress"):
    """
    runs one pooled download pass over a list of (target_link, local_dest) tuples.
    returns (permanent_failures, transient_failures, num_skipped) where each
    failure list holds (target_link, local_dest, error) so transient ones can be retried
    """
    permanent = []
    transient = []
    num_skipped = 0

    link_by_dest = {dest: link for link, dest in targets}

    with ThreadPoolExecutor(max_workers=min(run_data.num_jobs, max_threads)) as pool:
        futures = {
            pool.submit(download_one, link, dest): dest
            for link, dest in targets
        }

        if run_data.quiet:
            desc_buffer = "    "
            ncols = 78
        else:
            desc_buffer = "      "
            ncols = 70
        with tqdm(total=len(targets), desc=f"{desc_buffer}{desc}", unit=" file", ncols=ncols) as pbar:
            for future in as_completed(futures):
                dest, error, status = future.result()
                if status == "failed":
                    permanent.append((link_by_dest[dest], dest, error))
                elif status == "failed_transient":
                    transient.append((link_by_dest[dest], dest, error))
                elif status == "skipped":
                    num_skipped += 1
                pbar.update(1)

    return permanent, transient, num_skipped


def download_assemblies(run_data):

    if not run_data.quiet:
        print(color_text("    Downloading assemblies...\n", "yellow"))

    targets = []
    with open(run_data.ncbi_sub_table_path, "r") as f:
        header = f.readline().strip().split("\t")
        link_idx = header.index("target_link")
        dest_idx = header.index("local_destination")
        for line in f:
            fields = line.strip().split("\t")
            targets.append((fields[link_idx], fields[dest_idx]))

    permanent, transient, num_skipped = run_download_pass(targets, run_data)

    # second pass on transient-only failures
    if transient:
        retry_targets = [(link, dest) for link, dest, _ in transient]
        print(color_text(f"\n    {len(transient)} file(s) failed with transient error messages, doing another pass", "yellow"))
        print(color_text(f"    to see if we can grab them...\n", "yellow"))

        time.sleep(3)
        retry_permanent, retry_transient, retry_skipped = run_download_pass(
            retry_targets, run_data, desc="Progress"
        )
        num_skipped += retry_skipped
        # anything still failing after the retry is final, regardless of category
        permanent.extend(retry_permanent)
        permanent.extend(retry_transient)

    failed = [(dest, error) for _, dest, error in permanent]

    if failed:
        with open(run_data.not_downloaded_path, "w") as fh:
            fh.write("accession\terror\n")
            for dest, error in failed:
                acc = Path(dest).stem.split(".")[0]
                fh.write(f"{acc}\t{error}\n")
        run_data.num_not_downloaded = len(failed)
    else:
        Path(run_data.not_downloaded_path).unlink(missing_ok=True)

    run_data.num_skipped = num_skipped
    run_data.num_downloaded = len(targets) - len(failed)

    return run_data


def report_finish(run_data):

    skipped_note = ""
    if run_data.num_skipped > 0:
        skipped_note = f" ({run_data.num_skipped} already present, skipped)"

    if run_data.num_downloaded == run_data.num_wanted:
        if not run_data.quiet:
            print(color_text(f"\n    All {run_data.num_wanted} file(s) downloaded successfully!{skipped_note}\n", "green"))

    elif run_data.num_downloaded == run_data.num_found:
        if not run_data.quiet:
            print(color_text(f"\n    All {run_data.num_found} found file(s) downloaded successfully!{skipped_note}\n", "yellow"))

    elif run_data.num_not_downloaded > 0:
        print(color_text(f"\n{' ' * 34}NOTICE", "orange"))
        print(f"      {run_data.num_not_downloaded} file(s) failed to download from NCBI. They may not be available")
        print(f"      in the requested format, or it may have been a transient problem.")
        print(f"          See '{run_data.not_downloaded_path}'.")

        if run_data.num_downloaded > 0:
            if not run_data.quiet:
                print(color_text(f"\n\n    The remaining {run_data.num_downloaded} found file(s) downloaded successfully.{skipped_note}\n", "yellow"))
        else:
            print(color_text(f"\n\n    No files were successfully downloaded...{skipped_note}\n", "orange"))
            sys.exit(1)
