import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests # type: ignore
from dataclasses import dataclass
from tqdm import tqdm # type: ignore
from bit.modules.general import (color_text, check_files_are_found,
                                 attempt_to_make_dir, report_message)
from bit.modules.ncbi.parse_assembly_summary import parse_assembly_summary
from bit.modules.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data


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
        not_downloaded_path=Path(args.output_dir) / "ncbi-accessions-not-downloaded.tsv"
        )

    return run_data


def parse_main_assembly_table(run_data):

    print(color_text(f"\n    Targeting {run_data.num_wanted} accession(s) in {run_data.wanted_format} format...", "yellow"))

    assembly_summary_file = Path(os.environ['NCBI_assembly_data_dir']) / "ncbi-assembly-info.tsv"

    run_data = parse_assembly_summary(assembly_summary_file, run_data)

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
    num_not_downloaded: int = 0
    ncbi_sub_table_path: str = None
    not_found_path: str = None
    not_downloaded_path: str = None


def summarize_search(summary):

    if summary.num_found != summary.num_wanted:
        if summary.num_found > 0:
            print(color_text("\n\n                      NOTICE", "orange"))
            print(f"        {summary.num_not_found} accession(s) not found at NCBI.")
            print(f"        See '{summary.not_found_path}'.\n\n")
        else:
            print(color_text("\n\n                      NOTICE", "orange"))
            print(f"        None of the {summary.num_wanted} target accession(s) were found at NCBI...")
            print(f"        This is weird. Are the inputs assembly accessions?\n")
            os.remove(summary.not_found_path)
            exit(1)

        print(f"    Remaining total targets: {summary.num_found}")


def download_one(target_link, local_dest, retries=3):

    local_path = Path(local_dest)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(target_link, stream=True, timeout=60)
            resp.raise_for_status()

            # catch for if we got an XML/HTML error page instead of the file
            content_type = resp.headers.get("Content-Type", "")
            if "xml" in content_type.lower() or "html" in content_type.lower():
                return (local_dest, f"NCBI returned an error page (attempt {attempt}/{retries})")

            with open(local_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1024 * 64):
                    fh.write(chunk)

            # check we actually got something
            if local_path.stat().st_size == 0:
                local_path.unlink(missing_ok=True)
                return (local_dest, "Downloaded file was empty")

            return (local_dest, None)

        except (requests.RequestException, OSError) as e:
            if attempt == retries:
                local_path.unlink(missing_ok=True)
                return (local_dest, str(e))


def download_assemblies(run_data):

    print(color_text("\n    Downloading assemblies...\n", "yellow"))

    # reading the table into a list of (target_link, local_dest) tuples
    targets = []
    with open(run_data.ncbi_sub_table_path, "r") as f:
        header = f.readline().strip().split("\t")
        link_idx = header.index("target_link")
        dest_idx = header.index("local_destination")
        for line in f:
            fields = line.strip().split("\t")
            targets.append((fields[link_idx], fields[dest_idx]))

    failed = []

    with ThreadPoolExecutor(max_workers=min(run_data.num_jobs, 20)) as pool:
        futures = {
            pool.submit(download_one, link, dest): dest
            for link, dest in targets
        }


        with tqdm(total=len(targets), desc="      Progress", unit="file", ncols=70) as pbar:
            for future in as_completed(futures):
                dest, error = future.result()
                if error:
                    failed.append((dest, error))
                pbar.update(1)

    if failed:
        with open(run_data.not_downloaded_path, "w") as fh:
            fh.write("accession\terror\n")
            for dest, error in failed:
                acc = Path(dest).stem.split(".")[0]
                fh.write(f"{acc}\t{error}\n")
        run_data.num_not_downloaded = len(failed)

    run_data.num_downloaded = len(targets) - len(failed)

    return run_data


def report_finish(run_data):

    if run_data.num_downloaded == run_data.num_wanted:
        print(color_text(f"\n    All {run_data.num_wanted} file(s) downloaded successfully!\n", "green"))

    elif run_data.num_downloaded == run_data.num_found:
        print(color_text(f"\n    All {run_data.num_found} found file(s) downloaded successfully!\n", "yellow"))

    elif run_data.num_not_downloaded > 0:
        print(color_text("\n\n                      NOTICE", "orange"))
        print(f"        {run_data.num_not_downloaded} file(s) failed to download from NCBI.")
        print(f"        They may not be available in the requested format.")
        print(f"        See '{run_data.not_downloaded_path}'.\n")

        print(color_text(f"\n    The remaining {run_data.num_downloaded} found file(s) downloaded successfully.\n", "yellow"))
