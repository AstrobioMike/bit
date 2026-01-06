import time
import os
import sys
import subprocess
from bit.modules.general import report_message


def dl_test_data(args):

    max_retries = 3
    success = False

    if args.datatype == "metagenomics":

        page_url="https://figshare.com/account/projects/203736/articles/25750935"
        data_url = "https://figshare.com/ndownloader/files/46096083"
        dest = "test-metagenomics-reads.zip"
        start_message=("Downloading and unpacking 2 paired-end Illumina metagenomics test samples (4 files, ~800 MB total) from here: ")
        end_message=("Successfully pulled metagenomics (Illumina) reads for two test samples.")

    report_message(start_message, trailing_newline=True)
    print(f"    {page_url}\n\n")

    curl_cmd = [
        "curl",
        "-L",
        "--connect-timeout", "30",
        "-o", dest,
        data_url
    ]

    for i in range(max_retries):

        result = subprocess.run(curl_cmd)

        if result.returncode == 0:
            success = True
            break
        else:
            if i < max_retries - 1:
                report_message(f"Download attempt {i+1} failed. Trying again...", "yellow", initial_indent = "    ", trailing_newline=True)
                time.sleep(1)

    if not success:
        time.sleep(1)
        report_message("Failed to download data after multiple attempts :( Maybe the connection is bad or the site is blocked for you?", "red", trailing_newline=True)
        sys.exit(1)

    try:
        report_message("Unpacking data...")
        subprocess.run(["unzip", "-qo", dest], check=True)
        os.remove(dest)
    except subprocess.CalledProcessError:
        report_message("Failed to unzip the downloaded file for some reason :(", "red", trailing_newline=True)
        sys.exit(1)

    report_message(end_message, "green", trailing_newline=True)
