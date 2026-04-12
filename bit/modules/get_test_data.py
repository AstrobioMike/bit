import time
import os
import sys
import subprocess
from bit.modules.general import report_message


def dl_test_data(args):

    max_retries = 3
    success = False

    if args.datatype == "metagenome":

        page_url="https://github.com/AstrobioMike/test-data/releases/tag/test-metagenomics-reads-v1"
        data_url="https://github.com/AstrobioMike/test-data/releases/download/test-metagenomics-reads-v1/test-metagenomics-reads.zip"
        dest="test-metagenomics-reads.zip"
        start_message=("Downloading and unpacking a tiny paired-end Illumina metagenomics test sample (2 files, ~15 MB total) from here: ")
        end_message=("Successfully pulled in a tiny mock metagenome.")

    elif args.datatype == "genome":

        page_url="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/"
        data_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
        dest="GCF_000005845.2.fasta.gz"
        start_message=("Downloading an E. coli genome (GCF_000005845.2) from here: ")
        end_message=("Successfully pulled in the E. coli genome (GCF_000005845.2.fasta).")


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

    if args.datatype in ["metagenome"]:
        try:
            report_message("Unpacking data...")
            subprocess.run(["unzip", "-qo", dest], check=True)
            os.remove(dest)
        except subprocess.CalledProcessError:
            report_message("Failed to unzip the downloaded file for some reason :(", "red", trailing_newline=True)
            sys.exit(1)

    if args.datatype in ["genome"]:
        try:
            report_message("Unpacking data...")
            subprocess.run(["gunzip", "-f", dest], check=True)
        except subprocess.CalledProcessError:
            report_message("Failed to gunzip the downloaded file for some reason :(", "red", trailing_newline=True)
            sys.exit(1)

    report_message(end_message, "green", trailing_newline=True)
