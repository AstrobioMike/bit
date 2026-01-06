import time
import sys
import os
import subprocess
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help
from bit.modules.general import report_message


def build_parser():

    desc = """
        This is a program for downloading test data files.
        For version info, run `bit-version`.
        """

    parser = argparse.ArgumentParser(
        description=desc,
        epilog="Ex. usage: bit-get-test-data metagenomics",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    required = parser.add_argument_group("REQUIRED PARAMETERS")

    required.add_argument(
        "datatype",
        choices=["metagenomics"],
        help="The first positional argument should be what type of test data you'd like to download",
    )

    add_help(parser)

    return parser


def main():

    parser = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_test_data(args)


# def dl_test_data(args):

#     if args.datatype == "metagenomics":

#         report_message("Downloading and unpacking 2 paired-end Illumina metagenomics test samples (4 files, ~800 MB total; they are kinda large for test data so MAGs can be recovered):")
#         print("")

#         # getting the metagenomics test data
#         os.system("curl -L -o test-metagenomics-reads.zip https://figshare.com/ndownloader/files/46096083")

#         # extracting
#         os.system("unzip -qo test-metagenomics-reads.zip")

#         # removing archive
#         os.system("rm test-metagenomics-reads.zip")

#         report_message("Pulled metagenomics (Illumina) reads for two test samples from here:", "green")
#         print("    https://figshare.com/account/projects/203736/articles/25750935\n")

#     else:

#         report_message("The data type you requested is not currently available.", "red")

#         print("\n    Please check the currently available data types with 'bit-get-test-data --help'\n")

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
