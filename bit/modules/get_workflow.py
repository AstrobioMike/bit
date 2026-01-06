import os
import requests
import re
from packaging import version
from bit.modules.general import report_message


base_repo = "https://github.com/astrobiomike/bit/"
base_download_link = base_repo + "releases/download/"
releases_page_link = "https://api.github.com/repos/astrobiomike/bit/releases"
html_releases_page = base_repo + "releases"

workflow_dict = {
    "metagenomics": { "basename": "metagenomics-wf" },
    "genome-summarize": { "basename": "genome-summarize-wf" },
    "sra-download": { "basename": "sra-download-wf" }
}


def get_all_releases():
    releases = []
    headers = {"Accept": "application/vnd.github.v3+json"}
    params = {"per_page": 100, "page": 1}

    while True:
        # resp = requests.get(releases_page_link, headers=headers, params=params)
        resp = requests.get(releases_page_link, headers=headers, params=params.copy())
        resp.raise_for_status()
        page_data = resp.json()
        if not page_data:
            break
        releases.extend(page_data)
        params["page"] += 1

    return releases


def get_versions_available(target_wf):

    """
    This finds which versions are available on the release page(s), returns them as a list,
    and creates a dictionary holding the versions as keys and full-workflow names as values.
    """

    all_releases = get_all_releases()

    # getting all links available
    links = [item['html_url'] for item in all_releases if 'html_url' in item]

    # starting empty list of available versions
    available_versions = []

    # starting empty dictionary to hold version: full_workflow_basename
    dict_of_versions_and_basenames = {}

    for link in links:

        # getting just those for the requested workflow
        if re.search(target_wf, link):

            # getting just those that hold a tag which should help us get the version
            if re.search("releases/tag", link):

                # cutting out just the version info
                link_basename = os.path.basename(link)
                version_info_dash_index = link_basename.rfind('-')
                link_version = link_basename[version_info_dash_index + 1:].replace("v", "")

                # adding to list
                available_versions.append(link_version)

                # adding to dictionary of verisons
                dict_of_versions_and_basenames[link_version] = link_basename


    return(available_versions, dict_of_versions_and_basenames)


def check_version_available(workflow, wanted_version, available_versions):

    """ this checks if the requested version was found on the releases page """

    if wanted_version not in available_versions:

        report_message(f"The requested version, {wanted_version}, is not available.")
        print(f"\n    Those available for the {workflow} workflow include:\n")

        for version in available_versions:

            print(f"        {version}")

        print("\n  Exiting for now.\n")
        exit(1)


def check_if_dir_already_exists(dir_name):

    """ this checks if the output directory already exists, so as to not overwrite anything """

    if os.path.exists(dir_name):

        report_message(f"The output directory '{dir_name}/' already exists and we don't want to overwrite anything.")
        print("\n    Please rename or remove it first.")
        print("\n  Exiting for now.\n")
        exit(1)


def download_and_unzip(workflow, target_link):

    """ this downloads and unzips the target workflow """

    # python module zipfile does not retain permissions (see https://github.com/python/cpython/issues/59999 and https://github.com/python/cpython/pull/32289)
    # this was messing up places where files were expected to come in as executables, this was the old way
    # # downloading
    # target = requests.get(target_link)

    # # extracting zip
    # zip = zipfile.ZipFile(BytesIO(target.content))
    # zip.extractall()

    # now doing it with the system zip (that's part of bit) to retain permissions

    # getting base filename
    base_zip = os.path.basename(target_link)

    # downloading
    target = requests.get(target_link)

    # writing out
    with open(base_zip, "wb") as downloaded_zip:
        downloaded_zip.write(target.content)

    # extracting zip and removing archive
    os.system(f"unzip -q {base_zip}")
    os.remove(base_zip)

    downloaded_basename = base_zip.replace(".zip", "")

    report_message(f"The {workflow} workflow was downloaded to '{downloaded_basename}/'", "green")

    print(f"\n    It was pulled from this releases page:\n         {html_releases_page}\n")


def dl_wf(args):

    """ main download function """

    # getting the specific filename/directory name for the requested workflow
    wf_basename = workflow_dict[args.workflow]["basename"]

    # getting which versions are available for download
    available_versions, dict_of_versions_and_basenames = get_versions_available(wf_basename)

    # sorting (so we can grab the latest if needed)
    available_versions = sorted(available_versions, key = lambda x: version.Version(x), reverse = True)

    # just reporting available versions if requested
    if args.list_available_versions:

        print(f"\n    Versions available for the {args.workflow} workflow include:\n")

        for each_version in available_versions:

            print(f"        {each_version}")

        print(f"\n    You can also find these at the releases page:\n")
        print(f"        {html_releases_page}\n")

        exit(0)

    if args.wanted_version:

        # checking specified version exists and exiting if not
        check_version_available(args.workflow, args.wanted_version, available_versions)

        full_target_name = dict_of_versions_and_basenames[args.wanted_version]
        # full_target_name = wf_basename + "_" + args.wanted_version

    else:

        # full_target_name = wf_basename + "_" + available_versions[0]
        full_target_name = dict_of_versions_and_basenames[available_versions[0]]

    # building full download link
    full_link = base_download_link + full_target_name + "/bit-" + full_target_name + ".zip"

    # checking if output directory exists already, and exiting if yes
    check_if_dir_already_exists("bit-" + full_target_name)

    download_and_unzip(args.workflow, full_link)
