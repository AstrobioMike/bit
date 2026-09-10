import sys
import argparse
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg
from bit.modules.get_workflow import dl_wf


WORKFLOWS = {
    "amplicon":         "amplicon workflow",
    "genome-summarize": "genome-summarize workflow",
    "metagenomics":     "metagenomics workflow",
    "sra-download":     "sra-download workflow",
}


def build_parser(parent_subparsers=None):

    desc = """
        This is a helper program for downloading bit workflows.
        Workflow version is included with the downloaded workflow.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "get-workflow",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit get-workflow metagenomics`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    add_help(parser)

    add_version_arg(parser)

    workflow_subparsers = parser.add_subparsers(dest="workflow", required=True, metavar='')
    parser.subparsers = workflow_subparsers

    def add_workflow_common_args(group):
        group.add_argument(
            "-l",
            "--list-available-versions",
            help="Print the versions available for this workflow",
            action="store_true"
        )
        group.add_argument(
            "-w",
            "--wanted-version",
            metavar="VERSION",
            help="Specify the workflow version to download"
        )

    for workflow_name, workflow_desc in WORKFLOWS.items():

        workflow_parser = workflow_subparsers.add_parser(
            workflow_name,
            help=f"Download the {workflow_desc}",
            description=f"This subcommand downloads bit's {workflow_desc}.",
            epilog=f"Ex. usage: `bit get-workflow {workflow_name}`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

        workflow_optional = workflow_parser.add_argument_group("Optional Parameters")

        add_workflow_common_args(workflow_optional)

        add_help(workflow_optional)

        add_version_arg(workflow_optional)

    return parser

def main():

    parser = build_parser()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    dl_wf(args)
