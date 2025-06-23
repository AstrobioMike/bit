from rich_argparse import RichHelpFormatter

class CustomRichHelpFormatter(RichHelpFormatter):
    def start_section(self, heading):
        if heading == "positional arguments":
            heading = "AVAILABLE SUBCOMMANDS"
        elif heading == "options":
            heading = "OPTIONAL PARAMETERS"
        super().start_section(heading)


def add_common_snakemake_arguments(group):
    group.add_argument(
        "-j",
        "--jobs",
        help = "Number of jobs to run in parallel (default: 10)",
        metavar = "<NUM>",
        action = "store",
        default = 10,
        type = int,
    )

    group.add_argument(
        "--rerun-incomplete",
        help = "Re-run all jobs the output of which is recognized as incomplete",
        action = "store_true",
    )

    group.add_argument(
        "--dry-run",
        help = "Do not execute anything, only show what would be done",
        action = "store_true",
    )
