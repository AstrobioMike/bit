from rich_argparse import RichHelpFormatter

class CustomRichHelpFormatter(RichHelpFormatter):
    def start_section(self, heading):
        if heading == "positional arguments":
            heading = "AVAILABLE SUBCOMMANDS"
        elif heading == "options":
            heading = "OPTIONAL PARAMETERS"
        super().start_section(heading)

