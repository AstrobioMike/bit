from bit.cli.run_wf import command_build_parser, command_main

# a packaged workflow run as its own command; its settings are in bit/smk/workflows/characterize_genomes/params.py
NAME = "characterize-genomes"


def build_parser(parent_subparsers=None):
    return command_build_parser(NAME, parent_subparsers)


def main():
    command_main(NAME)
