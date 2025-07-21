import sys
import shlex
import argparse
from pathlib import Path
from rich_argparse import RichHelpFormatter

class CustomRichHelpFormatter(RichHelpFormatter):
    def start_section(self, heading):
        if heading == "positional arguments":
            heading = "AVAILABLE SUBCOMMANDS"
        elif heading == "options":
            heading = "OPTIONAL PARAMETERS"
        super().start_section(heading)


def add_help(group):
    group.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit"
    )


def add_common_snakemake_arguments(group):
    group.add_argument(
        "-j",
        "--jobs",
        help = "Max number of jobs to run in parallel (default: 10)",
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

def reconstruct_invocation(parser, args):
    """
    This is a helper for capturing the executed command-line (including default args) for logging.
    """
    cmd = [sys.executable, sys.argv[0]]

    # handling if there is a subcommand
    sub = getattr(args, "subcommand", None)
    if sub:
        cmd.append(sub)

    # gathering all actions
    all_actions = list(parser._actions)
    # including subparsers if they exist
    for a in parser._actions:
        if isinstance(a, argparse._SubParsersAction):
            subpar_action = a
            break
    else:
        subpar_action = None

    if sub and subpar_action:
        all_actions += subpar_action.choices[sub]._actions

    # walking through and getting settings for each item
    for action in all_actions:
        # skip help and the subparsers action itself
        if isinstance(action, (argparse._HelpAction,
                               argparse._SubParsersAction)):
            continue

        # named parameters handled here
        if action.option_strings:
            val = getattr(args, action.dest, None)
            if isinstance(action, argparse._StoreTrueAction):
                if val: cmd.append(action.option_strings[-1])
                continue
            if isinstance(action, argparse._StoreFalseAction):
                if not val: cmd.append(action.option_strings[-1])
                continue

            flag = next((o for o in action.option_strings if o.startswith("--")),
                        action.option_strings[0])

            # this little bit deals with cases where the value is a list or tuple
            if isinstance(val, (list, tuple)):
                for v in val:
                    cmd.append(flag)
                    cmd.append(str(v))
            else:
                cmd.append(flag)
                cmd.append(str(val))

        # positional parameters handled here
        else:
            val = getattr(args, action.dest, None)
            if val is None:
                continue
            if isinstance(val, (list, tuple)):
                cmd.extend(str(v) for v in val)
            else:
                cmd.append(str(val))

    return shlex.join(cmd)
