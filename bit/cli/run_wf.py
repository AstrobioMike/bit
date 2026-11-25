import sys
import shlex
import argparse
from pathlib import Path
from bit.cli.common import (CustomRichHelpFormatter,
                            add_help,
                            add_version_arg,
                            add_common_snakemake_arguments)
from bit.modules.workflows import WORKFLOWS, load_workflow, invocation, RUN_GROUP

# one-liners for the `bit run-wf -h` listing (the full description of each is in its params.py)
WORKFLOW_SUMMARIES = {
    "metagenomics": "assembly-based processing of short-read metagenomes (QC through to annotation and MAGs)",
    "magrec": "MAG recovery from reads and already-made assemblies (mapping, binning, quality, taxonomy, dereplication)",
    "amplicon": "amplicon processing of Illumina 16S, 18S or ITS data (primer trimming through to ASVs and taxonomy)",
}

DETAILED_HELP_FLAGS = ("-H", "--detailed-help")


def build_parser(parent_subparsers=None, show_fine=False):

    desc = """
        This program runs the workflows packaged with bit. Run one with `-h` to see its
        main parameters, or with `-H` to see every setting it has.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "run-wf",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            epilog="Ex. usage: `bit run-wf metagenomics -r reads-dir/`",
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )

    add_help(parser)
    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="workflow", metavar="")
    parser.subparsers = subparsers

    for name in WORKFLOWS:
        add_workflow_parser(subparsers, name, show_fine, help=WORKFLOW_SUMMARIES[name])

    return parser


def add_workflow_parser(subparsers, name, show_fine=False, **parser_kwargs):
    """
    Builds the parser of one packaged workflow, as a subparser of `subparsers`
    or, when that is None, as a parser of its own.
    """

    workflow = load_workflow(name)
    # the workflows that are their own commands run as packaged: no saved copies or snakemake pass-through
    as_command = name not in WORKFLOWS

    desc = " ".join(workflow.description.split())
    if not show_fine:
        desc += " Run with `-H` to see the detailed help menu with every setting of the workflow."

    kwargs = dict(
        description=desc,
        epilog=f"Ex. usage: `{workflow.example}`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False,
    )
    if subparsers is not None:
        parser = subparsers.add_parser(name, **kwargs, **parser_kwargs)
    else:
        parser = argparse.ArgumentParser(**kwargs)

    # the workflow's own settings, from its parameter table; those not marked as main are still
    # registered when the concise menu is built (so they parse), just hidden from it. Settings in
    # the Run Parameters group are added further down, after the shared run options.
    groups = {}
    for p in workflow.params:
        if p.group == RUN_GROUP:
            continue
        visible = p.main or show_fine
        if visible and p.group not in groups:
            groups[p.group] = parser.add_argument_group(p.group)
            if p.group == "General Parameters":
                # -o leads the General section
                add_output_dir_argument(groups[p.group], name)
        target = groups[p.group] if visible else parser
        add_param_argument(target, p, visible)

    if "General Parameters" not in groups:
        groups["General Parameters"] = parser.add_argument_group("General Parameters")
        add_output_dir_argument(groups["General Parameters"], name)
    general = groups["General Parameters"]
    if workflow.ref_dbs_needed:
        general.add_argument(
            "--setup-dbs-only",
            help="Only download and set up the reference databases this run would use",
            action="store_true",
        )
    else:
        parser.set_defaults(setup_dbs_only=False)
    general.add_argument(
        *DETAILED_HELP_FLAGS,
        action="store_true",
        help="Show the detailed help, including every setting of the workflow",
    )
    add_help(general)
    add_version_arg(general)

    # Run Parameters: listed only in the detailed menu (-H); in the concise one they still parse, just unlisted
    def shown(text):
        return text if show_fine else argparse.SUPPRESS

    run = parser.add_argument_group(RUN_GROUP) if show_fine else parser
    run.add_argument(
        "-c",
        "--config-file",
        metavar="<FILE>",
        help=shown("Config file with settings for the workflow (as written by a previous run"
                   + ("" if as_command else " or by --get-workflow")
                   + "); settings given as flags take priority over it"),
        type=Path,
    )
    run.add_argument(
        "--memory",
        metavar="<GB>",
        help=shown("Max memory any jobs running simultaneously may use, in GB (default: all of this machine's memory); jobs are "
                   "only started while they fit, and a run with a step needing more than this stops before anything runs "
                   "(not used when the jobs go to a cluster through a snakemake executor or profile)"),
        type=float,
        default=None,
    )
    run.add_argument(
        "--conda-prefix",
        metavar="<DIR>",
        help=shown("Where the workflow's conda environments are created and re-used (default: ${CONDA_PREFIX}/envs)"),
        type=Path,
        default=None,
    )
    if as_command:
        parser.set_defaults(get_workflow=None, workflow_dir=None, snakemake_args=[])
    else:
        add_workflow_copy_arguments(run, shown)
    # the workflow's own run settings (e.g., the memory of each step), after the shared ones
    for p in workflow.params:
        if p.group == RUN_GROUP:
            add_param_argument(run, p, show_fine)

    snakemake = parser.add_argument_group("Snakemake Parameters")
    add_common_snakemake_arguments(snakemake)
    if not as_command:
        snakemake.add_argument(
            "snakemake_args",
            nargs="*",
            metavar="-- <ARGS>",
            help="Anything given after a lone `--` is passed to snakemake as is",
        )

    return parser


def add_workflow_copy_arguments(run, shown):
    run.add_argument(
        "--get-workflow",
        metavar="<DIR>",
        help=shown("Don't run anything; write a standalone copy of the workflow (Snakefile, rules, environments, "
                   "scripts, and a config.yaml holding the settings given here) into this directory (new or empty) for you to modify as wanted"),
        type=Path,
    )
    run.add_argument(
        "--workflow-dir",
        metavar="<DIR>",
        help=shown("Run the copy of the workflow in this directory (made with --get-workflow) rather than the packaged one; "
                   "its config.yaml is used unless --config-file is given"),
        type=Path,
    )


def add_output_dir_argument(target, name):
    target.add_argument(
        "-o",
        "--output-dir",
        metavar="<DIR>",
        help=f"Directory for the outputs, also where the workflow is run (default: '{name}-outputs/')",
        type=Path,
        default=Path(f"{name}-outputs"),
    )


def add_param_argument(target, p, visible):
    names = ([p.short] if p.short else []) + [p.flag]
    help_text = " ".join(p.help.split())
    # (no note for an empty default: the help says what happens without the setting)
    if p.default not in (None, "") and p.kind != "flag":
        # merged into a closing parenthesis the help already ends with
        help_text = help_text[:-1] + f"; default: {p.default})" if help_text.endswith(")") else help_text + f" (default: {p.default})"
    # (argparse formats help text with %, so a literal one is doubled)
    help_text = help_text.replace("%", "%%") if visible else argparse.SUPPRESS

    # argparse defaults are all None so we can tell what was given on the command line
    # (defaults and config-file values are filled in afterwards, see modules/workflows.py)
    if p.kind == "flag":
        target.add_argument(*names, dest=p.name, action="store_const", const=True, default=None, help=help_text)
    else:
        metavar = p.metavar or {"file": "<FILE>", "dir": "<DIR>", "path": "<PATH>", "paths": "<PATH>"}.get(p.kind, "<STR>")
        kwargs = dict(dest=p.name, default=None, metavar=metavar, help=help_text)
        if p.kind == "paths":
            kwargs["nargs"] = "+"
        if p.kind == "value":
            kwargs["type"] = p.type
            if p.choices:
                kwargs["choices"] = p.choices
                kwargs.pop("metavar")
        target.add_argument(*names, **kwargs)


def attach_dash_values(parser, argv):
    """
    Joins a value that starts with "-" to its option (`-x -assembly.fasta` -> `-x=-assembly.fasta`), which argparse would
    otherwise take for another option (quoting doesn't help: the shell removes the quotes). Only for options that take a
    single value, and only when the value isn't one of the parser's own options; nothing after `--` is touched.
    """

    options = parser._option_string_actions
    fixed, i = [], 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--":
            fixed += argv[i:]
            break
        action = options.get(arg)
        if isinstance(action, argparse._StoreAction) and action.nargs is None and i + 1 < len(argv) \
                and argv[i + 1].startswith("-") and argv[i + 1] not in options and argv[i + 1] != "--":
            fixed.append(f"{arg}={argv[i + 1]}")
            i += 2
            continue
        fixed.append(arg)
        i += 1
    return fixed


def main():

    argv = sys.argv[1:]
    wanted = next((a for a in argv if a in WORKFLOWS), None)

    if wanted and any(flag in argv for flag in DETAILED_HELP_FLAGS):
        build_parser(show_fine=True).subparsers.choices[wanted].print_help(sys.stderr)
        sys.exit(0)

    parser = build_parser()

    if not argv or wanted is None:
        if argv and argv[0] in ("-v", "--version"):
            parser.parse_args()
        parser.print_help(sys.stderr)
        sys.exit(0)

    # a workflow with nothing else given gets its help menu
    if argv == [wanted]:
        parser.subparsers.choices[wanted].print_help(sys.stderr)
        sys.exit(0)

    position = argv.index(wanted) + 1
    argv = argv[:position] + attach_dash_values(parser.subparsers.choices[wanted], argv[position:])
    args = parser.parse_args(argv)
    run_parsed(args.workflow, args, parser, argv)


def command_build_parser(name, parent_subparsers=None, show_fine=False):
    """The parser of a workflow that is its own bit command (see WORKFLOW_COMMANDS)."""
    return add_workflow_parser(parent_subparsers, name, show_fine)


def command_main(name):
    """Entry point of a workflow that is its own bit command."""

    argv = sys.argv[1:]

    if any(flag in argv for flag in DETAILED_HELP_FLAGS):
        command_build_parser(name, show_fine=True).print_help(sys.stderr)
        sys.exit(0)

    parser = command_build_parser(name)

    if not argv:
        parser.print_help(sys.stderr)
        sys.exit(0)

    argv = attach_dash_values(parser, argv)
    run_parsed(name, parser.parse_args(argv), parser, argv)


def run_parsed(name, args, parser, argv):

    from bit.modules.workflows import resolve_settings, get_workflow, run_workflow
    from bit.modules.general import report_message

    workflow = load_workflow(name)
    # the command as typed; every setting it resolves to is recorded in the config file written for the run
    command = shlex.join([*invocation(name).split(), *[a for a in argv if a != name or name not in WORKFLOWS]])

    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater")
    if args.memory is not None and args.memory <= 0:
        parser.error("--memory must be greater than 0")

    config_file = args.config_file
    if config_file is None and args.workflow_dir is not None:
        config_file = Path(args.workflow_dir) / "config.yaml"

    cli_values = {p.name: getattr(args, p.name) for p in workflow.params if getattr(args, p.name) is not None}
    settings = resolve_settings(workflow.params, cli_values, config_file,
                                only_require=workflow.db_setup_settings if args.setup_dbs_only else None)

    if args.get_workflow is not None:
        dest = get_workflow(workflow, settings, args.get_workflow, command)
        report_message(f"The {workflow.name} workflow was written to: {dest}/", color="green")
        print("    Modify it as wanted, then run it with either of:")
        print(f"      {invocation(workflow.name)} --workflow-dir {dest}")
        print(f"      snakemake --software-deployment-method conda --conda-prefix ${{CONDA_PREFIX}}/envs --cores <N> -p   (from within {dest}/)\n")
        sys.exit(0)

    run_workflow(workflow, settings, args, command)
