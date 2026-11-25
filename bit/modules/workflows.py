"""
Shared machinery behind `bit run-wf` and the workflows that are their own bit commands.

Each workflow lives in bit/smk/workflows/<name>/ with a Snakefile and a params.py.
The params.py holds one table (PARAMS) listing every setting of the workflow once;
the command-line flags, both help menus (-h / -H), the config file and its checks
are all generated from that table here, so they cannot drift apart.
"""

import os
import re
import shlex
import subprocess
import shutil
import textwrap
import importlib
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
import yaml # type: ignore
from bit.modules.general import get_package_path, report_message, notify_premature_exit, color_text
from bit.modules.snakemake import build_snakemake_cmd, run_snakemake, planned_conda_envs, conda_env_urls


# cli name -> directory under bit/smk/workflows/
# run with `bit run-wf <name>`
WORKFLOWS = {
    "metagenomics": "metagenomics",
    "magrec": "magrec",
    "amplicon": "amplicon",
}
# run as their own commands, `bit <name>`
WORKFLOW_COMMANDS = {
    "characterize-genomes": "characterize_genomes",
    "dl-sra-reads": "dl_sra_reads",
    "classify-seqs": "classify_seqs",
}
WORKFLOW_DIRS = {**WORKFLOWS, **WORKFLOW_COMMANDS}


def invocation(name):
    """How a packaged workflow is called on the command line."""
    return f"bit run-wf {name}" if name in WORKFLOWS else f"bit {name}"


RUN_GROUP = "Run Parameters"

@dataclass(frozen=True)
class Param:
    """
    One workflow setting.

    name:     the key in the config file; the flag is the same with hyphens (--like-this)
    default:  default value (None means there is none, see `required`)
    help:     help text, used in the help menus and as the comment in the config file
    group:    heading it is listed under (RUN_GROUP -> the Run Parameters section, only in the detailed help)
    main:     True -> shown in the concise help (-h); False -> only in the detailed help (-H)
    kind:     "value", "flag" (True/False switch), "file" or "dir" (must exist; stored as absolute path),
              or "path" (stored as absolute path, need not exist)
    type:     python type of a "value" (str, int, float)
    """
    name: str
    default: object
    help: str
    group: str = "General Parameters"
    main: bool = False
    short: str = None
    kind: str = "value"
    type: type = str
    choices: tuple = None
    metavar: str = None
    required: bool = False

    @property
    def flag(self):
        return "--" + self.name.replace("_", "-")


@dataclass(frozen=True)
class Workflow:
    name: str                # as given on the command line
    directory: Path          # its directory in the package
    description: str
    params: tuple
    shared_files: tuple      # files under bit/smk/ it uses (e.g., "rules/common.smk", "envs/x.yaml")
    example: str
    outputs_note: str
    db_setup_settings: tuple # the settings that matter to --setup-dbs-only
    ref_dbs_needed: object   # settings -> the reference databases a run uses (bit.modules.ref_dbs), or None for none
    ref_dbs_in_output_dir: bool  # without a reference-database directory, they go in the output directory
    runs_gtdbtk: object      # settings -> whether a run uses GTDB-Tk (where it runs), or None for never
    downloads_from: dict     # rule -> the URLs its jobs download from (checked before a run with such jobs starts)


def load_workflow(name):
    module = importlib.import_module(f"bit.smk.workflows.{WORKFLOW_DIRS[name]}.params")
    return Workflow(
        name=name,
        directory=Path(get_package_path(f"smk/workflows/{WORKFLOW_DIRS[name]}")),
        description=module.DESCRIPTION,
        params=tuple(module.PARAMS),
        shared_files=tuple(module.SHARED_FILES),
        example=module.EXAMPLE,
        outputs_note=getattr(module, "OUTPUTS_NOTE", ""),
        db_setup_settings=tuple(getattr(module, "DB_SETUP_SETTINGS", ())),
        ref_dbs_needed=getattr(module, "ref_dbs_needed", None),
        ref_dbs_in_output_dir=getattr(module, "REF_DBS_IN_OUTPUT_DIR", False),
        runs_gtdbtk=getattr(module, "runs_gtdbtk", None),
        downloads_from=dict(getattr(module, "DOWNLOADS_FROM", {})),
    )


############################################################
################## building the settings ###################
############################################################

def _fail(message):
    report_message(message, initial_indent="    ", subsequent_indent="    ")
    notify_premature_exit()


def read_config_file(path, params):
    path = Path(path)
    if not path.is_file():
        _fail(f"The specified config file '{path}' was not found.")
    with open(path) as f:
        loaded = yaml.safe_load(f) or {}
    if not isinstance(loaded, dict):
        _fail(f"The config file '{path}' doesn't look like a table of settings.")
    known = {p.name for p in params}
    unknown = sorted(k for k in loaded if k not in known)
    if unknown:
        _fail(f"The config file '{path}' holds settings this workflow doesn't have: {', '.join(unknown)}. "
              "The detailed help (-H) lists all of them.")
    return loaded


def _check_and_convert(param, value, source, base_dir):
    """Returns the value in the form it is stored in the config; exits with a message if it can't be used."""
    if param.kind == "flag":
        if isinstance(value, bool):
            return value
        if isinstance(value, str) and value.strip().lower() in ("true", "false", "yes", "no"):
            return value.strip().lower() in ("true", "yes")
        _fail(f"The setting '{param.name}' ({source}) needs to be true or false, but is '{value}'.")

    if param.kind == "paths":
        # one or more existing files or directories
        given = value if isinstance(value, (list, tuple)) else [value]
        paths = []
        for item in given:
            path = Path(str(item)).expanduser()
            if not path.is_absolute():
                path = Path(base_dir) / path
            path = Path(os.path.abspath(path))
            if not path.exists():
                _fail(f"A path given for '{param.name}' ({source}) was not found: {path}")
            paths.append(str(path))
        return paths

    if param.kind in ("file", "dir", "path"):
        path = Path(str(value)).expanduser()
        if not path.is_absolute():
            path = Path(base_dir) / path
        path = Path(os.path.abspath(path))
        if param.kind == "file" and not path.is_file():
            _fail(f"The file given for '{param.name}' ({source}) was not found: {path}")
        if param.kind == "dir" and not path.is_dir():
            _fail(f"The directory given for '{param.name}' ({source}) was not found: {path}")
        return str(path)

    try:
        value = param.type(value)
    except (TypeError, ValueError):
        _fail(f"The setting '{param.name}' ({source}) needs to be of type {param.type.__name__}, but is '{value}'.")
    if param.choices and value not in param.choices:
        _fail(f"The setting '{param.name}' ({source}) needs to be one of: {', '.join(str(c) for c in param.choices)}.")
    return value


def resolve_settings(params, cli_values, config_file=None, only_require=None):
    """
    Later sources win: defaults, then the config file, then flags given on the command line.
    cli_values: {name: value} holding only what was given on the command line.
    only_require: when given, the only required settings that still need to be given (for --setup-dbs-only).
    Relative paths are taken relative to the config file's directory (config file) or the current directory (flags).
    """
    from_file = read_config_file(config_file, params) if config_file else {}
    config_dir = Path(config_file).resolve().parent if config_file else Path.cwd()

    settings = {}
    for p in params:
        if p.name in cli_values:
            settings[p.name] = _check_and_convert(p, cli_values[p.name], "command line", Path.cwd())
        elif p.name in from_file and from_file[p.name] is not None:
            settings[p.name] = _check_and_convert(p, from_file[p.name], "config file", config_dir)
        elif p.default is not None:
            settings[p.name] = _check_and_convert(p, p.default, "default", Path.cwd())
        elif p.required and (only_require is None or p.name in only_require):
            flags = f"{p.short}/{p.flag}" if p.short else p.flag
            _fail(f"A required setting was not given: {flags} (or '{p.name}' in a config file).")
        else:
            settings[p.name] = None
    return settings


def write_config(workflow, settings, path, command=None):
    """Writes a commented config file holding every setting of the workflow."""
    lines = [
        f"# Config for the bit {workflow.name} workflow (bit v{version('bit')})",
        "# https://github.com/AstrobioMike/bit",
    ]
    if command:
        lines += ["#", "# Written by:"] + [f"#   {l}" for l in textwrap.wrap(command, 96, break_long_words=False, break_on_hyphens=False)]
    lines += [
        "#",
        f"# Use it with: {invocation(workflow.name)} --config-file <this file>",
        "# (settings given as flags on the command line take priority over the ones here)",
        "",
    ]
    groups = []
    for p in workflow.params:
        if p.group not in groups:
            groups.append(p.group)
    for group in groups:
        lines += ["#" * 80, f"## {group}", "#" * 80, ""]
        for p in (p for p in workflow.params if p.group == group):
            lines += [f"# {l}" for l in textwrap.wrap(" ".join(p.help.split()), 96)]
            lines.append(yaml.safe_dump({p.name: settings[p.name]}, default_flow_style=False, width=1000).rstrip("\n"))
            lines.append("")
    Path(path).write_text("\n".join(lines))


############################################################
############### saving a standalone workflow ###############
############################################################

def workflow_files(workflow, workflow_dir=None):
    """
    (source, relative destination) for every file of the workflow. In the package the shared files
    sit in bit/smk/; in a saved workflow they sit next to the Snakefile, at the same relative paths.
    """
    workflow_dir = Path(workflow_dir) if workflow_dir else workflow.directory
    saved = workflow_dir != workflow.directory
    shared_root = workflow_dir if saved else Path(get_package_path("smk"))
    pairs = [(workflow_dir / "Snakefile", Path("Snakefile"))]
    scripts_dir = workflow_dir / "scripts"
    if scripts_dir.is_dir():
        shared_names = {Path(f).name for f in workflow.shared_files if f.startswith("scripts/")}
        pairs += [(f, Path("scripts") / f.name) for f in sorted(scripts_dir.iterdir())
                  if f.is_file() and f.name not in shared_names]
    pairs += [(shared_root / f, Path(f)) for f in workflow.shared_files]
    return pairs


def get_workflow(workflow, settings, dest, command=None):
    """ Writes the copy into `dest`, which has to be new or empty: nothing already there is ever written over or removed. """
    dest = Path(dest)
    if dest.exists() and not dest.is_dir():
        _fail(f"The path given to --get-workflow is an existing file: {dest}. Please name a different directory.")
    if dest.is_dir() and any(dest.iterdir()):
        _fail(f"The directory given to --get-workflow already exists and isn't empty: {dest}. "
              "Please name a different directory, or remove that one yourself first.")
    for src, rel in workflow_files(workflow):
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, target)
    write_config(workflow, settings, dest / "config.yaml", command)
    return dest


############################################################
######################### running ##########################
############################################################

CLUSTER_ARGS = ("--executor", "-e", "--profile", "--workflow-profile", "--cluster", "--cluster-generic-submit-cmd", "--slurm")


def runs_on_cluster(snakemake_args):
    """ whether the snakemake arguments passed through hand the jobs to a cluster scheduler """
    return any(str(a).split("=")[0] in CLUSTER_ARGS for a in (snakemake_args or []))


def machine_memory_mb():
    """ total physical memory in MB, or None where it can't be read """
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 1_000_000
    except (ValueError, OSError, AttributeError):
        return None


def dry_run_report(dry_run_cmd):
    """ what a dry run of the workflow prints; None when the dry run itself fails (the real run then reports why) """
    try:
        process = subprocess.run(dry_run_cmd, capture_output=True, text=True)
    except OSError:
        return None
    if process.returncode != 0:
        return None
    return process.stdout + process.stderr


def planned_jobs_memory(report):
    """
    The steps a dry run of the workflow lists (its dry_run_report), with the memory (mem_mb) each asks for:
    {rule: mem_mb or None}. None when the dry run failed.
    """
    if report is None:
        return None

    steps, rule = {}, None
    for line in report.splitlines():
        found = re.match(r"(?:local)?(?:rule|checkpoint) (\S+):$", line.strip())
        if found:
            rule = found.group(1)
            steps.setdefault(rule, None)
            continue
        found = re.search(r"\bmem_mb=(\d+)", line)
        if rule and found and line.strip().startswith("resources:"):
            steps[rule] = max(steps[rule] or 0, int(found.group(1)))
    return steps


def check_steps_can_run(report, memory_mb, memory_given):
    """ stops before anything runs if a step planned in the dry run (its dry_run_report) can't run on this machine """

    steps = planned_jobs_memory(report)
    if not steps:
        return

    if memory_mb:
        where = f"the {memory_mb / 1000:g} GB given with --memory" if memory_given else \
                f"this machine has ({memory_mb / 1000:g} GB)"
        too_big = {r: mb for r, mb in steps.items() if mb and mb > memory_mb}
        if too_big:
            listed = ", ".join(f"'{r}' ({mb / 1000:g} GB)" for r, mb in sorted(too_big.items()))
            _fail(f"Some steps of this run need more memory than {where}: {listed}. Please run it on a machine with more "
                  "memory. If you know a step needs less than its setting says, lower that step's memory setting "
                  "(listed under 'Run Parameters' in the detailed help, -H).")


def run_workflow(workflow, settings, args, command):
    """
    The run directory is the output directory: the full config is written there as config.yaml
    (the Snakefile reads it), next to the outputs and a record of the command.
    """
    workflow_dir = Path(args.workflow_dir).resolve() if args.workflow_dir else workflow.directory
    snakefile = workflow_dir / "Snakefile"
    if not snakefile.is_file():
        _fail(f"No Snakefile was found in the workflow directory: {workflow_dir}")

    output_dir = Path(args.output_dir).resolve()

    # the reference databases: where they are (recorded in the run's config file), and which this run needs
    ref_dbs, ref_root, refs_in_output_dir = [], None, False
    if workflow.ref_dbs_needed:
        from bit.modules.ref_dbs import RefDbsDirNotSet, ref_dbs_root
        ref_dbs = workflow.ref_dbs_needed(settings)
        try:
            ref_root = ref_dbs_root(settings.get("ref_dbs_dir"))
        except RefDbsDirNotSet as e:
            if not workflow.ref_dbs_in_output_dir:
                _fail(str(e))
            ref_root = output_dir / "ref-dbs"
            refs_in_output_dir = True
        settings = {**settings, "ref_dbs_dir": str(ref_root)}

    # where GTDB-Tk can't run (macOS), the genomes' taxonomy comes from bit-tax-db instead: said up front
    if workflow.runs_gtdbtk and workflow.runs_gtdbtk(settings):
        from bit.modules.ref_dbs import GTDBTK_RUNS_HERE, GTDBTK_NOT_HERE_NOTE
        if not GTDBTK_RUNS_HERE:
            report_message("Note: " + GTDBTK_NOT_HERE_NOTE, color="yellow", initial_indent="    ", subsequent_indent="    ",
                           trailing_newline=True)

    # the reference databases not set up yet are downloaded before the run starts, so whether that can be done is
    # checked up front (bit.modules.connectivity), along with whatever else the run downloads
    missing_dbs = []
    if ref_dbs:
        from bit.modules.ref_dbs import missing, download_urls
        missing_dbs = missing(ref_root, ref_dbs)

    def require_internet_for(urls, what, env_note=None):
        from bit.modules.connectivity import require_internet
        require_internet(urls, what=what, missing_dbs=missing_dbs, ref_dbs_dir=ref_root,
                         ref_dbs_in_output_dir=refs_in_output_dir, env_note=env_note)

    if args.setup_dbs_only:
        from bit.modules.ref_dbs import set_up_missing
        if missing_dbs:
            require_internet_for(download_urls(missing_dbs), "the db")
        set_up_missing(ref_root, ref_dbs, check_internet=False)
        report_message("REFERENCE DATABASES SET UP!", color="green")
        print(f"    They are in: {color_text(str(ref_root) + '/', 'green')}\n")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    write_config(workflow, settings, output_dir / "config.yaml", command)

    threads = settings.get("threads") or 1

    # on one machine, snakemake is told how much memory there is, so it only starts jobs that fit together, and
    # rules whose memory is a ceiling (megahit) are kept within it; through a cluster executor or profile each
    # job's mem_mb is its request to the scheduler instead
    memory_mb = None
    if not runs_on_cluster(args.snakemake_args):
        memory_mb = int(args.memory * 1000) if args.memory else machine_memory_mb()

    # bit sets up the reference databases itself (below), so the workflow doesn't check for them
    config_pairs = (["ref_dbs_set_up_by_bit=True"] if ref_dbs else []) + \
                   ([f"memory_total_mb={memory_mb}"] if memory_mb else [])
    cmd_args = dict(
        snakefile=snakefile,
        directory=output_dir,
        cores=args.jobs * threads,
        use_conda=True,
        conda_prefix=args.conda_prefix,
        rerun_incomplete=args.rerun_incomplete,
        extra_args=args.snakemake_args,
        config_pairs=config_pairs or None,
    )

    # a dry run first, for what snakemake itself plans: each step's memory, and the conda environments it still
    # has to create
    planned = dry_run_report(build_snakemake_cmd(**cmd_args, dry_run=True))

    if memory_mb:
        check_steps_can_run(planned, memory_mb, bool(args.memory))

    # what this run would download, and from where: the missing reference databases, the conda environments of its
    # steps that aren't there yet, and what the steps it still has to do download themselves (where they run on
    # this machine)
    envs_to_create = planned_conda_envs(planned)
    step_urls = []
    if workflow.downloads_from and not runs_on_cluster(args.snakemake_args):
        to_do = planned_jobs_memory(planned) or {}
        step_urls = [url for rule, rule_urls in workflow.downloads_from.items() if rule in to_do for url in rule_urls]
    urls = (download_urls(missing_dbs) if missing_dbs else []) + conda_env_urls(envs_to_create) + step_urls
    if urls and not args.dry_run:
        env_note = None
        if envs_to_create:
            from bit.modules.connectivity import shared_conda_prefix_note
            env_note = shared_conda_prefix_note(
                [Path(env).stem for env in envs_to_create],
                # (`bit run-wf` passes what follows `--` to snakemake)
                "adding `-- --conda-create-envs-only` to the same command there creates just them"
                if workflow.name in WORKFLOWS else None)
        require_internet_for(urls, "the db" if missing_dbs else ("the data" if step_urls else None), env_note)

    # any reference database missing is set up now, before snakemake starts (so the progress bars have the terminal)
    if ref_dbs and not args.dry_run:
        from bit.modules.ref_dbs import set_up_missing
        set_up_missing(ref_root, ref_dbs, check_internet=False)

    cmd = build_snakemake_cmd(**cmd_args, dry_run=args.dry_run,
                              resources=[f"mem_mb={memory_mb}"] if memory_mb else None)
    run_snakemake(cmd)

    if not args.dry_run:
        with open(output_dir / "command-execution-info.txt", "w") as f:
            f.write(f"bit version: v{version('bit')}\n\n")
            f.write(f"Rendered command:\n{command}\n\n")
            f.write(f"Snakemake command:\n{shlex.join(cmd)}\n\n")

    border = "-" * 80
    print(f"\n{border}")
    if args.dry_run:
        report_message("DRY-RUN COMPLETE!", color="green")
        print("    No workflow steps were run, but you can see what would have been done above.")
        if missing_dbs:
            print(f"    These reference databases would be set up first in {ref_root}/: {', '.join(missing_dbs)}")
        if envs_to_create:
            print(f"    These conda environments would be created first: {', '.join(Path(e).stem for e in envs_to_create)}")
        if urls:
            from bit.modules.connectivity import note_for_a_dry_run
            note_for_a_dry_run(urls)
    else:
        report_message("DONE!", color="green")
        print(f"    All outputs are in: {color_text(str(output_dir) + '/', 'green')}")
        if workflow.outputs_note:
            print(f"    {workflow.outputs_note}")
    print(f"\n{border}\n")
