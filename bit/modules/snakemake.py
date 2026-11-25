"""
The one place where bit builds and runs snakemake commands (used by `bit assemble`,
`bit ez-screen reads` and `bit run-wf`), so anything specific to the snakemake
version lives here.
"""

import json
import os
import re
import subprocess
from subprocess import run
from bit.modules.general import report_failure

# where conda gets a channel given by name, unless its configuration says otherwise (channel_alias, custom_channels)
CONDA_CHANNEL_ALIAS = "https://conda.anaconda.org"
# where pip gets the packages of an environment file's `pip:` section, unless PIP_INDEX_URL says otherwise
PIP_INDEX = "https://pypi.org"


def default_conda_prefix():
    """Where per-rule conda environments are created (and re-used across runs)."""
    return os.path.expandvars("${CONDA_PREFIX}/envs")


def build_snakemake_cmd(snakefile, directory, cores, config_pairs=None, use_conda=False,
                        conda_prefix=None, dry_run=False, rerun_incomplete=False,
                        extra_args=None, resources=None):
    """
    Returns the snakemake command as a list.

    resources:    `name=value` strings passed with --resources (totals the jobs running at once must fit in)
    config_pairs: `key=value` strings passed with --config
    use_conda:    have snakemake create/activate the environments named in the rules
    extra_args:   anything else to hand to snakemake as given
    """
    # extra_args go first: a target given there (e.g. a rule name) would otherwise be read as another value of an
    # option that takes several (--resources, --config); bit's own options that follow all start with "--"
    cmd = [
        "snakemake", *[str(a) for a in extra_args or []],
        "--snakefile", str(snakefile),
        "--cores", str(cores),
        "--printshellcmds",
        "--directory", str(directory),
    ]

    if use_conda:
        cmd += ["--software-deployment-method", "conda",
                "--conda-prefix", str(conda_prefix or default_conda_prefix())]

    if dry_run:
        cmd.append("--dry-run")
    if rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if resources:
        cmd += ["--resources", *resources]

    # --config takes all following key=value items, so it goes last
    if config_pairs:
        cmd += ["--config", *config_pairs]

    return cmd


def planned_conda_envs(dry_run_report):
    """
    The environment files of the conda environments snakemake says (in a dry run's output) it would still have to
    create: those of the run's steps that aren't in the conda prefix yet, or are there incomplete.
    """

    found = re.findall(r"^(?:Incomplete )?Conda environment (.+) will be (?:re)?created\.\s*$", dry_run_report or "",
                       flags=re.MULTILINE)
    return list(dict.fromkeys(found))


def _conda_channel_settings():
    """ (channel_alias, custom_channels) of the conda configuration; the defaults where it can't be read """

    try:
        shown = subprocess.run(["conda", "config", "--show", "channel_alias", "custom_channels", "--json"],
                               capture_output=True, text=True, timeout=60)
        settings = json.loads(shown.stdout)
        alias = settings.get("channel_alias") or CONDA_CHANNEL_ALIAS
        custom = settings.get("custom_channels") or {}
        if isinstance(alias, str) and isinstance(custom, dict):
            return alias, custom
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    return CONDA_CHANNEL_ALIAS, {}


def conda_env_urls(env_files):
    """
    Where creating the conda environments of the given environment files downloads from: the channels each lists
    (by name: at conda's channel_alias, conda.anaconda.org unless configured otherwise), and PyPI for a `pip:` section.
    """

    import yaml # type: ignore

    channels, pip = [], False
    for env_file in env_files:
        try:
            with open(env_file) as f:
                env = yaml.safe_load(f) or {}
        except (OSError, yaml.YAMLError):
            continue
        channels += [str(c) for c in env.get("channels") or []]
        pip = pip or any(isinstance(d, dict) and "pip" in d for d in env.get("dependencies") or [])

    urls = []
    if channels:
        alias, custom = _conda_channel_settings()
        for channel in dict.fromkeys(channels):
            if "://" in channel:
                urls.append(channel)
            elif isinstance(custom.get(channel), str) and "://" in custom[channel]:
                urls.append(custom[channel])
            else:
                urls.append(alias)
    if pip:
        urls.append(os.environ.get("PIP_INDEX_URL") or PIP_INDEX)
    return list(dict.fromkeys(urls))


def run_snakemake(cmd):
    process = run(cmd)
    if process.returncode != 0:
        message = "Snakemake failed. Hopefully its output above can help you spot why."
        report_failure(message)
