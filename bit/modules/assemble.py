import os
import json
from dataclasses import dataclass, field
from pathlib import Path
from subprocess import run
import numpy as np
from bit.modules.input_parsing import (get_input_reads_dict_from_dir,
                                       get_input_reads_dict_from_paths)
from bit.modules.general import (report_message,
                       report_failure,
                       get_package_path,
                       color_text,
                       log_command_run)
import sys
import platform


def run_assembly(args, full_cmd_executed):
    reads_dict = gen_reads_dict(args)
    config = RunConfig.from_args(args)
    config = check_conda_env(config)
    run_snakemake(reads_dict, config)
    log_command_run(full_cmd_executed, config.output_dir)
    report_finished(args)

def gen_reads_dict(args):
    if args.reads_dir:
        reads_dict = get_input_reads_dict_from_dir(args.reads_dir)
    else:
        reads_dict = get_input_reads_dict_from_paths(args.read_1, args.read_2)

    return reads_dict

@dataclass
class RunConfig:
    output_dir: Path = field(init=None)
    base_output_dir: str = field(init=None)
    log_files_dir: Path = field(init=None)
    base_log_files_dir: str = field(init=None)
    reads_dir: Path = field(init=None)
    assembler: str = field(init=None)
    threads: int = field(init=None)
    run_fastp: bool = field(init=False)
    run_bbnorm: bool = field(init=False)
    memory: str = field(init=None)
    min_contig_len: int = field(init=None)
    isolate: bool = field(init=False)
    num_cores: int = field(init=None)
    rerun_incomplete: bool = field(init=False)
    dry_run: bool = field(init=False)
    conda_env: str = field(init=False)

    @classmethod
    def from_args(cls, args):
        config = cls()
        config.populate_run_config(args)
        return config

    def populate_run_config(self, args):
        self.output_dir = Path(args.output_dir).resolve()
        self.base_output_dir = self.output_dir.name
        self.log_files_dir = self.output_dir / "log-files"
        self.base_log_files_dir = f"{self.base_output_dir}/log-files"
        self.reads_dir = Path(args.reads_dir).absolute() if args.reads_dir else None
        self.assembler = args.assembler
        self.threads = args.threads
        self.run_fastp = args.run_fastp
        self.run_bbnorm = args.run_bbnorm
        self.memory = args.memory
        self.min_contig_len = args.min_contig_len
        self.isolate = args.isolate
        self.num_cores = args.jobs
        self.rerun_incomplete = args.rerun_incomplete
        self.dry_run = args.dry_run

    @property
    def key_value_pairs(self):
        return [f"{key}={str(value)}" for key, value in vars(self).items()]


def check_conda_env(config):

    conda_env_area = os.path.expandvars("${CONDA_PREFIX}/envs")
    conda_env = f"bit-assemble"
    conda_env_path = f"{conda_env_area}/{conda_env}"
    if sys.platform.startswith("linux"):
        yaml_path = str(get_package_path("smk/envs/assemble.yaml"))
    elif sys.platform.startswith("darwin"):
        yaml_path = str(get_package_path("smk/envs/assemble-osx-64.yaml"))
    else:
        message = f"Unsupported platform: {platform.system()}"
        report_message(message, "red")
        report_failure("Please use linux-64 or osx-64.", "yellow")

    if not os.path.exists(conda_env_path):
        try:
            report_message(f"Creating conda environment for bit-assemble at {conda_env_path}...")
            print("")
            run(["conda", "env", "create", "--file", yaml_path, "-p", conda_env_path])
            report_message("Conda environment created successfully!\nMoving on :)\n", "green")
            print("")
        except Exception as e:
            message = f"Failed to create conda environment for bit-assemble at {conda_env_path}."
            report_message(message, "red")
            report_message("Please check the error below:")
            print("\n" + str(e) + "\n")
            exit(1)

    config.conda_env = conda_env
    return config


def run_snakemake(reads_dict, config):
    reads_json = json.dumps(reads_dict)
    conda_env_area = os.path.expandvars("${CONDA_PREFIX}/envs")

    cmd = [
        "snakemake",
        "--snakefile", str(get_package_path("smk/assemble.smk")),
        "--cores", str(config.num_cores),
        "--printshellcmds",
        "--use-conda",
        "--conda-prefix", str(conda_env_area),
        "--directory", config.output_dir,
        "--config", f'reads_json={reads_json}',
        *config.key_value_pairs,
    ]

    if config.dry_run:
        cmd.append("--dry-run")
    if config.rerun_incomplete:
        cmd.append("--rerun-incomplete")

    process = run(cmd)
    if process.returncode != 0:
        message = "Snakemake failed. Hopefully its output above can help you spot why."
        report_failure(message)


def report_finished(args):
    border = "-" * 80
    print(f"\n{border}")
    if not args.dry_run:
        report_message("DONE!", color = "green")
        out_dir = str(args.output_dir) if str(args.output_dir).endswith("/") else f"{args.output_dir}/"
        out_file = Path(args.output_dir) / "assembly-summaries.tsv"
        print(f"    All outputs are in: {color_text(out_dir, 'green')}")
        print(f"    Summary table written to: {color_text(out_file, 'green')}")
    else:
        report_message("DRY-RUN COMPLETE!", color = "green")
        print("    No files were written, but you can see what would have been done above.")
    print(f"\n{border}\n")
