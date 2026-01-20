import os
import subprocess
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from Bio import SeqIO #type: ignore
import numpy as np
from colorama import Fore, init
from bit.modules.general import check_files_are_found
from bit.modules.get_mapped_reads_pid import get_mapped_reads_pids


init(autoreset=True)


def get_cov_stats(args):

    preflight_checks(args)
    refs = parse_refs(args.reference_fastas)

    ref_mean_pids, ref_mapped_counts = compute_ref_pid_stats(args.bam, refs)

    if args.bed:
        refs = parse_bed_file(refs, args.bed)
    else:
        bed = run_mosdepth(args.bam, args.output_prefix)
        refs = parse_bed_file(refs, bed)

    generate_output(refs, args.output_prefix, ref_mean_pids, ref_mapped_counts)


def preflight_checks(args):

    paths_list = list(args.reference_fastas)
    if args.bed:
        paths_list.append(args.bed)
    if args.bam:
        paths_list.append(args.bam)
    check_files_are_found(paths_list)


def parse_refs(reference_fastas):
    refs = []
    for fasta in reference_fastas:
        ref = RefData(fasta)
        ref.load_fasta()
        refs.append(ref)

    return refs


def compute_ref_pid_stats(bam_path, refs):

    if not bam_path:
        return None, None

    try:
        ref_read_pids, _all_pids = get_mapped_reads_pids(bam_path)
    except Exception:
        return None, None

    ref_mean_pids = {}
    ref_mapped_counts = {}

    for ref in refs:
        pids = []
        mapped_count = 0
        for header in ref.headers:
            if header in ref_read_pids:
                header_list = ref_read_pids[header]
                mapped_count += len(header_list)
                pids.extend(pid for (_rid, pid) in header_list)

        ref_mapped_counts[ref.path] = mapped_count
        ref_mean_pids[ref.path] = round(float(np.mean(pids)), 2) if pids else None

    return ref_mean_pids, ref_mapped_counts


@dataclass
class RefData:
    path: str
    headers: set = field(default_factory=set)
    total_length: int = 0
    total_coverage_count: int = 0
    total_bases_detected_at_all: int = 0
    total_bases_detected_at_10x: int = 0

    def load_fasta(self):
        with open(self.path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.headers.add(record.id)
                self.total_length += len(record.seq)

    def update_from_bed_line(self, header: str, start: int, end: int, num_reads: int):
        if header in self.headers:
            cur_range_covered = end - start
            self.total_coverage_count += num_reads * cur_range_covered
            if num_reads > 0:
                self.total_bases_detected_at_all += cur_range_covered
            if num_reads >= 10:
                self.total_bases_detected_at_10x += cur_range_covered

    def compute_metrics(self):
        detection = round(self.total_bases_detected_at_all / self.total_length, 4)
        detection_at_10x = round(self.total_bases_detected_at_10x / self.total_length, 4)
        mean_coverage = round(self.total_coverage_count / self.total_length, 4)
        return detection, detection_at_10x, mean_coverage


def parse_bed_file(refs, bed_file):

    print(f"\n  {Fore.YELLOW}Parsing coverage info...")

    with gzip.open(bed_file, "rt") as f:
        for line in f:
            header, start, end, num_reads = line.strip().split("\t")
            for ref in refs:
                ref.update_from_bed_line(header, int(start), int(end), int(num_reads))

    return(refs)


def run_mosdepth(bam_file, output_prefix):
    check_bam_file_is_indexed(bam_file)
    mosdepth_output_dir = f"{output_prefix}-mosdepth-files"
    mosdepth_prefix = str(Path(mosdepth_output_dir) / Path(bam_file).stem)
    os.makedirs(mosdepth_output_dir, exist_ok=True)
    cmd = f"mosdepth {mosdepth_prefix} {bam_file}"
    print(f"\n  {Fore.YELLOW}Running mosdepth to generate the required per-base coverage file...")
    subprocess.run(cmd, shell=True)
    bed_file = f"{mosdepth_prefix}.per-base.bed.gz"
    print(f"\n    Mosdepth outputs can be found in: {Fore.YELLOW}{mosdepth_output_dir}/\n")

    return bed_file


def check_bam_file_is_indexed(bam_file):
    if not Path(bam_file + ".bai").is_file():
        cmd = f"samtools index {bam_file}"
        subprocess.run(cmd, shell=True)


def generate_output(refs, output_prefix, ref_mean_pids=None, ref_mapped_counts=None):

    primary_out_path = f"{output_prefix}.tsv"

    with open(primary_out_path, "w") as f:

        if ref_mean_pids is not None:
            f.write("ref\tdetection\tdetection_at_10x\tmean_coverage\tmean_pid\tnum_mapped_reads\n")
        else:
            f.write("ref\tdetection\tdetection_at_10x\tmean_coverage\n")

        for ref in refs:
            detection, detection_at_10x, mean_coverage = ref.compute_metrics()

            if ref_mean_pids is not None:
                mean_pid = ref_mean_pids.get(ref.path)
                mapped_reads = ref_mapped_counts.get(ref.path, 0) if ref_mapped_counts is not None else 0
                mapped_reads_str = f"{mapped_reads:,}"
                mean_pid_str = "NA" if mean_pid is None else f"{mean_pid}"
                f.write(f"{ref.path}\t{detection}\t{detection_at_10x}\t{mean_coverage}\t{mean_pid_str}\t{mapped_reads_str}\n")
            else:
                f.write(f"{ref.path}\t{detection}\t{detection_at_10x}\t{mean_coverage}\n")

    print(f"\n    Coverage stats written to: {Fore.YELLOW}{primary_out_path}\n")
