import os
import subprocess
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from Bio import SeqIO #type: ignore
from collections import defaultdict
from colorama import Fore, init #type: ignore
from tqdm import tqdm #type: ignore
from bit.modules.general import check_files_are_found, report_message
from bit.modules.get_mapped_reads_pid import get_mapped_reads_pids


init(autoreset=True)


def get_cov_stats(args):

    preflight_checks(args)

    refs, contigs, header_to_ref_dict = parse_refs(args.reference_fastas, args.skip_per_contig)

    (ref_mean_pids, ref_mapped_counts,
     contig_mean_pids, contig_mapped_counts) = compute_pid_stats(args.bam, refs, contigs,
                                                                 args.include_non_primary,
                                                                 args.skip_per_contig)

    bed = args.bed if args.bed else run_mosdepth(args.bam, args.output_prefix)

    refs = parse_bed_file(refs, bed, contigs, header_to_ref_dict)

    generate_output(refs, args.output_prefix, ref_mean_pids, ref_mapped_counts,
                    contigs, contig_mean_pids, contig_mapped_counts, args.skip_per_contig)


def preflight_checks(args):

    paths_list = list(args.reference_fastas)
    if args.bed:
        paths_list.append(args.bed)
    if args.bam:
        paths_list.append(args.bam)
    check_files_are_found(paths_list)


@dataclass
class CoverageStats:
    length: int = 0
    total_coverage_count: int = 0
    total_bases_detected_at_all: int = 0
    total_bases_detected_at_10x: int = 0
    depth_hist: dict = field(default_factory=lambda: defaultdict(int))

    def update_from_bed_line(self, span, num_reads):
        self.total_coverage_count += num_reads * span
        if num_reads > 0:
            self.total_bases_detected_at_all += span
        if num_reads >= 10:
            self.total_bases_detected_at_10x += span

        self.depth_hist[num_reads] += span

    def compute_metrics(self):
        if self.length == 0:
            return 0.0, 0.0, 0.0
        detection = round(self.total_bases_detected_at_all / self.length, 2)
        detection_at_10x = round(self.total_bases_detected_at_10x / self.length, 2)
        mean_coverage = round(self.total_coverage_count / self.length, 2)
        median_coverage = self._median_from_hist()

        return detection, detection_at_10x, mean_coverage, median_coverage

    def _median_from_hist(self):
        if self.length == 0:
            return 0.0

        # getting two (possible) middle indices
        i1 = (self.length - 1) // 2
        i2 = self.length // 2

        # handling potential missing bases as 0s
        total_counted = sum(self.depth_hist.values())
        needed_extra_zeros = self.length - total_counted

        cumulative_count = 0
        v1 = v2 = None

        depths = sorted(self.depth_hist)
        if 0 not in self.depth_hist:
            depths = [0] + depths

        for depth in depths:
            count = self.depth_hist[depth]

            if depth == 0:
                count += needed_extra_zeros

            next_cumulative_count = cumulative_count + count

            if v1 is None and i1 < next_cumulative_count:
                v1 = depth
            if v2 is None and i2 < next_cumulative_count:
                v2 = depth
                break

            cumulative_count = next_cumulative_count

        return round((float(v1) + float(v2)) / 2.0, 2)


@dataclass
class ContigData:
    path: str
    header: str
    stats: CoverageStats = field(default_factory=CoverageStats)


@dataclass
class RefData:
    path: str
    headers: set = field(default_factory=set)
    contig_order: list = field(default_factory=list)
    stats: CoverageStats = field(default_factory=CoverageStats)


def parse_refs(reference_fastas, skip_per_contig=False):

    refs = []
    contigs_dict = {}
    header_to_ref_dict = {}

    for fasta in reference_fastas:
        ref = RefData(fasta)
        with open(fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                ref.headers.add(record.id)
                ref.contig_order.append(record.id)
                rec_len = len(record.seq)
                ref.stats.length += rec_len
                if not skip_per_contig:
                    contigs_dict[record.id] = ContigData(path=fasta, header=record.id, stats=CoverageStats(length=rec_len))
                header_to_ref_dict[record.id] = ref
        refs.append(ref)

    return refs, contigs_dict, header_to_ref_dict


def compute_pid_stats(bam_path, refs, contigs_dict, include_non_primary=False, skip_per_contig=False):

    if not bam_path:
        return None, None, None, None

    try:
        ref_read_pids, _all_pids = get_mapped_reads_pids(bam_path, include_non_primary)
    except Exception:
        return None, None, None, None

    contig_mean_pids = {}
    contig_mapped_counts = {}

    if not skip_per_contig:
        for contig_name in contigs_dict.keys():
            if contig_name in ref_read_pids:
                read_data = ref_read_pids[contig_name]
                contig_sum = sum(pid for (_, pid) in read_data)
                contig_count = len(read_data)

                contig_mapped_counts[contig_name] = contig_count
                contig_mean_pids[contig_name] = round(contig_sum / contig_count, 2) if contig_count > 0 else None
            else:
                contig_mapped_counts[contig_name] = 0
                contig_mean_pids[contig_name] = None

    ref_mean_pids = {}
    ref_mapped_counts = {}

    for ref in refs:
        ref_sum = 0.0
        ref_count = 0
        for contig_name in ref.headers:
            if contig_name in ref_read_pids:
                read_data = ref_read_pids[contig_name]
                ref_sum += sum(pid for (_, pid) in read_data)
                ref_count += len(read_data)

        ref_mapped_counts[ref.path] = ref_count
        ref_mean_pids[ref.path] = round(ref_sum / ref_count, 2) if ref_count > 0 else None

    return ref_mean_pids, ref_mapped_counts, contig_mean_pids, contig_mapped_counts


def parse_bed_file(refs, bed_file, contigs, header_to_ref_dict):

    print(f"\n  {Fore.YELLOW}Parsing coverage info...")

    with gzip.open(bed_file, "rt") as f:
        for line in tqdm(f, unit=" lines", desc="    Processing bed file", leave=False):
            header, start, end, num_reads = line.strip().split("\t")
            span = int(end) - int(start)
            n_reads = int(num_reads)

            if header in contigs:
                contigs[header].stats.update_from_bed_line(span, n_reads)

            if header in header_to_ref_dict:
                header_to_ref_dict[header].stats.update_from_bed_line(span, n_reads)

    return refs


def run_mosdepth(bam_file, output_prefix):
    check_bam_file_is_indexed(bam_file)
    mosdepth_output_dir = f"{output_prefix}-mosdepth-files"
    mosdepth_prefix = str(Path(mosdepth_output_dir) / Path(bam_file).stem)
    os.makedirs(mosdepth_output_dir, exist_ok=True)
    cmd = f"mosdepth -x {mosdepth_prefix} {bam_file}"
    print(f"\n  {Fore.YELLOW}Running mosdepth to generate the required per-base coverage file...")
    subprocess.run(cmd, shell=True)
    bed_file = f"{mosdepth_prefix}.per-base.bed.gz"
    print(f"\n    Mosdepth outputs can be found in: {Fore.YELLOW}{mosdepth_output_dir}/\n")

    return bed_file


def check_bam_file_is_indexed(bam_file):
    if not Path(bam_file + ".bai").is_file():
        cmd = f"samtools index {bam_file}"
        subprocess.run(cmd, shell=True)
        message = """
                  We indexed the BAM for you. Why? Because it's common courtesy.
                  All the programs that don't do this for us when it's needed are just big jerks!
                  """
        report_message(message, color="orange", initial_indent="    ", subsequent_indent="    ")


def generate_output(refs, output_prefix, ref_mean_pids=None, ref_mapped_counts=None,
                    contigs=None, contig_mean_pids=None, contig_mapped_counts=None,
                    skip_per_contig=False):

    primary_out_path = f"{output_prefix}-per-ref.tsv"

    with open(primary_out_path, "w") as f:

        cols = ["ref", "detection", "detection_at_10x", "mean_coverage", "median_coverage"]
        if ref_mean_pids is not None:
            cols.extend(["mean_pid", "num_mapped_reads"])
        f.write("\t".join(cols) + "\n")

        for ref in refs:
            detection, detection_at_10x, mean_coverage, median_coverage = ref.stats.compute_metrics()

            if ref_mean_pids is not None:
                mean_pid = ref_mean_pids.get(ref.path)
                mapped_reads = ref_mapped_counts.get(ref.path, 0)
                mean_pid_str = "NA" if mean_pid is None else f"{mean_pid}"
                f.write(f"{ref.path}\t{detection}\t{detection_at_10x}\t{mean_coverage}\t{median_coverage}\t{mean_pid_str}\t{mapped_reads:,}\n")
            else:
                f.write(f"{ref.path}\t{detection}\t{detection_at_10x}\t{mean_coverage}\t{median_coverage}\n")

    print(f"\n    Ref-level coverage stats written to: {Fore.YELLOW}{primary_out_path}\n")

    if not skip_per_contig:

        contig_out_path = f"{output_prefix}-per-contig.tsv"

        with open(contig_out_path, "w") as cf:

            cols = ["ref", "contig", "length", "detection", "detection_at_10x", "mean_coverage", "median_coverage"]

            if contig_mean_pids is not None:
                cols.extend(["mean_pid", "num_mapped_reads"])
            cf.write("\t".join(cols) + "\n")

            written_contigs = set()
            for ref in sorted(refs, key=lambda r: r.path):
                for contig_name in ref.contig_order:
                    if contig_name not in contigs:
                        continue
                    contig_data = contigs[contig_name]
                    detection, detection_at_10x, mean_coverage, median_coverage = contig_data.stats.compute_metrics()

                    if contig_mean_pids is not None:
                        mean_pid = contig_mean_pids.get(contig_name)
                        mapped_reads = contig_mapped_counts.get(contig_name, 0)
                        mean_pid_str = "NA" if mean_pid is None else f"{mean_pid}"
                        cf.write(f"{contig_data.path}\t{contig_name}\t{contig_data.stats.length}\t{detection}\t{detection_at_10x}\t{mean_coverage}\t{median_coverage}\t{mean_pid_str}\t{mapped_reads:,}\n")
                    else:
                        cf.write(f"{contig_data.path}\t{contig_name}\t{contig_data.stats.length}\t{detection}\t{detection_at_10x}\t{mean_coverage}\t{median_coverage}\n")
                    written_contigs.add(contig_name)

        print(f"    Contig-level coverage stats written to: {Fore.YELLOW}{contig_out_path}\n")
