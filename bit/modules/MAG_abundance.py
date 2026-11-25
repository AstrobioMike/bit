"""
The metagenomics workflow's MAG-by-sample tables: the MAGs of all samples are dereplicated (dRep), each sample's
reads are mapped to the dereplicated set together, and `bit cov-stats` of each sample's mapping gives per MAG its
detection, the reads it recruited, and its mean coverage.
"""

import csv
import os

from bit.modules.processing_overview import total_reads, percent


def read_cov_stats(per_ref_path):
    """ {MAG ID: (detection %, mean coverage, reads mapped)} from `bit cov-stats`' per-ref table (refs are MAG fastas) """

    stats = {}
    with open(per_ref_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            MAG = os.path.basename(row["ref"])
            MAG = MAG[:-len(".fasta")] if MAG.endswith(".fasta") else MAG
            stats[MAG] = (float(row["detection (%)"]), float(row["mean_coverage"]),
                          int(row["num_mapped_reads"].replace(",", "")))
    return stats


def write_MAG_tables(sample_IDs, cov_stats_paths, mapping_info_paths, paired, detection_out, percent_out, coverage_out):
    """
    Three MAG-by-sample tables (rows the dereplicated MAGs, columns the samples): detection (percent of the MAG's
    bases covered by at least one read), the percent of the sample's reads (mates counted separately) that mapped
    to it, and its mean coverage. Each says so if there were no MAGs.
    """

    per_sample = {s: read_cov_stats(p) for s, p in zip(sample_IDs, cov_stats_paths)}
    reads = {s: total_reads(p, paired) for s, p in zip(sample_IDs, mapping_info_paths)}
    MAGs = sorted({m for stats in per_sample.values() for m in stats})

    for out_path, value in ((detection_out, lambda s, d: f"{d[0]:.2f}"),
                            (percent_out, lambda s, d: percent(d[2], reads[s])),
                            (coverage_out, lambda s, d: f"{d[1]:.2f}")):
        with open(out_path, "w") as out:
            if not MAGs:
                out.write("There were no MAGs recovered.\n")
                continue
            out.write("\t".join(["MAG_ID", *sample_IDs]) + "\n")
            for MAG in MAGs:
                out.write("\t".join([MAG, *(value(s, per_sample[s][MAG]) for s in sample_IDs)]) + "\n")
