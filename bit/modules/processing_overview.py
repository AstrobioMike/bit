"""
The metagenomics workflow's per-sample processing overview: whether each step produced anything, and how much
of each sample's reads its assembly, bins and MAGs account for (the share of its reads that mapped to them).

Reads are counted as single reads (each mate of a pair on its own), as bowtie2's overall alignment rate does.
The total comes from bowtie2's summary, the mapped reads per contig from `samtools idxstats` of the sample's
mapping to its own assembly.
"""

import glob
import os


def total_reads(mapping_info_path, paired):
    """
    reads given to bowtie2 (mates counted separately), from its summary's "N reads; of these:" line (any warnings bowtie2
    wrote about particular reads come before it, in the same file)
    """

    with open(mapping_info_path) as f:
        for line in f:
            fields = line.split()
            if line.rstrip().endswith("reads; of these:") and fields[0].isdigit():
                return int(fields[0]) * (2 if paired else 1)
    return 0


def mapped_reads_per_contig(idxstats_path):
    """ {contig: reads mapped to it} from `samtools idxstats` output (name, length, mapped, unmapped) """

    mapped = {}
    with open(idxstats_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] != "*":
                mapped[fields[0]] = int(fields[2])
    return mapped


def contigs_in(fasta_paths):
    """ the sequence names (up to the first space) in the given fasta files """

    names = set()
    for path in fasta_paths:
        with open(path) as f:
            names.update(line[1:].split()[0] for line in f if line.startswith(">") and len(line) > 2)
    return names


def percent(part, whole):
    return f"{100 * part / whole:.2f}" if whole else "NA"


def write_processing_overview(sample_IDs, out_path, assembly_of, mapping_dir, paired,
                              genes_dir=None, bins_dir=None, MAGs_dir=None):
    """
    One row per sample.
    assembly_of: {sample: (name, fasta path) of the assembly its reads were mapped to} (its own, or a co-assembly)
    genes_dir is given when genes were called, bins_dir and MAGs_dir when binning was run (the bins and MAGs of an
    assembly are named <assembly name>-bin.<n> and <assembly name>-MAG-<n>); their columns are left out otherwise.
    """

    genes = genes_dir is not None
    binning = bins_dir is not None
    header = ["Sample_ID", "assembly_produced"] + (["gene_calls_identified"] if genes else []) + ["read_mapping_successful"]
    if binning:
        header += ["bins_recovered", "MAGs_recovered"]
    header += ["reads", "percent_of_reads_mapped_to_assembly"]
    if binning:
        header += ["percent_of_reads_mapped_to_bins", "percent_of_reads_mapped_to_MAGs"]

    rows = []
    for sample in sample_IDs:
        asm_name, asm_path = assembly_of[sample]
        bam = os.path.join(mapping_dir, f"{sample}.bam")
        row = [sample, "Yes" if os.path.getsize(asm_path) else "No"]
        if genes:
            row.append("Yes" if os.path.getsize(os.path.join(genes_dir, f"{sample}-genes.faa")) else "No")
        row.append("Yes" if os.path.getsize(bam) else "No")

        bins = sorted(glob.glob(os.path.join(bins_dir, f"{glob.escape(asm_name)}-bin.[0-9]*.fasta"))) if binning else []
        MAGs = sorted(glob.glob(os.path.join(MAGs_dir, f"{glob.escape(asm_name)}-MAG-[0-9]*.fasta"))) if binning else []
        if binning:
            row += [str(len(bins)), str(len(MAGs))]

        reads = total_reads(os.path.join(mapping_dir, f"{sample}-mapping-info.txt"), paired)
        mapped = mapped_reads_per_contig(os.path.join(mapping_dir, f"{sample}-mapped-reads-per-contig.tsv"))
        row += [str(reads), percent(sum(mapped.values()), reads)]
        if binning:
            for genomes in (bins, MAGs):
                row.append(percent(sum(mapped.get(c, 0) for c in contigs_in(genomes)), reads))
        rows.append(row)

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
