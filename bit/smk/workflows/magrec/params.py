from bit.modules.workflows import Param, RUN_GROUP
from bit.modules.ref_dbs import REF_DBS_DIR_HELP
from bit.smk.workflows.classify_seqs.params import labelling_params

DESCRIPTION = """
    This workflow recovers metagenome-assembled genomes (MAGs) from reads and assemblies made
    already: each sample's reads are mapped to its assembly (bowtie2), or every sample's to one
    co-assembly, the assemblies are binned (metabat2), the bins estimated for quality (checkm2),
    and those passing the cutoffs become MAGs, which get taxonomy (GTDB-Tk). The MAGs of all
    assemblies are dereplicated (dRep), and each sample's reads mapped to them for tables of each
    MAG's detection, share of the reads, and coverage per sample. The reference databases are set
    up the first time they are needed.
    """

EXAMPLE = "bit run-wf magrec -r reads-dir/ -a assemblies-dir/"

OUTPUTS_NOTE = "The MAGs and their tables are in the 'MAGs/' subdirectory, the bins in 'bins/'."
DB_SETUP_SETTINGS = ()


def ref_dbs_needed(settings):
    """ the reference databases (bit.modules.ref_dbs) a run uses """
    return binning_ref_dbs()


def binning_ref_dbs():
    """ those of binning and MAG recovery (rules/magrec.smk): checkm2, and GTDB-Tk or (where it doesn't run) bit-tax-db """
    from bit.modules.ref_dbs import GTDBTK_RUNS_HERE
    return ["checkm2-db", "gtdbtk-db" if GTDBTK_RUNS_HERE else "tax-db"]



def runs_gtdbtk(settings):
    """ whether a run with these settings would use GTDB-Tk (where it runs) """
    return True

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "rules/reads.smk",
    "rules/ref-dbs.smk",
    "rules/mapping.smk",
    "rules/magrec.smk",
    "rules/genome-quality.smk",
    "rules/genome-tax.smk",
    "scripts/drep-clusters.py",
    "envs/bit.yaml",
    "envs/checkm2.yaml",
    "envs/classify-seqs.yaml",
    "envs/drep.yaml",
    "envs/gtdb-tk.yaml",
    "envs/mapping.yaml",
    "envs/metabat.yaml",
)


def binning_params(group, run_group, tax_group):
    """ the settings of binning and MAG recovery (rules/magrec.smk), also used by the metagenomics workflow """

    return [
        Param("min_completion", 90, type=float, group=group, metavar="<NUM>",
              help="Minimum estimated percent completion (checkm2) for a bin to be considered a MAG"),
        Param("max_redundancy", 10, type=float, group=group, metavar="<NUM>",
              help="Maximum estimated percent redundancy (checkm2) for a bin to be considered a MAG"),
        Param("derep_ani", 95, type=float, group=group, metavar="<NUM>",
              help="Percent average nucleotide identity at which the MAGs of all samples are dereplicated (dRep) before each "
                   "sample's reads are mapped to them for the MAG-by-sample tables (95 is about the species level; at 99, "
                   "strains of one species stay apart and the reads they share are split between them at random)"),
        Param("MAG_min_read_pid", 95, type=float, group=group, metavar="<NUM>",
              help="Minimum percent identity of a read mapped to the dereplicated MAGs for it to count in the MAG-by-sample "
                   "tables (reads of organisms without a MAG otherwise land on their closest relative)"),
        Param("gtdbtk_cpus", 8, type=int, group=run_group, metavar="<INT>",
              help="Number of cpus for GTDB-Tk"),
        Param("gtdbtk_scratch_dir", "", group=run_group, metavar="<DIR>",
              help="Directory GTDB-Tk's pplacer step can use for scratch space instead of memory (slower, but needs far less memory; not used by default; lower --gtdbtk-mem to match)"),

        Param("mapping_mem", 25000, type=int, group=run_group, metavar="<INT>", help="bowtie2 read-mapping memory in MB"),
        Param("checkm2_mem", 16000, type=int, group=run_group, metavar="<INT>", help="checkm2 memory in MB"),
        Param("gtdbtk_mem", 175000, type=int, group=run_group, metavar="<INT>", help="GTDB-Tk memory in MB"),
    ] + ([] if tax_group is None else labelling_params(tax_group, run_group, "only used when GTDB-Tk isn't available", genomes_only=True))


REQ = "Required Parameters"
GEN = "General Parameters"
MAG = "Binning and MAG Parameters"
TAX = "Taxonomic-assignment Parameters"

PARAMS = [
    Param("reads_dir", None, required=True, kind="dir", group=REQ, main=True, short="-r",
          help="Directory holding the gzipped fastq files the assemblies were made from, already quality-filtered "
               "(not normalized)"),
    Param("assemblies", None, required=True, kind="path", group=REQ, main=True, short="-a", metavar="<PATH>",
          help="Either a directory holding one assembly per sample, named <sample ID><extension> (see -x), each to be binned "
               "with its sample's reads, or a single co-assembly fasta to be binned with the reads of all the samples"),

    Param("assembly_extension", ".fasta", group=GEN, main=True, short="-x", metavar="<STR>",
          help="Extension of the per-sample assemblies in the -a directory, including the period (e.g., '.fasta', "
               "'-assembly.fasta', '.fa')"),
    Param("ref_dbs_dir", None, kind="path", group=GEN, main=True, short="-d", metavar="<DIR>", help=REF_DBS_DIR_HELP),
    Param("sample_ids_file", None, kind="file", group=GEN, main=True, short="-s",
          help="Single-column file with the sample IDs to process if not wanting all of those in the reads directory "
               "(a sample ID is what comes before the R1/R2 designation in the read-file names, e.g., before _R1.fastq.gz, "
               "or before the extension when the reads are single-end)"),
    Param("single_end", False, kind="flag", group=GEN, main=True,
          help="The reads are single-end rather than paired-end"),

    Param("threads", 8, type=int, group=GEN, main=True, short="-t", metavar="<INT>",
          help="Number of threads per job for the steps that use multiple (may be multiplied by --jobs)"),
    Param("output_prefix", "", group=GEN, metavar="<STR>",
          help="Prefix added to the names of the output tables (a '-' is automatically added unless it ends in '-' or '_')"),

    *binning_params(MAG, RUN_GROUP, TAX),
]
