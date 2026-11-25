from bit.modules.workflows import Param, RUN_GROUP
from bit.smk.workflows.classify_seqs.params import labelling_params
from bit.modules.ref_dbs import REF_DBS_DIR_HELP
from bit.smk.workflows.magrec.params import binning_params, binning_ref_dbs

DESCRIPTION = """
    This is an assembly-based metagenomics workflow for short-read Illumina data. It does quality
    trimming/filtering (bbduk), assembly of each sample (megahit), gene-calling (prodigal),
    functional annotation (KOFamScan), taxonomic classification of genes and contigs, read-mapping
    and coverage (bowtie2, bit cov-stats), and by default binning and MAG recovery (metabat2,
    checkm2, GTDB-Tk, KEGG-Decoder; the same steps as `bit run-wf magrec`). With --binning-only
    it goes from QC and assembly straight to binning and MAG recovery. The reference databases are
    set up the first time they are needed.
    """

EXAMPLE = "bit run-wf metagenomics -r reads-dir/"

OUTPUTS_NOTE = "The combined tables are in the 'combined-outputs/' subdirectory."
DB_SETUP_SETTINGS = ("skip_binning", "binning_only")


def ref_dbs_needed(settings):
    """ the reference databases (bit.modules.ref_dbs) a run with these settings uses """
    annotation = [] if settings["binning_only"] else ["tax-db", "kofamscan-db"]
    binning = [] if settings["skip_binning"] else binning_ref_dbs()
    return list(dict.fromkeys(annotation + binning))


def runs_gtdbtk(settings):
    """ whether a run with these settings would use GTDB-Tk (where it runs) """
    return not settings["skip_binning"]

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "rules/read-qc.smk",
    "rules/reads.smk",
    "rules/genome-quality.smk",
    "rules/genome-tax.smk",
    "rules/ref-dbs.smk",
    "rules/mapping.smk",
    "rules/magrec.smk",
    "rules/classify-seqs.smk",
    "scripts/multiqc.config",
    "scripts/drep-clusters.py",
    "envs/bit.yaml",
    "envs/checkm2.yaml",
    "envs/classify-seqs.yaml",
    "envs/drep.yaml",
    "envs/gtdb-tk.yaml",
    "envs/keggdecoder.yaml",
    "envs/kofamscan.yaml",
    "envs/mapping.yaml",
    "envs/megahit.yaml",
    "envs/metabat.yaml",
    "envs/prodigal.yaml",
    "envs/qc.yaml",
)

REQ = "Required Parameters"
GEN = "General Parameters"
TAX = "Taxonomic-assignment Parameters"
MAG = "Binning and MAG Parameters"

PARAMS = [
    Param("reads_dir", None, required=True, kind="dir", group=REQ, main=True, short="-r",
          help="Directory holding the input gzipped fastq files"),
    Param("ref_dbs_dir", None, kind="path", group=GEN, main=True, short="-d", metavar="<DIR>", help=REF_DBS_DIR_HELP),

    Param("sample_ids_file", None, kind="file", group=GEN, main=True, short="-s",
          help="Single-column file with the sample IDs to process if not wanting all of those in the reads directory "
               "(a sample ID is what comes before the R1/R2 designation in the read-file names, e.g., before _R1.fastq.gz, "
               "or before the extension when the reads are single-end)"),
    Param("single_end", False, kind="flag", group=GEN, main=True,
          help="The reads are single-end rather than paired-end"),

    Param("threads", 8, type=int, group=GEN, main=True, short="-t", metavar="<INT>",
          help="Number of threads per job for the steps that use multiple (may be multiplied by --jobs)"),
    Param("skip_binning", False, kind="flag", group=GEN, main=True,
          help="Don't do binning and MAG recovery"),
    Param("binning_only", False, kind="flag", group=GEN, main=True,
          help="Only QC, assembly, and binning and MAG recovery (for starting from assemblies already made, see `bit run-wf magrec`)"),
    Param("output_prefix", "", group=GEN, metavar="<STR>",
          help="Prefix added to the names of the combined output files (a '-' is automatically added unless it ends in '-' or '_')"),

    Param("swift_1s", False, kind="flag", group=GEN,
          help="The libraries were prepared with the Swift 1S kit (its adaptase tails are trimmed off)"),


    *labelling_params(TAX, RUN_GROUP, "only used when GTDB-Tk isn't available"),

    *binning_params(MAG, RUN_GROUP, None),

    Param("megahit_mem", 32000, type=int, group=RUN_GROUP, metavar="<INT>", help="megahit memory in MB (the most megahit is allowed to use per assembly; a failed assembly is retried with double this, then quadruple, up to the memory the run has)"),
    Param("pileup_mem", 5000, type=int, group=RUN_GROUP, metavar="<INT>", help="Coverage calculation (bbtools' pileup.sh) memory in MB"),
    Param("kofamscan_mem", 5000, type=int, group=RUN_GROUP, metavar="<INT>", help="KOFamScan memory in MB"),
]

# the run settings come from several tables above; list them in the order of the steps
_RUN_ORDER = ["gtdbtk_cpus", "gtdbtk_scratch_dir", "megahit_mem", "mapping_mem", "pileup_mem", "kofamscan_mem",
              "search_mem", "classify_seqs_mem", "checkm2_mem", "gtdbtk_mem"]
assert sorted(_RUN_ORDER) == sorted(p.name for p in PARAMS if p.group == RUN_GROUP)
PARAMS.sort(key=lambda p: _RUN_ORDER.index(p.name) if p.group == RUN_GROUP else -1)
