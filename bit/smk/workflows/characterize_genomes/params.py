from bit.modules.workflows import Param, RUN_GROUP
from bit.modules.ref_dbs import REF_DBS_DIR_HELP
from bit.smk.workflows.classify_seqs.params import labelling_params

DESCRIPTION = """
    This workflow summarizes genome assemblies giving you assembly summary stats, estimated completeness and redundancy,
    and taxonomy. For bacteria and archaea we use checkm2 and GTDB-Tk, for eukarya (specified with --euk) we use EukCC and taxonomy is
    done by `bit classify-seqs`. The reference databases are set up the first time they are needed.
    """

EXAMPLE = "bit characterize-genomes -g genomes-dir/"

OUTPUTS_NOTE = "The summary table is 'genome-summaries.tsv' (following the output prefix if one was given)."
DB_SETUP_SETTINGS = ("euk",)


def ref_dbs_needed(settings):
    """ the reference databases (bit.modules.ref_dbs) a run with these settings uses """
    from bit.modules.ref_dbs import GTDBTK_RUNS_HERE
    if settings["euk"]:
        return ["eukcc-db", "tax-db"]
    return ["checkm2-db", "gtdbtk-db" if GTDBTK_RUNS_HERE else "tax-db"]


def runs_gtdbtk(settings):
    """ whether a run with these settings would use GTDB-Tk (where it runs) """
    return not settings["euk"]

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "rules/genome-quality.smk",
    "rules/ref-dbs.smk",
    "rules/classify-seqs.smk",
    "rules/genome-tax.smk",
    "envs/bit.yaml",
    "envs/checkm2.yaml",
    "envs/classify-seqs.yaml",
    "envs/eukcc.yaml",
    "envs/gtdb-tk.yaml",
)

REQ = "Required Parameters"
GEN = "General Parameters"
TAX = "Taxonomic-assignment Parameters"

PARAMS = [
    Param("genomes_dir", None, required=True, kind="dir", group=REQ, main=True, short="-g",
          help="Directory holding the genome fasta files (not gzipped)"),
    Param("ref_dbs_dir", None, kind="path", group=GEN, main=True, short="-d", metavar="<DIR>", help=REF_DBS_DIR_HELP),

    Param("assembly_extension", ".fasta", group=GEN, main=True, short="-x", metavar="<STR>",
          help="Extension of the genome fasta files"),
    Param("euk", False, kind="flag", group=GEN, main=True,
          help="Specify this flag if the input genomes are eukaryotic"),

    Param("threads", 8, type=int, group=GEN, main=True, short="-t", metavar="<INT>",
          help="Number of threads per job for the steps that use multiple (may be multiplied by --jobs)"),
    Param("output_prefix", "", group=GEN, metavar="<STR>",
          help="Prefix added to the name of the output table (a '-' is automatically added unless it ends in '-' or '_')"),
    Param("keep_all_files", False, kind="flag", group=GEN,
          help="Keep the full outputs of the programs run rather than only the summary table"),

    Param("gtdbtk_scratch_dir", "", group=RUN_GROUP, metavar="<DIR>",
          help="Directory GTDB-Tk's pplacer step can use for scratch space instead of memory (slower, but needs far less memory; not used by default; lower --gtdbtk-mem to match)"),

    Param("checkm2_mem", 16000, type=int, group=RUN_GROUP, metavar="<INT>", help="checkm2 memory in MB"),
    Param("gtdbtk_mem", 175000, type=int, group=RUN_GROUP, metavar="<INT>", help="GTDB-Tk memory in MB"),
    Param("eukcc_mem", 50000, type=int, group=RUN_GROUP, metavar="<INT>", help="EukCC memory in MB"),
    *labelling_params(TAX, RUN_GROUP, "used with --euk, and when GTDB-Tk isn't available", genomes_only=True),
]
