from bit.modules.workflows import Param, RUN_GROUP

DESCRIPTION = """
    This workflow downloads reads from NCBI's Sequence Read Archive (SRA) given run accessions
    (those starting with SRR, ERR, or DRR), using prefetch and fasterq-dump, and gzips the
    resulting fastq files.
    """

EXAMPLE = "bit dl-sra-reads -a target-sra-accs.txt"

OUTPUTS_NOTE = "The gzipped fastq files are in the 'fastq-files/' subdirectory."

# what the steps of a run download from, by rule; checked before a run with those steps to do starts
# (bit.modules.connectivity). prefetch asks NCBI's locator where a run's data is (it is then served from NCBI or one
# of its cloud mirrors)
SRA_LOCATOR = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
DOWNLOADS_FROM = {"prefetch": (SRA_LOCATOR,)}

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "envs/sra-dl.yaml",
)

REQ = "Required Parameters"
GEN = "General Parameters"

PARAMS = [
    Param("accessions_file", None, required=True, kind="file", group=REQ, main=True, short="-a", metavar="<FILE>",
          help="File with the target SRA run accessions, one per line (these should start with SRR, ERR, or DRR). To have runs that "
               "belong to the same sample combined, make it two columns: the sample name, then the run accession "
               "(sample names can't start with SRR, ERR, or DRR)"),

    Param("threads", 8, type=int, group=GEN, main=True, short="-t", metavar="<INT>",
          help="Number of threads per job for the steps that use multiple (may be multiplied by --jobs)"),

    Param("prefetch_max_size", "500G", group=RUN_GROUP, metavar="<SIZE>",
          help="Largest SRA object prefetch will download (its --max-size setting)"),
    Param("keep_sra_files", False, kind="flag", group=RUN_GROUP,
          help="Keep the downloaded SRA objects (in an 'sra-files/' subdirectory) rather than removing them once converted to fastq"),
]
