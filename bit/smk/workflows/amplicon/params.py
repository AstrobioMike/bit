from bit.modules.workflows import Param, RUN_GROUP

DESCRIPTION = """
    This is an amplicon workflow for Illumina data (16S, 18S or ITS). It removes primers
    (cutadapt), then quality-filters, infers amplicon sequence variants (ASVs), removes
    chimeras and assigns taxonomy with DADA2 and DECIPHER, giving ASV sequences, count and
    taxonomy tables, and fastqc/multiqc reports of the raw and filtered reads.
    """

EXAMPLE = "bit run-wf amplicon -r reads-dir/ --target-region 16S --f-primer ^GTGCCAGCMGCCGCGGTAA --r-primer ^GGACTACHVGGGTWTCTAA"

OUTPUTS_NOTE = "The ASV sequences and the count and taxonomy tables are in the 'final-outputs/' subdirectory."
DB_SETUP_SETTINGS = ()
# without a reference-database directory (-d/--ref-dbs-dir or BIT_REF_DBS_DIR), the references go in the output directory
REF_DBS_IN_OUTPUT_DIR = True


def ref_dbs_needed(settings):
    """ the reference databases (bit.modules.ref_dbs) a run uses """
    return ["amplicon-refs"]

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "rules/read-qc.smk",
    "rules/reads.smk",
    "rules/ref-dbs.smk",
    "scripts/multiqc.config",
    "envs/R.yaml",
    "envs/cutadapt.yaml",
    "envs/qc.yaml",
)

REQ = "Required Parameters"
PRI = "Primer-trimming Parameters"
FIL = "Quality-filtering and Merging Parameters (DADA2)"
GEN = "General Parameters"

PARAMS = [
    Param("reads_dir", None, required=True, kind="dir", group=REQ, main=True, short="-r",
          help="Directory holding the input gzipped fastq files"),

    Param("sample_ids_file", None, kind="file", group=GEN, main=True, short="-s",
          help="Single-column file with the sample IDs to process if not wanting all of those in the reads directory "
               "(a sample ID is what comes before the R1/R2 designation in the read-file names, e.g., before _R1.fastq.gz, "
               "or before the extension when the reads are single-end)"),
    Param("single_end", False, kind="flag", group=GEN, main=True,
          help="The reads are single-end rather than paired-end"),
    Param("ref_dbs_dir", None, kind="path", group=GEN, main=True, short="-d", metavar="<DIR>",
          help="Directory to keep the taxonomy references in (~575 MB), so they are downloaded once and re-used (default: "
               "the one in the BIT_REF_DBS_DIR variable, set with `bit data locations set`; without either, they are "
               "downloaded into the output directory of each run)"),

    Param("target_region", "16S", choices=["16S", "18S", "ITS"], group=GEN, main=True,
          help="The target region, which sets the reference used for taxonomy (16S: SILVA; 18S: PR2; ITS: UNITE)"),

    Param("f_primer", "^GTGCCAGCMGCCGCGGTAA", group=PRI, main=True, metavar="<SEQ>",
          help="Forward primer as cutadapt takes it (a leading ^ anchors it to the start of the read; degenerate bases are fine)"),
    Param("r_primer", "^GGACTACHVGGGTWTCTAA", group=PRI, main=True, metavar="<SEQ>",
          help="Reverse primer as cutadapt takes it (with single-end reads, it is looked for towards the ends of the reads)"),
    Param("skip_primer_trimming", False, kind="flag", group=PRI, main=True,
          help="Don't trim primers (e.g., they were already removed)"),
    Param("primers_not_linked", False, kind="flag", group=PRI,
          help="Only look for each read's own primer at its start. By default each read is also checked for the reverse complement "
               "of the other primer towards its end (cutadapt's linked adapters), as happens when reads are longer than the amplicon"),
    Param("keep_untrimmed", False, kind="flag", group=PRI,
          help="Keep reads in which the primers weren't found (they are discarded by default)"),
    Param("min_cutadapt_len", 150, type=int, group=PRI, metavar="<INT>",
          help="Reads shorter than this after primer trimming are discarded"),

    Param("left_trunc", 0, type=int, group=FIL, main=True, metavar="<INT>",
          help="Length the forward reads are truncated to (dada2's truncLen; 0 for no truncation)"),
    Param("right_trunc", 0, type=int, group=FIL, main=True, metavar="<INT>",
          help="Length the reverse reads are truncated to (0 for no truncation; not used with single-end reads)"),
    Param("left_maxEE", 1, type=int, group=FIL, metavar="<INT>",
          help="Maximum expected errors allowed in a forward read (dada2's maxEE)"),
    Param("right_maxEE", 1, type=int, group=FIL, metavar="<INT>",
          help="Maximum expected errors allowed in a reverse read (not used with single-end reads)"),
    Param("concatenate_reads_only", False, kind="flag", group=FIL,
          help="Concatenate the forward and reverse reads rather than merging them (for when they don't overlap; not used with single-end reads)"),

    Param("output_prefix", "", group=GEN, metavar="<STR>",
          help="Prefix added to the names of the output files (a '-' is automatically added unless it ends in '-' or '_')"),
    Param("R_mem", 16000, type=int, group=RUN_GROUP, metavar="<INT>",
          help="R step (DADA2 and DECIPHER's IdTaxa) memory in MB (runs with many samples or deep sequencing may need more)"),
]
