from bit.modules.workflows import Param, RUN_GROUP
from bit.modules.classify_seqs import (DEFAULT_MIN_CONTAINMENT, DEFAULT_GENOME_MIN_CONTAINMENT,
                                      DEFAULT_IDENTITY_CUTOFFS)
from bit.modules.ref_dbs import REF_DBS_DIR_HELP

DESCRIPTION = """
    This program assigns taxonomy to nucleotide sequences, e.g., the contigs of metagenome assemblies
    or high-quality long reads such as HiFi, with a reference built from the GTDB (bacteria and
    archaea) and NCBI (eukarya and viruses). A species is assigned when a sequence shares enough of
    its sampled k-mers with a reference genome; otherwise genes are called and searched against the
    reference proteins, and the sequence gets the lowest taxon most of that evidence agrees on (genes
    are only searched for sequences their k-mers didn't assign). When several fasta files are given
    they are searched together, and each gets its own output table. With --genomes, each fasta file
    is instead classified as a whole (e.g. a genome, MAG or bin), all its contigs' k-mers and genes
    together, giving one table with a row per file. The fasta files can't be gzipped.
    """

EXAMPLE = "bit classify-seqs -f assembly-1.fasta assembly-2.fasta"

OUTPUTS_NOTE = "There is one '<input name>-tax.tsv' table per input fasta file (with --genomes, one 'genome-taxonomies.tsv')."

# the settings that matter to --setup-dbs-only
DB_SETUP_SETTINGS = ()


def ref_dbs_needed(settings):
    """ the reference databases (bit.modules.ref_dbs) a run uses """
    return ["tax-db"]

# files under bit/smk/ that this workflow uses
SHARED_FILES = (
    "rules/common.smk",
    "rules/ref-dbs.smk",
    "rules/classify-seqs.smk",
    "envs/classify-seqs.yaml",
)


def labelling_params(group, run_group, genome_note, genomes_only=False):
    """ the settings of the classification (rules/classify-seqs.smk, rules/genome-tax.smk), also used by the other workflows
        genome_note:  when --genome-min-containment applies in the workflow (in its help)
        genomes_only: the workflow only classifies whole genomes, so --min-containment (for single sequences) isn't offered """

    if genomes_only:
        containment = [
            Param("genome_min_containment", DEFAULT_GENOME_MIN_CONTAINMENT, type=float, group=group, metavar="<FLOAT>",
                  help="Minimum fraction of a genome's sampled k-mers a reference genome must hold for the genome to be "
                       f"assigned to that reference's species from k-mers ({genome_note})"),
        ]
    else:
        containment = [
            Param("min_containment", DEFAULT_MIN_CONTAINMENT, type=float, group=group, metavar="<FLOAT>",
                  help="Minimum fraction of a sequence's sampled k-mers a reference genome must hold for the sequence to be "
                       "assigned to that reference's species from k-mers"),
            Param("genome_min_containment", DEFAULT_GENOME_MIN_CONTAINMENT, type=float, group=group, metavar="<FLOAT>",
                  help="The same as --min-containment, except here for a genome to be taxonomically classified at the species "
                       f"level by k-mers ({genome_note})"),
        ]

    return containment + [
        Param("identity_cutoffs", ",".join(str(c) for c in DEFAULT_IDENTITY_CUTOFFS), group=group, metavar="<STR>",
              help="Minimum percent identity of a gene's best protein hit for it to support a species, genus, or family "
                   "taxonomic assignment (comma-separated)"),

        Param("search_mem", 32000, type=int, group=run_group, metavar="<INT>", help="DIAMOND protein search memory in MB"),
        Param("classify_seqs_mem", 40000, type=int, group=run_group, metavar="<INT>", help="k-mer lookup and taxonomic-assignment memory in MB"),
    ]


REQ = "Required Parameters"
GEN = "General Parameters"
LAB = "Taxonomic-assignment Parameters"

PARAMS = [
    Param("fasta_files", None, required=True, kind="paths", group=REQ, main=True, short="-f",
          help="The fasta file(s) of the sequences to classify: one, several (space-delimited), or a directory "
               "holding them (files ending with .fasta, .fa, .fna, .fas, or .fsa are taken from a directory)"),
    Param("genomes", False, kind="flag", group=GEN, main=True,
          help="Classify each fasta file as a whole (a genome, MAG or bin): the k-mers of all its contigs pooled, and all its "
               "genes voting together (its k-mer threshold is in the -H menu). Writes one table, 'genome-taxonomies.tsv', with a row "
               "per file"),
    Param("ref_dbs_dir", None, kind="path", group=GEN, main=True, short="-d", metavar="<DIR>", help=REF_DBS_DIR_HELP),

    Param("threads", 8, type=int, group=GEN, main=True, short="-t", metavar="<INT>",
          help="Number of threads per job for the steps that use multiple (may be multiplied by --jobs)"),
    Param("keep_all_files", False, kind="flag", group=GEN,
          help="Keep the intermediate files (called genes, search results, sketches) rather than only the output tables"),

    *labelling_params(LAB, RUN_GROUP, "used with --genomes"),
]
