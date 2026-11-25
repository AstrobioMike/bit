# The reference databases: where they are. They are set up by bit (`bit data get <name>`, or by `bit run-wf` and the
# other workflow commands before snakemake starts), not by these workflows; see bit.modules.ref_dbs.
# Pulled in with `include:` after rules/common.smk, by the workflows that use reference databases. Defines ref_root
# (the reference-database directory, ending with "/"), the paths of the databases' contents, and require_ref_dbs().

import sys
from bit.modules.ref_dbs import REF_DBS, GTDBTK_RUNS_HERE, GTDBTK_NOT_HERE_NOTE, RefDbsDirNotSet, missing, ref_dbs_root

def _ref_root(fallback = None):
    """ as bit resolved it (in the config), or else the BIT_REF_DBS_DIR variable (a saved copy run with snakemake) """
    try:
        return str(ref_dbs_root(config.get("ref_dbs_dir"))) + "/"
    except RefDbsDirNotSet as problem:
        if fallback:
            return fallback
        exit_with_message(str(problem))

def require_ref_dbs(names):
    """ stops with what to run when a database the run uses isn't set up (bit sets them up itself before starting) """
    if config.get("ref_dbs_set_up_by_bit"):
        return
    gone = missing(ref_root, names)
    if gone:
        exit_with_message(f"These reference databases aren't set up in {ref_root}: {', '.join(gone)}\n\n    Set them up with:\n" +
                          "\n".join(f"      bit data get {name} -d {ref_root}" for name in gone))

def ref_db_dir(name):
    return str(REF_DBS[name].path(ref_root)) + "/"

def note_gtdbtk_skipped():
    """ said by bit when it starts the run; here only for a saved copy run with snakemake directly """
    if not GTDBTK_RUNS_HERE and not config.get("ref_dbs_set_up_by_bit"):
        print(color_text(f"\n  Note: {GTDBTK_NOT_HERE_NOTE}\n", "yellow"), file = sys.stderr)

ref_root = _ref_root(globals().get("ref_root_fallback"))

checkm2_db_file = ref_db_dir("checkm2-db") + "uniref100.KO.1.dmnd"
gtdbtk_db_dir = ref_db_dir("gtdbtk-db")
