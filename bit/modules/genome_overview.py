"""
Builds the overview table of a set of genomes (input genomes, bins, or MAGs) for the workflows: assembly summary
stats, estimated completeness and redundancy, and taxonomy when there is any, one row per genome.
"""

import csv

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# the columns kept from `bit summarize-assembly -t`
STATS_COLUMNS = ["Total contigs", "Total length", "GC content", "Maximum contig length", "N50", "L50",
                 "Num. contigs >= 10000", "Num. contigs >= 50000", "Num. contigs >= 100000"]

QUALITY_COLUMNS = ["est. completeness", "est. redundancy"]


def _rows(path):
    with open(path, newline="") as f:
        return [row for row in csv.reader(f, delimiter="\t") if row and any(cell.strip() for cell in row)]


def read_assembly_stats(path):
    """
    Reads the table written by `bit summarize-assembly -t` (one row per assembly).
    Returns (ordered assembly IDs, {ID: {column: value}}), or (None, the file's text) if it isn't such a table,
    as when a workflow wrote a note that nothing was recovered instead.
    """
    rows = _rows(path)
    if not rows or rows[0][0] != "Assembly":
        with open(path) as f:
            return None, f.read()
    header = rows[0]
    missing = [col for col in STATS_COLUMNS if col not in header]
    if missing:
        raise ValueError(f"'{path}' doesn't have the expected column(s): {', '.join(missing)}")
    stats = {row[0]: dict(zip(header, row)) for row in rows[1:]}
    return [row[0] for row in rows[1:]], stats


def read_checkm2(path):
    """checkm2's quality_report.tsv -> {genome: (completeness, contamination)}"""
    rows = _rows(path)
    if not rows or rows[0][:3] != ["Name", "Completeness", "Contamination"]:
        return {}
    return {row[0]: (row[1], row[2]) for row in rows[1:]}


def read_eukcc(path):
    """the workflow's combined EukCC table (Assembly, Est. Comp., Est. Redund.) -> {genome: (completeness, redundancy)}"""
    return {row[0]: (row[1], row[2]) for row in _rows(path)[1:] if len(row) >= 3}


def split_gtdb_lineage(lineage):
    """'d__Bacteria;p__Bacillota;...;s__' -> the seven rank names, 'NA' where there is none"""
    names = []
    for part in lineage.strip().split(";"):
        part = part.strip()
        name = part[3:] if len(part) >= 3 and part[1:3] == "__" else ""
        names.append(name if name else "NA")
    # anything that isn't a lineage (e.g. "Unclassified", "Unclassified Bacteria")
    if not any(len(part.strip()) >= 3 and part.strip()[1:3] == "__" for part in lineage.split(";")):
        words = lineage.split()
        names = [words[1]] if len(words) == 2 and words[0] == "Unclassified" else []
    return (names + ["NA"] * len(RANKS))[:len(RANKS)]


def read_gtdbtk(paths):
    """GTDB-Tk's summary tables (or their first two columns) -> {genome: the seven rank names}"""
    tax = {}
    for path in paths:
        for row in _rows(path):
            if row[0] == "user_genome" or len(row) < 2:
                continue
            tax[row[0]] = split_gtdb_lineage(row[1])
    return tax


def read_rank_table(path):
    """a table of genome and then the seven rank names (as characterize-genomes writes for eukaryotic genomes) -> {genome: rank names}"""
    tax = {}
    for row in _rows(path)[1:]:
        names = [name.strip() if name.strip() not in ("", "NA", "Not Assigned") else "NA" for name in row[1:8]]
        tax[row[0]] = (names + ["NA"] * len(RANKS))[:len(RANKS)]
    return tax


def build_overview(assembly_ids, stats, quality, taxonomy=None):
    """
    Returns the header and rows of the overview table. With taxonomy the rows are ordered by it (then by ID);
    without, they keep the order of the assembly stats.
    """
    header = ["Assembly"] + STATS_COLUMNS + QUALITY_COLUMNS + (RANKS if taxonomy is not None else [])
    rows = []
    for ID in assembly_ids:
        row = [ID] + [stats[ID][col] for col in STATS_COLUMNS] + list(quality.get(ID, ("NA", "NA")))
        if taxonomy is not None:
            row += taxonomy.get(ID, ["NA"] * len(RANKS))
        rows.append(row)
    if taxonomy is not None:
        rows.sort(key=lambda row: (row[-len(RANKS):], row[0]))
    return header, rows


def write_overview(stats_path, out_path, checkm2_path=None, eukcc_path=None, gtdbtk_paths=None, rank_table_path=None):
    """
    Writes the overview table. If the assembly-stats file is a note that there was nothing to summarize rather
    than a table, the note is written instead.
    """
    assembly_ids, stats = read_assembly_stats(stats_path)
    if assembly_ids is None:
        with open(out_path, "w") as out:
            out.write(stats)
        return

    quality = read_checkm2(checkm2_path) if checkm2_path else read_eukcc(eukcc_path) if eukcc_path else {}

    taxonomy = None
    if gtdbtk_paths is not None:
        taxonomy = read_gtdbtk(gtdbtk_paths)
    elif rank_table_path is not None:
        taxonomy = read_rank_table(rank_table_path)

    header, rows = build_overview(assembly_ids, stats, quality, taxonomy)
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(cell) for cell in row) + "\n")
