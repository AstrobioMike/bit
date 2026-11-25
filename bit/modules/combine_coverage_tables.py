"""
The metagenomics workflow's combined tables: coverages from each sample's gene or contig table summed by KO
or by lineage, one column per sample, plus a version scaled to coverage per million (CPM: like a percent,
out of 1,000,000 instead of 100, so each sample's column sums to 1,000,000).
"""

from bit.modules.classify_seqs import RANKS, NO_LABEL

NOT_ANNOTATED = "Not annotated"


def read_table(path):
    """ the rows of a tab-delimited table with a header, as dicts """

    with open(path) as table:
        header = table.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in table if line.strip()]


def sum_coverages(tables, key_of):
    """
    tables: {sample name: rows}; key_of: row -> the key its coverage is summed under
    returns {key: [summed coverage per sample, in the order of tables]}
    """

    sums = {}
    for i, rows in enumerate(tables.values()):
        for row in rows:
            sums.setdefault(key_of(row), [0.0] * len(tables))[i] += float(row["coverage"])
    return sums


def per_million(sums, n_samples):
    """ each sample's coverages scaled to sum to 1,000,000 (left at 0 for a sample with no coverage) """

    totals = [sum(values[i] for values in sums.values()) for i in range(n_samples)]
    return {key: [v / t * 1_000_000 if t else 0.0 for v, t in zip(values, totals)] for key, values in sums.items()}


def write_combined(path, key_columns, samples, sums, order):
    with open(path, "w") as out:
        out.write("\t".join(list(key_columns) + list(samples)) + "\n")
        for key in sorted(sums, key=order):
            out.write("\t".join(list(key) + [f"{v:g}" for v in sums[key]]) + "\n")


def lineage_of(row):
    return tuple(row[rank] for rank in RANKS)


def lineage_order(lineage):
    """ alphabetical by lineage, with a rank left unassigned after the assigned ones (so unclassified is last) """

    return [(name == NO_LABEL, name) for name in lineage]


def KO_of(row):
    KO = row["KO_ID"]
    return (KO, row["KO_function"]) if KO.startswith("K") else (NOT_ANNOTATED, NOT_ANNOTATED)


def KO_order(key):
    return (key[0] == NOT_ANNOTATED, key[0])


def combine(tables, key_columns, key_of, order, out_path, cpm_out_path):
    """ tables: {sample name: table path}; samples whose table has no rows (no assembly or no genes) are left out """

    read = {sample: read_table(path) for sample, path in tables.items()}
    read = {sample: rows for sample, rows in read.items() if rows}
    sums = sum_coverages(read, key_of)
    write_combined(out_path, key_columns, read, sums, order)
    write_combined(cpm_out_path, key_columns, read, per_million(sums, len(read)), order)


def combine_gene_tables(tables, KO_out, KO_cpm_out, tax_out, tax_cpm_out):
    """ each sample's gene table (gene_ID, coverage, KO_ID, KO_function, then the ranks) summed by KO and by lineage """

    combine(tables, ("KO_ID", "KO_function"), KO_of, KO_order, KO_out, KO_cpm_out)
    combine(tables, RANKS, lineage_of, lineage_order, tax_out, tax_cpm_out)


def combine_contig_tables(tables, tax_out, tax_cpm_out):
    """ each sample's contig table (contig_ID, coverage, then the ranks) summed by lineage """

    combine(tables, RANKS, lineage_of, lineage_order, tax_out, tax_cpm_out)
