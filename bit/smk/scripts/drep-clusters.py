#!/usr/bin/env python
"""
Writes each MAG and the representative standing for it ("MAG_ID<tab>representative", representatives first) from a
dRep dereplicate output directory: its clusters (data_tables/Cdb.csv) and their representatives (data_tables/Wdb.csv).
"""

import csv
import os
import sys


def ID(file_name):
    return file_name[:-len(".fasta")] if file_name.endswith(".fasta") else file_name


def clusters(drep_dir):
    tables = os.path.join(drep_dir, "data_tables")
    with open(os.path.join(tables, "Wdb.csv")) as f:
        representative_of = {row["cluster"]: row["genome"] for row in csv.DictReader(f)}
    with open(os.path.join(tables, "Cdb.csv")) as f:
        return {ID(row["genome"]): ID(representative_of[row["secondary_cluster"]]) for row in csv.DictReader(f)}


if __name__ == "__main__":
    found = clusters(sys.argv[1])
    for MAG in sorted(found, key=lambda m: (found[m], m != found[m], m)):
        print(f"{MAG}\t{found[MAG]}")
