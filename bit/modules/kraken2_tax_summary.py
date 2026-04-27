#!/usr/bin/env python
import pandas as pd
import numpy as np
import re
from bit.modules.general import check_files_are_found, report_failure


RANK_MAP = {
    'D': 'domain',
    'P': 'phylum',
    'C': 'class',
    'O': 'order',
    'F': 'family',
    'G': 'genus',
    'S': 'species',
}

STD_RANKS = ['domain','phylum','class','order','family','genus','species']


def kraken2_to_taxon_summaries(input_report, output_tsv):

    preflight_checks(input_report)
    df = parse_report(input_report)
    df = refine_df(df)
    df.to_csv(output_tsv, sep='\t', index=False, float_format="%.2f")


def preflight_checks(input_report):

    check_files_are_found([input_report])


def parse_report(input_report):

    rows = []
    current_rank = {r: 'NA' for r in STD_RANKS}

    with open(input_report) as in_report:
        for line in in_report:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            rec = parse_report_line(line)
            if rec is None:
                continue

            taxid = rec['taxid']
            taxon_reads = rec['taxon_reads']
            name = rec['name']
            rank_code = normalize_rank_code(rec['rank'], name)

            # handling unclassified
            if taxid == 0 or rank_code == 'U':
                rows.append({
                    'taxid': 0,
                    **{r: 'Unclassified' for r in STD_RANKS},
                    'read_counts': taxon_reads,
                })
                continue

            # if standard rank, update that rank and clear lower ones
            if rank_code in RANK_MAP:
                rank_name = RANK_MAP[rank_code]
                current_rank[rank_name] = name
                lower = False
                for r in STD_RANKS:
                    if r == rank_name:
                        lower = True
                        continue
                    if lower:
                        current_rank[r] = 'NA'
            else:
                # if non-standard rank (e.g., strains, subspecies, or '-')
                pass

            # create a row for every taxon line
            rows.append({
                'taxid': taxid,
                **{r: current_rank[r] for r in STD_RANKS},
                'read_counts': taxon_reads,
            })

    return pd.DataFrame(rows)


def parse_report_line(line):
    # kraken2 report columns (space-separated with name possibly containing spaces)
    # percent clade_reads taxon_reads rank taxid name (with leading spaces for depth in front of name)
    parts = re.split(r'\s+', line.strip(), maxsplit=5)
    if len(parts) < 6:
        # trying tab split as a fallback
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 6:
            report_failure(f"Could not parse kraken2 report line: {line}")

    clade_reads, taxon_reads, rank, taxid, name = parts[1:6]
    m = re.match(r'^(\s*)(.*)$', name)
    clean_name = m.group(2)

    return {
        'clade_reads': int(clade_reads),
        'taxon_reads': int(taxon_reads),
        'rank': rank,
        'taxid': int(taxid),
        'name': clean_name,
    }


def normalize_rank_code(rank_code, name):
    # this is here because gtdb kraken2 outputs seem to have 'R1' for domain instead of 'D'
    DOMAIN_LIKE_NAMES = {'Bacteria', 'Archaea', 'Eukaryota', 'Viruses'}
    if rank_code == 'R1' and name in DOMAIN_LIKE_NAMES:
        return 'D'
    return rank_code

def refine_df(df):

    # combining any duplicate taxids
    agg_cols = ['read_counts']
    keep_cols = ['taxid'] + STD_RANKS
    df = (df
          .groupby(keep_cols, as_index=False)[agg_cols].sum())

    # calculating overall percent_of_reads (kraken ones are clade percents)
    total_reads = df['read_counts'].sum()
    if total_reads > 0:
        df['percent_of_reads'] = df['read_counts'] / total_reads * 100
    else:
        df['percent_of_reads'] = 0.0

    df['percent_of_reads'] = df['percent_of_reads'].round(2)

    # dropping rows if they hold no reads
    df = df[df['read_counts'] > 0]

    # sorting by taxid
    df.sort_values(by=['taxid'], inplace=True)

    # re-ordering columns
    out_cols = ['taxid'] + STD_RANKS + ['read_counts', 'percent_of_reads']

    df = sort_df(df)

    return df[out_cols]


def sort_df(df):

    # sorting so unclassified is first, then any all NAs, then the rest based on lineage

    # creating groups
    is_unclassified = df['taxid'].eq(0)
    all_na_lineage = df[STD_RANKS].eq('NA').all(axis=1)

    # making column holding group ID: 0 = Unclassified, 1 = all-NAs, 2 = the rest
    df['sort_group'] = np.select(
        [is_unclassified, all_na_lineage],
        [0, 1],
        default=2
    )

    df.sort_values(
        by=['sort_group'] + STD_RANKS + ['taxid'],
        inplace=True,
        kind='mergesort'
    )

    df.drop(columns=['sort_group'], inplace=True)

    return df
