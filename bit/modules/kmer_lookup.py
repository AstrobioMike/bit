"""
Finding the reference genomes that share sampled k-mers with query sequences, for `bit classify-seqs`.

The reference is a flat table made from sourmash sketches (FracMinHash: the k-mer hashes under a
fraction of the hash range), two files sharing a prefix:
  <prefix>.hashes.u64    every genome's hashes, each genome's block sorted (raw little-endian uint64)
  <prefix>.genomes.tsv   genome name, offset of its block (in hashes), number of hashes; in block order

All the queries' hashes are pooled into one sorted array and each genome's block is intersected with
it, so the cost is one read of the table plus the matches, whatever the number of queries.
"""

import gzip
import json
import zipfile
from concurrent.futures import ThreadPoolExecutor
import numpy as np # type: ignore

GENOMES_PER_TASK = 500


def read_sketches(sig_zip):
    """ the names and hash arrays of the sketches in a sourmash signature zip """

    names, hashes = [], []
    with zipfile.ZipFile(sig_zip) as archive:
        for member in archive.namelist():
            if not member.startswith("signatures/"):
                continue
            raw = archive.read(member)
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            for record in json.loads(raw):
                for signature in record["signatures"]:
                    names.append(record["name"])
                    hashes.append(np.array(signature["mins"], dtype=np.uint64))

    return names, hashes


def read_genome_index(table_prefix):
    """ genome names, block offsets and block sizes of a k-mer table """

    names, offsets, counts = [], [], []
    with open(f"{table_prefix}.genomes.tsv") as index:
        for line in index:
            name, offset, count = line.rstrip("\n").split("\t")
            names.append(name)
            offsets.append(int(offset))
            counts.append(int(count))

    return names, offsets, counts


def write_table(table_prefix, genomes):
    """ genomes: (name, hashes) pairs; writes the two files of a k-mer table """

    offset = 0
    with open(f"{table_prefix}.hashes.u64", "wb") as table, open(f"{table_prefix}.genomes.tsv", "w") as index:
        for name, hashes in genomes:
            block = np.sort(np.asarray(hashes, dtype=np.uint64))
            block.astype("<u8").tofile(table)
            index.write(f"{name}\t{offset}\t{len(block)}\n")
            offset += len(block)


def find_shared_kmers(table_prefix, query_hashes, threads=1):
    """
    query_hashes: one array of hashes per query
    Yields (query index, genome name, number of shared hashes) for every query and genome sharing any.
    """

    sizes = [len(h) for h in query_hashes]
    if not sum(sizes):
        return

    owner = np.repeat(np.arange(len(query_hashes), dtype=np.uint32), sizes)
    pooled = np.concatenate([np.asarray(h, dtype=np.uint64) for h in query_hashes])
    order = np.argsort(pooled, kind="stable")
    pooled, owner = pooled[order], owner[order]
    # a hash may belong to several queries: unique[i] is held by owner[first[i]:last[i]]
    unique, first = np.unique(pooled, return_index=True)
    last = np.append(first[1:], len(pooled))
    lowest, highest = unique[0], unique[-1]

    names, offsets, counts = read_genome_index(table_prefix)
    table = np.memmap(f"{table_prefix}.hashes.u64", dtype="<u8", mode="r") if sum(counts) else np.empty(0, dtype="<u8")

    def scan(genome_indices):
        found = []
        for g in genome_indices:
            block = table[offsets[g]:offsets[g] + counts[g]]
            block = block[np.searchsorted(block, lowest):np.searchsorted(block, highest, side="right")]
            position = np.searchsorted(unique, block)
            position[position == len(unique)] = 0
            hit = position[unique[position] == block]
            if len(hit):
                if (last[hit] - first[hit]).max() > 1:
                    owners = np.concatenate([owner[first[p]:last[p]] for p in hit])
                else:
                    owners = owner[first[hit]]
                queries, shared = np.unique(owners, return_counts=True)
                found.append((g, queries, shared))
        return found

    tasks = [range(i, min(i + GENOMES_PER_TASK, len(names))) for i in range(0, len(names), GENOMES_PER_TASK)]
    with ThreadPoolExecutor(max(1, threads)) as pool:
        for found in pool.map(scan, tasks):
            for g, queries, shared in found:
                for query, count in zip(queries.tolist(), shared.tolist()):
                    yield query, names[g], count
