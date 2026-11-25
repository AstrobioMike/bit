import gzip
import json
import zipfile
import numpy as np # type: ignore
import pytest # type: ignore
from bit.modules.kmer_lookup import read_sketches, read_genome_index, write_table, find_shared_kmers


@pytest.fixture
def table(tmp_path):
    prefix = str(tmp_path / "ref")
    write_table(prefix, [("gA", [50, 10, 30, 20, 40]), ("gB", [30, 60, 70]), ("gC", []), ("gD", [2 ** 63 + 5, 7])])
    return prefix


def test_write_table_sorts_each_block_and_indexes_it(table):
    assert read_genome_index(table) == (["gA", "gB", "gC", "gD"], [0, 5, 8, 8], [5, 3, 0, 2])
    assert np.fromfile(table + ".hashes.u64", dtype="<u8").tolist() == [10, 20, 30, 40, 50, 30, 60, 70, 7, 2 ** 63 + 5]


@pytest.mark.parametrize("threads", [1, 4])
def test_find_shared_kmers(table, threads):
    queries = [np.array([10, 20, 99], dtype=np.uint64),          # two with gA
               np.array([30, 60, 2 ** 63 + 5], dtype=np.uint64), # 30 is in gA and gB; 60 in gB; the large one in gD
               np.array([1, 2, 3], dtype=np.uint64),             # nothing
               np.array([], dtype=np.uint64)]
    found = sorted(find_shared_kmers(table, queries, threads=threads))
    assert found == [(0, "gA", 2), (1, "gA", 1), (1, "gB", 2), (1, "gD", 1)]


def test_find_shared_kmers_counts_a_hash_for_every_query_holding_it(table):
    queries = [np.array([30, 40], dtype=np.uint64), np.array([30], dtype=np.uint64), np.array([30, 70], dtype=np.uint64)]
    found = sorted(find_shared_kmers(table, queries))
    assert found == [(0, "gA", 2), (0, "gB", 1), (1, "gA", 1), (1, "gB", 1), (2, "gA", 1), (2, "gB", 2)]


def test_find_shared_kmers_without_query_hashes(table):
    assert list(find_shared_kmers(table, [np.array([], dtype=np.uint64)])) == []
    assert list(find_shared_kmers(table, [])) == []


def test_read_sketches(tmp_path):
    path = tmp_path / "q.sig.zip"
    records = [{"name": "contig1", "signatures": [{"mins": [5, 9], "ksize": 31}]}]
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("SOURMASH-MANIFEST.csv", "x")
        archive.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(records).encode()))
        archive.writestr("signatures/b.sig", json.dumps([{"name": "contig2", "signatures": [{"mins": [], "ksize": 31}]}]))
    names, hashes = read_sketches(path)
    assert names == ["contig1", "contig2"]
    assert [h.tolist() for h in hashes] == [[5, 9], []]
    assert hashes[0].dtype == np.uint64
