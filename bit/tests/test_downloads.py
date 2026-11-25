"""
Tests for bit.modules.downloads, against a local HTTP server that serves byte ranges (or, for one test, ignores them).
"""

import gzip
import hashlib
import io
import tarfile
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest # type: ignore

from bit.modules import downloads
from bit.modules.downloads import DownloadError, fetch, unpack
from bit.tests.utils import require_local_connections


class RangeHandler(SimpleHTTPRequestHandler):
    """ serves files from its directory, honouring "Range: bytes=a-b"; counts requests per path """

    ranges = True
    requests = {}

    def log_message(self, *args):
        pass

    def do_GET(self):
        path = Path(self.directory) / self.path.lstrip("/")
        RangeHandler.requests[self.path] = RangeHandler.requests.get(self.path, 0) + 1
        if not path.is_file():
            self.send_error(404)
            return
        data = path.read_bytes()
        wanted = self.headers.get("Range")
        if wanted and self.ranges:
            start, end = wanted.split("=")[1].split("-")
            start, end = int(start), min(int(end), len(data) - 1)
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
            body = data[start:end + 1]
        else:
            self.send_response(200)
            body = data
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@pytest.fixture
def server(tmp_path, monkeypatch):
    require_local_connections()
    monkeypatch.setattr(downloads, "PIECE", 1000)
    monkeypatch.setattr(downloads, "BLOCK", 256)
    monkeypatch.setattr(downloads, "QUIET_BELOW", 0)
    served = tmp_path / "served"
    served.mkdir()
    RangeHandler.requests = {}
    RangeHandler.ranges = True
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), partial(RangeHandler, directory=str(served)))
    # (shutdown() waits for the serving loop's next look, half a second apart unless told otherwise)
    threading.Thread(target=httpd.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True).start()
    yield served, f"http://127.0.0.1:{httpd.server_address[1]}/"
    httpd.shutdown()


def write_file(path, size, seed):
    data = bytes((i * seed + seed) % 251 for i in range(size))
    path.write_bytes(data)
    return hashlib.md5(data).hexdigest()


def test_several_files_share_the_connections(server, tmp_path):
    served, base = server
    sums = {name: write_file(served / name, size, seed)
            for name, size, seed in (("a.part-000", 4500, 3), ("a.part-001", 1000, 5), ("a.part-002", 7, 7))}
    out = tmp_path / "out"

    notes = fetch([(base + name, str(out / name), md5) for name, md5 in sums.items()], "test", connections=3)

    for name in sums:
        assert (out / name).read_bytes() == (served / name).read_bytes()
    assert sorted(p.name for p in out.iterdir()) == sorted(sums)        # no .partial or .done files left
    assert RangeHandler.requests["/a.part-000"] == 1 + 5                # the size lookup, then 5 pieces of 1000
    assert "3 file(s) match their md5" in notes


def test_resumes_from_the_finished_pieces(server, tmp_path):
    served, base = server
    md5 = write_file(served / "big.tar.gz", 5000, 11)
    data = (served / "big.tar.gz").read_bytes()
    out = tmp_path / "out"
    out.mkdir()
    # an earlier try finished pieces 0 and 3 (bytes 0-999 and 3000-3999)
    partial_file = bytearray(5000)
    partial_file[0:1000], partial_file[3000:4000] = data[0:1000], data[3000:4000]
    (out / "big.tar.gz.partial").write_bytes(bytes(partial_file))
    (out / "big.tar.gz.partial.done").write_text("0\n3\n")

    notes = fetch([(base + "big.tar.gz", str(out / "big.tar.gz"), md5)], "test")

    assert (out / "big.tar.gz").read_bytes() == data
    assert RangeHandler.requests["/big.tar.gz"] == 1 + 3                # only pieces 1, 2 and 4 fetched
    assert any("already there from an earlier try" in note for note in notes)


def test_every_piece_already_there_is_moved_into_place(server, tmp_path):
    served, base = server
    md5 = write_file(served / "big.tar.gz", 3000, 17)
    out = tmp_path / "out"
    out.mkdir()
    # an earlier try recorded all three pieces but stopped before renaming the file
    (out / "big.tar.gz.partial").write_bytes((served / "big.tar.gz").read_bytes())
    (out / "big.tar.gz.partial.done").write_text("0\n1\n2\n")

    fetch([(base + "big.tar.gz", str(out / "big.tar.gz"), md5)], "test")

    assert (out / "big.tar.gz").read_bytes() == (served / "big.tar.gz").read_bytes()
    assert sorted(p.name for p in out.iterdir()) == ["big.tar.gz"]
    assert RangeHandler.requests["/big.tar.gz"] == 1                    # just the size lookup


def test_a_complete_file_is_only_checked(server, tmp_path):
    served, base = server
    md5 = write_file(served / "f.bin", 2000, 13)
    (tmp_path / "f.bin").write_bytes((served / "f.bin").read_bytes())
    fetch([(base + "f.bin", str(tmp_path / "f.bin"), md5)], "test")
    assert RangeHandler.requests == {}


def test_md5_mismatch_removes_the_file_and_fails(server, tmp_path):
    served, base = server
    write_file(served / "f.bin", 2000, 17)
    with pytest.raises(DownloadError, match="didn't match their md5"):
        fetch([(base + "f.bin", str(tmp_path / "f.bin"), "0" * 32)], "test")
    assert not (tmp_path / "f.bin").exists()


def test_a_server_without_ranges_and_no_md5(server, tmp_path):
    served, base = server
    RangeHandler.ranges = False
    write_file(served / "content", 3500, 19)
    notes = fetch([(base + "content", str(tmp_path / "db.tar.gz"), None)], "test")
    assert (tmp_path / "db.tar.gz").read_bytes() == (served / "content").read_bytes()
    assert RangeHandler.requests["/content"] == 2                       # the size lookup, then the whole file at once
    assert "no md5 published for db.tar.gz, so only the size was checked" in notes


def test_missing_file_fails(server, tmp_path):
    _, base = server
    with pytest.raises(DownloadError, match="failed"):
        fetch([(base + "nothing-here", str(tmp_path / "x"), None)], "test")


def test_unpack_feeds_the_parts_through_the_decompressor_into_tar(tmp_path, capsys):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as tar:
        for name, text in (("top/a.txt", b"alpha"), ("top/sub/b.txt", b"beta" * 1000)):
            info = tarfile.TarInfo(name)
            info.size = len(text)
            tar.addfile(info, io.BytesIO(text))
    packed = gzip.compress(buffer.getvalue())
    parts = []
    for i in range(3):
        parts.append(tmp_path / f"x.tar.gz.part-{i}")
        parts[-1].write_bytes(packed[i * len(packed) // 3:(i + 1) * len(packed) // 3])

    unpack(parts, tmp_path / "out", "Unpacking test", ["gzip", "-dc"], strip_components=1)

    assert (tmp_path / "out" / "a.txt").read_bytes() == b"alpha"
    assert (tmp_path / "out" / "sub" / "b.txt").read_bytes() == b"beta" * 1000


def test_unpack_failure_says_so(tmp_path):
    bad = tmp_path / "bad.tar.gz"
    bad.write_bytes(b"not a gzip file")
    with pytest.raises(DownloadError, match="unpacking bad.tar.gz failed"):
        unpack([bad], tmp_path / "out", "Unpacking test", ["gzip", "-dc"])


def test_progress_prints_plain_lines_when_not_a_terminal(capsys, monkeypatch):
    monkeypatch.setattr(downloads, "QUIET_BELOW", 0)
    progress = downloads.Progress("Downloading x", 1000)
    for _ in range(10):
        progress.add(100)
    progress.close()
    lines = capsys.readouterr().err.splitlines()
    assert lines[0].strip() == "Downloading x: 1000 B"
    assert lines[-1].strip().startswith("Downloading x: 100%")
    assert len(lines) == 1 + 20
