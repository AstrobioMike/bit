"""
Downloads of large files with one progress bar for the lot (used to set up the reference databases, see
bit.modules.ref_dbs).

fetch() downloads several files, or one big one, over several connections at once: every file is fetched in 64 MB
pieces, so a single big file gets them all and several files share them. An interrupted download picks up where it
stopped (the finished pieces of <file> are listed in <file>.partial.done), and each file is checked against its md5
when one is given. unpack() feeds files through a decompressor into tar, with a bar of how much has been read.

The bars are tqdm's on a terminal; otherwise a plain line is printed at every 5%.
"""

import hashlib
import os
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm

PIECE = 64 * 1024 * 1024
BLOCK = 1024 * 1024
QUIET_BELOW = 20 * 1024 * 1024    # no bar for less than this (e.g. an md5 list)
TRIES = 6
TIMEOUT = 60
USER_AGENT = "bit"


class DownloadError(Exception):
    pass


class Progress:
    """ a tqdm bar on a terminal, otherwise a plain line at every 5%; thread-safe """

    def __init__(self, label, total, indent="    "):
        self.label, self.total, self.done = label, total, 0
        self.indent = indent
        self.lock = threading.Lock()
        self.quiet = total < QUIET_BELOW
        self.tty = sys.stderr.isatty()
        self.next_step = 5
        self.start = time.time()
        self.bar = None
        if self.quiet:
            return
        if self.tty:
            self.bar = tqdm(total=total, desc=f"{indent}{label}", unit="B", unit_scale=True, unit_divisor=1024,
                            ncols=90, file=sys.stderr, miniters=1, mininterval=0.5)
        else:
            print(f"{indent}{label}: {human(total)}", file=sys.stderr, flush=True)

    def add(self, n):
        with self.lock:
            self.done += n
            if self.bar is not None:
                self.bar.update(n)
            elif not self.quiet and self.total:
                while self.next_step <= 100 and self.done * 100 >= self.next_step * self.total:
                    print(f"{self.indent}{self.label}: {self.next_step}% ({human(self.done)} of {human(self.total)})",
                          file=sys.stderr, flush=True)
                    self.next_step += 5

    def close(self):
        if self.bar is not None:
            self.bar.close()
        return time.time() - self.start


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def _request(url, headers=None):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return urllib.request.urlopen(req, timeout=TIMEOUT)


def _retrying(what, action):
    """ runs action(), retrying network errors (not a 404) with a growing wait """

    for attempt in range(1, TRIES + 1):
        try:
            return action()
        except urllib.error.HTTPError as e:
            if e.code == 404 or attempt == TRIES:
                raise DownloadError(f"{what} failed: {e}")
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            if attempt == TRIES:
                raise DownloadError(f"{what} failed: {e}")
        time.sleep(min(60, 5 * attempt))


def read_url(url):
    """ the contents of a small file (e.g. an md5 list or a record of a release) """

    def get():
        with _request(url) as resp:
            return resp.read()

    return _retrying(f"Getting {url}", get)


def _remote_size(url):
    """ (the file's size, whether the server sends byte ranges), from a request for its first byte """

    def ask():
        with _request(url, {"Range": "bytes=0-0"}) as resp:
            if resp.status == 206:
                return int(resp.headers["Content-Range"].rsplit("/", 1)[1]), True
            length = resp.headers.get("Content-Length")
            if length is None:
                raise DownloadError(f"{url} doesn't say how big it is")
            return int(length), False

    return _retrying(f"Looking up {url}", ask)


class _Target:
    """ one file being downloaded: <path>.partial, with the finished pieces listed in <path>.partial.done """

    def __init__(self, url, path):
        self.url, self.path = url, path
        self.partial, self.done_list = path + ".partial", path + ".partial.done"
        self.lock = threading.Lock()
        self.fd = None

    def prepare(self, size, ranges):
        self.size, self.ranges = size, ranges
        self.pieces = max(1, -(-size // PIECE)) if ranges else 1
        self.finished = set()
        if ranges and os.path.exists(self.done_list) and os.path.exists(self.partial) \
                and os.path.getsize(self.partial) == size:
            with open(self.done_list) as f:
                self.finished = {int(line) for line in f if line.strip().isdigit()}
        else:
            with open(self.partial, "wb") as f:
                f.truncate(size)
            open(self.done_list, "w").close()
        self.fd = os.open(self.partial, os.O_WRONLY)
        # (every piece already there: an earlier try stopped just before moving the file into place)
        self._finish_if_complete()

    def piece_range(self, i):
        if not self.ranges:
            return 0, self.size - 1
        start = i * PIECE
        return start, min(self.size, start + PIECE) - 1

    def bytes_done(self):
        return sum(end - start + 1 for start, end in map(self.piece_range, self.finished))

    def piece_done(self, i):
        with self.lock:
            self.finished.add(i)
            with open(self.done_list, "a") as f:
                f.write(f"{i}\n")
            self._finish_if_complete()

    def _finish_if_complete(self):
        if len(self.finished) == self.pieces:
            os.close(self.fd)
            self.fd = None
            os.replace(self.partial, self.path)
            os.remove(self.done_list)


def _fetch_piece(target, i, progress):
    start, end = target.piece_range(i)

    def attempt():
        written = 0
        try:
            headers = {"Range": f"bytes={start}-{end}"} if target.ranges else {}
            with _request(target.url, headers) as resp:
                if target.ranges and resp.status != 206:
                    raise OSError("the server ignored the byte range")
                offset = start
                while block := resp.read(BLOCK):
                    os.pwrite(target.fd, block, offset)
                    offset += len(block)
                    written += len(block)
                    progress.add(len(block))
            if written != end - start + 1:
                raise OSError(f"got {written} of {end - start + 1} bytes")
        except BaseException:
            progress.add(-written)
            raise

    _retrying(f"Downloading {os.path.basename(target.path)}", attempt)
    target.piece_done(i)


def _md5_of(path, progress):
    digest = hashlib.md5()
    with open(path, "rb") as f:
        while block := f.read(8 * BLOCK):
            digest.update(block)
            progress.add(len(block))
    return digest.hexdigest()


def fetch(downloads, label, connections=8, indent="    "):
    """
    downloads:   [(url, destination path, md5 or None), ...]
    label:       what the progress bars call them
    connections: how many pieces are fetched at once
    Returns a list of lines saying what was done (for a log). A file that is already there is only checked. A file
    that doesn't match its md5 is removed and DownloadError is raised (fetching again downloads it again).
    """

    notes = []
    targets = [_Target(url, path) for url, path, _ in downloads]
    for target in targets:
        os.makedirs(os.path.dirname(os.path.abspath(target.path)), exist_ok=True)
    to_fetch = [t for t in targets if not os.path.exists(t.path)]

    with ThreadPoolExecutor(max_workers=min(8, len(to_fetch) or 1)) as pool:
        looked_up = list(pool.map(lambda t: _remote_size(t.url), to_fetch))
    for target, (size, ranges) in zip(to_fetch, looked_up):
        if size == 0:
            open(target.path, "wb").close()
        else:
            target.prepare(size, ranges)
    to_fetch = [t for t in to_fetch if not os.path.exists(t.path)]

    if to_fetch:
        total = sum(t.size for t in to_fetch)
        already = sum(t.bytes_done() for t in to_fetch)
        notes.append(f"downloading {len(to_fetch)} file(s), {human(total)}" +
                     (f" ({human(already)} of it already there from an earlier try)" if already else ""))
        progress = Progress(f"Downloading {label}", total, indent)
        progress.add(already)
        jobs = [(t, i) for t in to_fetch for i in range(t.pieces) if i not in t.finished]
        pool = ThreadPoolExecutor(max_workers=connections)
        try:
            for future in as_completed([pool.submit(_fetch_piece, t, i, progress) for t, i in jobs]):
                future.result()
        except BaseException:
            pool.shutdown(wait=True, cancel_futures=True)
            progress.close()
            for target in to_fetch:
                if target.fd is not None:
                    os.close(target.fd)
            raise
        pool.shutdown()
        seconds = progress.close()
        took = f"{seconds:.0f} s" if seconds < 120 else f"{seconds / 60:.1f} min"
        notes.append(f"  done in {took} ({human((total - already) / max(seconds, 1e-9))}/s)")

    checked = [(t, md5) for t, (_, _, md5) in zip(targets, downloads) if md5]
    if checked:
        progress = Progress(f"Checking md5 of {label}", sum(os.path.getsize(t.path) for t, _ in checked), indent)
        with ThreadPoolExecutor(max_workers=min(4, len(checked))) as pool:
            sums = list(pool.map(lambda c: _md5_of(c[0].path, progress), checked))
        progress.close()
        bad = [t.path for (t, md5), got in zip(checked, sums) if got != md5.lower()]
        for path in bad:
            os.remove(path)
        if bad:
            names = ", ".join(os.path.basename(p) for p in bad)
            raise DownloadError(f"{len(bad)} download(s) didn't match their md5 and were removed: {names} "
                                "(running the same command again downloads them again)")
        notes.append(f"{len(checked)} file(s) match their md5")
    unchecked = [os.path.basename(t.path) for t, (_, _, md5) in zip(targets, downloads) if not md5]
    if unchecked:
        notes.append(f"no md5 published for {', '.join(unchecked)}, so only the size was checked")

    return notes


def unpack(paths, into, label, decompress, strip_components=0, indent="    "):
    """
    Feeds the files (the parts of one archive, or a single one) through `decompress` (e.g. ["zstd", "-dc"]) into
    tar, which unpacks into `into`, with a bar of how much has been read.
    """

    os.makedirs(into, exist_ok=True)
    total = sum(os.path.getsize(p) for p in paths)
    tar_cmd = ["tar", "-xf", "-", "-C", str(into)] + ([f"--strip-components={strip_components}"] if strip_components else [])
    # (tar's messages go to a file, not a pipe, so a lot of them can't stall it while the parts are still being fed)
    tar_messages = tempfile.TemporaryFile()
    unzip = subprocess.Popen(decompress, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    tar = subprocess.Popen(tar_cmd, stdin=unzip.stdout, stderr=tar_messages)
    unzip.stdout.close()
    progress = Progress(label, total, indent)
    try:
        for path in paths:
            with open(path, "rb") as f:
                while block := f.read(8 * BLOCK):
                    unzip.stdin.write(block)
                    progress.add(len(block))
        unzip.stdin.close()
    except BrokenPipeError:
        pass
    finally:
        progress.close()
    unzip_code, tar_code = unzip.wait(), tar.wait()
    tar_messages.seek(0)
    tar_err = tar_messages.read()
    tar_messages.close()
    if unzip_code or tar_code:
        raise DownloadError(f"unpacking {', '.join(os.path.basename(p) for p in paths)} failed: "
                            f"{tar_err.decode(errors='replace').strip() or f'{decompress[0]} exited with {unzip_code}'}")
