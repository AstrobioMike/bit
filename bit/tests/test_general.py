import os
import sys
import gzip
import socket
import urllib.error
from unittest import mock
import pytest # type: ignore
import bit.modules.general as general_mod
from bit.modules.general import (tee, report_failure, color_text, wprint,
                                 report_message, check_files_are_found,
                                 check_if_output_dir_exists, notify_premature_exit,
                                 is_gzipped, sniff_delimiter, colnames,
                                 log_command_run, report_version,
                                 download_with_tqdm, _TooSlow)


def test_tee_removes_ansi_and_logs(tmp_path):

    log_path = tmp_path / "test.log"
    msg_with_color = "\033[1mBold text\033[0m and normal"

    with mock.patch("sys.stdout", new_callable=mock.MagicMock()) as mock_stdout:
        tee(msg_with_color, log_path)

    calls = [mock.call(msg_with_color), mock.call("\n")]
    mock_stdout.write.assert_has_calls(calls)

    logged = log_path.read_text()
    assert "Bold text and normal" in logged
    assert "\033" not in logged

    log_path.unlink()
    tee("No newline", log_path, end="")
    assert log_path.read_text() == "No newline"


def test_report_failure_output_and_exit(capsys):

    with pytest.raises(SystemExit) as e:
        report_failure("Test error message", color="red")

    assert e.value.code == 1
    captured = capsys.readouterr()
    assert "Test error message" in captured.out
    assert "Exiting for now :(" in captured.out


def test_color_text_wraps_when_tty(monkeypatch):

    monkeypatch.setattr(sys.stdout, "isatty", lambda: True)
    out = color_text("hello", "red")
    assert "hello" in out
    assert out.startswith("\033")          # ANSI prefix applied
    assert out.endswith("\033[0m")


def test_color_text_plain_when_not_tty(monkeypatch):

    monkeypatch.setattr(sys.stdout, "isatty", lambda: False)
    assert color_text("hello", "red") == "hello"


def test_color_text_none_color_is_plain(monkeypatch):

    monkeypatch.setattr(sys.stdout, "isatty", lambda: True)
    assert color_text("hello", "none") == "hello"


def test_wprint_indents_and_wraps(capsys):

    wprint("word " * 40, width=40)
    out = capsys.readouterr().out
    assert out.startswith("  ")                       # initial indent
    assert all(len(line) <= 40 for line in out.splitlines())


def test_report_message_collapses_whitespace(capsys):

    report_message("a    b\n   c", color="none", leading_newline=False)
    out = capsys.readouterr().out
    assert "a b c" in out


def test_report_message_leading_and_trailing_newlines(capsys):

    report_message("hi", color="none", leading_newline=True, trailing_newline=True)
    out = capsys.readouterr().out
    assert out.startswith("\n")
    assert out.endswith("\n\n")


def test_check_files_are_found_passes_for_existing(tmp_path):

    f = tmp_path / "exists.txt"
    f.write_text("x")
    # should not raise / exit
    check_files_are_found([str(f)])


def test_check_files_are_found_exits_on_missing(tmp_path, capsys):

    missing = tmp_path / "nope.txt"
    with pytest.raises(SystemExit) as e:
        check_files_are_found([str(missing)])
    assert e.value.code == 1
    assert "not able to find" in capsys.readouterr().out


def test_check_if_output_dir_exists_no_force_exits(tmp_path, capsys):

    d = tmp_path / "outdir"
    d.mkdir()
    with pytest.raises(SystemExit) as e:
        check_if_output_dir_exists(str(d), force_overwrite=False)
    assert e.value.code == 1
    assert "already exists" in capsys.readouterr().out


def test_check_if_output_dir_exists_force_removes(tmp_path):

    d = tmp_path / "outdir"
    d.mkdir()
    (d / "sentinel.txt").write_text("x")
    check_if_output_dir_exists(str(d), force_overwrite=True)
    assert not d.exists()


def test_check_if_output_dir_exists_absent_is_noop(tmp_path):

    d = tmp_path / "not_there"
    check_if_output_dir_exists(str(d), force_overwrite=False)   # no raise
    assert not d.exists()


def test_notify_premature_exit(capsys):

    with pytest.raises(SystemExit) as e:
        notify_premature_exit()
    assert e.value.code == 1
    assert "Exiting for now" in capsys.readouterr().err


def test_is_gzipped_true(tmp_path):

    f = tmp_path / "f.gz"
    with gzip.open(f, "wt") as g:
        g.write("hello")
    assert is_gzipped(str(f)) is True


def test_is_gzipped_false(tmp_path):

    f = tmp_path / "f.txt"
    f.write_text("hello")
    assert is_gzipped(str(f)) is False


@pytest.mark.parametrize("line,expected", [
    ("a\tb\tc", "\t"),
    ("a,b,c", ","),
    ("a|b|c", "|"),
])
def test_sniff_delimiter_detects(line, expected):

    assert sniff_delimiter(line) == expected


def test_sniff_delimiter_returns_none_on_failure():

    # a single token with no delimiter candidates
    assert sniff_delimiter("singletoken") is None


def test_colnames_from_path(tmp_path, capsys):

    f = tmp_path / "table.tsv"
    f.write_text("colA\tcolB\tcolC\n1\t2\t3\n")

    class Args:
        input_file = str(f)

    colnames(Args())
    out = capsys.readouterr().out
    assert "colA" in out and "colB" in out and "colC" in out
    assert "1  colA" in out          # 1-indexed numbering


def test_colnames_empty_file_exits(tmp_path):

    f = tmp_path / "empty.tsv"
    f.write_text("")

    class Args:
        input_file = str(f)

    with pytest.raises(SystemExit) as e:
        colnames(Args())
    assert e.value.code == 1


def test_colnames_from_stream(capsys):

    import io

    class Args:
        input_file = io.StringIO("x\ty\n1\t2\n")

    colnames(Args())
    out = capsys.readouterr().out
    assert "x" in out and "y" in out


def test_log_command_run_writes_file(tmp_path):

    log_command_run("bit some-cmd --flag", str(tmp_path))
    log_file = tmp_path / "command-execution-info.txt"
    assert log_file.is_file()
    contents = log_file.read_text()
    assert "bit version: v" in contents
    assert "bit some-cmd --flag" in contents


def test_log_command_run_custom_path(tmp_path):

    custom = tmp_path / "mylog.txt"
    log_command_run("bit x", str(tmp_path), log_file=str(custom))
    assert custom.is_file()
    assert "bit x" in custom.read_text()


def test_report_version_runs(capsys):

    report_version()
    out = capsys.readouterr().out
    assert "bit" in out
    assert "github.com/AstrobioMike/bit" in out
    assert "F1000Research" in out


def test_tee_removes_ansi_and_logs(tmp_path):

    log_path = tmp_path / "test.log"
    msg_with_color = "\033[1mBold text\033[0m and normal"

    with mock.patch("sys.stdout", new_callable=mock.MagicMock()) as mock_stdout:
        tee(msg_with_color, log_path)

    # stdout write calls come as two separate writes: msg, then end
    calls = [mock.call(msg_with_color), mock.call("\n")]
    mock_stdout.write.assert_has_calls(calls)

    logged = log_path.read_text()
    assert "Bold text and normal" in logged
    assert "\033" not in logged

    log_path.unlink()
    tee("No newline", log_path, end="")
    assert log_path.read_text() == "No newline"


def test_report_failure_output_and_exit(capsys):

    with pytest.raises(SystemExit) as e:
        report_failure("Test error message", color="red")

    assert e.value.code == 1

    captured = capsys.readouterr()

    assert "Test error message" in captured.out
    assert "Exiting for now :(" in captured.out


################################################################################
# download_with_tqdm (retries + speed-gate now live inside this function)
################################################################################

def _fake_stream(fail_times, exc=None, dest_content=b"ok"):
    """
    _stream_once side effect: raise `exc` on the first `fail_times` calls, then
    "download" by writing dest_content to `filename`. fail_times=None fails every
    call. Signature mirrors _stream_once(url, filename, desc, total, leave,
    floor_bytes_per_s, probe_seconds).
    """
    if exc is None:
        exc = socket.timeout("handshake timed out")
    seen = {"count": 0}

    def _stream(url, filename, desc, total, leave, floor_bytes_per_s, probe_seconds):
        if fail_times is None or seen["count"] < fail_times:
            seen["count"] += 1
            raise exc
        with open(filename, "wb") as out:
            out.write(dest_content)
    return _stream


def test_download_with_tqdm_returns_path_on_success(tmp_path):
    """A clean download returns the destination path (callers feed it to pandas)."""
    dest = str(tmp_path / "out.bin")
    with mock.patch.object(general_mod, "_stream_once",
                           side_effect=_fake_stream(fail_times=0)), \
         mock.patch("bit.modules.general.urllib.request.urlopen"):
        result = download_with_tqdm("https://example.test/x", "label", dest)
    assert result == dest
    assert os.path.exists(dest)


def test_download_with_tqdm_urlopen_returns_response():
    """urlopen=True short-circuits and returns the response object."""
    sentinel = object()
    with mock.patch("bit.modules.general.urllib.request.urlopen",
                    return_value=sentinel), \
         mock.patch("bit.modules.general.urllib.request.build_opener"), \
         mock.patch("bit.modules.general.urllib.request.install_opener"):
        result = download_with_tqdm("https://example.test/x", "label",
                                    urlopen=True)
    assert result is sentinel


def test_download_with_tqdm_retries_then_succeeds(tmp_path):
    """A transient failure is retried the expected number of times, sleeping
    between attempts, then succeeds and returns the path."""
    dest = str(tmp_path / "out.bin")
    with mock.patch.object(general_mod, "_stream_once",
                           side_effect=_fake_stream(fail_times=2)), \
         mock.patch("bit.modules.general.urllib.request.urlopen"), \
         mock.patch("bit.modules.general.time.sleep") as mock_sleep:
        result = download_with_tqdm("https://example.test/x", "label", dest)
    assert result == dest
    assert mock_sleep.call_count == 2   # slept between the 2 failures, not after success


def test_download_with_tqdm_does_not_retry_404(tmp_path):
    """A 404 is not transient and must not be retried; it raises at once and
    sleep is never called."""
    dest = str(tmp_path / "out.bin")
    err = urllib.error.HTTPError("https://example.test/x", 404, "Not Found",
                                 hdrs=None, fp=None)
    with mock.patch.object(general_mod, "_stream_once",
                           side_effect=_fake_stream(fail_times=None, exc=err)), \
         mock.patch("bit.modules.general.urllib.request.urlopen"), \
         mock.patch("bit.modules.general.time.sleep") as mock_sleep:
        with pytest.raises(urllib.error.HTTPError):
            download_with_tqdm("https://example.test/x", "label", dest,
                               attempts=4)
    mock_sleep.assert_not_called()


def test_download_with_tqdm_exhausts_retries_and_raises(tmp_path):
    """When every attempt fails transiently, the last error propagates after
    exhausting the attempts budget."""
    dest = str(tmp_path / "out.bin")
    with mock.patch.object(general_mod, "_stream_once",
                           side_effect=_fake_stream(fail_times=None)), \
         mock.patch("bit.modules.general.urllib.request.urlopen"), \
         mock.patch("bit.modules.general.time.sleep"):
        with pytest.raises(socket.timeout):
            download_with_tqdm("https://example.test/x", "label", dest,
                               attempts=3)


def test_download_with_tqdm_no_retries_is_single_shot(tmp_path):
    """retries=False collapses to a single attempt: a transient error is raised
    immediately without sleeping."""
    dest = str(tmp_path / "out.bin")
    with mock.patch.object(general_mod, "_stream_once",
                           side_effect=_fake_stream(fail_times=None)), \
         mock.patch("bit.modules.general.urllib.request.urlopen"), \
         mock.patch("bit.modules.general.time.sleep") as mock_sleep:
        with pytest.raises(socket.timeout):
            download_with_tqdm("https://example.test/x", "label", dest,
                               retries=False)
    mock_sleep.assert_not_called()


def test_download_with_tqdm_speed_gate_rerolls_then_completes(tmp_path):
    """With speed_gate on, a _TooSlow on a non-final attempt rerolls (no sleep)
    and a subsequent good attempt completes."""
    dest = str(tmp_path / "out.bin")
    calls = {"n": 0}

    def _stream(url, filename, desc, total, leave, floor_bytes_per_s, probe_seconds):
        calls["n"] += 1
        if calls["n"] == 1:
            # first attempt is gated (floor enforced) and judged too slow
            assert floor_bytes_per_s > 0
            raise _TooSlow(0.5)
        with open(filename, "wb") as out:
            out.write(b"ok")

    with mock.patch.object(general_mod, "_stream_once", side_effect=_stream), \
         mock.patch("bit.modules.general.urllib.request.urlopen"), \
         mock.patch("bit.modules.general.time.sleep") as mock_sleep:
        result = download_with_tqdm("https://example.test/x", "label", dest,
                                    speed_gate=True)
    assert result == dest
    assert calls["n"] == 2
    mock_sleep.assert_not_called()   # slow-route reroll waits no time
