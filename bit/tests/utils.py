import subprocess
import pytest


def run_cli(cmd, **kwargs):

    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        pytest.fail(
            f"CLI failed with exit code {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    return result


def run_cli_here(cmd, cwd=None, check=True):
    """
    Runs a `bit` command in this process, through the same entry point (bit.cli.bit.main), and returns what
    subprocess.run(cmd, capture_output=True, text=True) would: the exit code and what was written to stdout and stderr.
    This saves starting python and importing bit for each run, which is most of what a small run takes. With check
    (as run_cli), a run that doesn't exit with 0 fails the test.

    Each command keeps at least one test run through run_cli, as a program of its own.
    """
    import contextlib
    import io
    import os
    import sys
    from bit.cli.bit import main

    assert cmd[0] == "bit"
    stdout, stderr = io.StringIO(), io.StringIO()
    argv, started_in, environ = sys.argv, os.getcwd(), dict(os.environ)
    sys.argv = [str(arg) for arg in cmd]
    returncode = 0
    try:
        if cwd is not None:
            os.chdir(cwd)
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                main()
            except SystemExit as e:
                if isinstance(e.code, int) or e.code is None:
                    returncode = e.code or 0
                else:
                    print(e.code, file=sys.stderr)
                    returncode = 1
    finally:
        sys.argv = argv
        os.chdir(started_in)
        os.environ.clear()
        os.environ.update(environ)

    result = subprocess.CompletedProcess(cmd, returncode, stdout.getvalue(), stderr.getvalue())
    if check and returncode != 0:
        pytest.fail(
            f"CLI failed with exit code {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    return result


_LOCAL_CONNECTIONS_WORK = None


def require_local_connections():
    """
    Skips the test calling it when connections to a server on this machine (127.0.0.1) don't get through, as where
    security software silently drops them: the tests that serve files locally would otherwise each wait out their
    timeouts. Checked once per test session, in well under a second when they work.
    """
    global _LOCAL_CONNECTIONS_WORK
    if _LOCAL_CONNECTIONS_WORK is None:
        import socket
        import threading
        listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            listener.bind(("127.0.0.1", 0))
            listener.listen(1)
            listener.settimeout(2)

            def answer():
                try:
                    conn, _ = listener.accept()
                    with conn:
                        conn.sendall(b"ok")
                except OSError:
                    pass

            threading.Thread(target=answer, daemon=True).start()
            with socket.create_connection(listener.getsockname(), timeout=2) as client:
                client.settimeout(2)
                _LOCAL_CONNECTIONS_WORK = client.recv(2) == b"ok"
        except OSError:
            _LOCAL_CONNECTIONS_WORK = False
        finally:
            listener.close()
    if not _LOCAL_CONNECTIONS_WORK:
        pytest.skip("connections to a local server (127.0.0.1) don't get through on this machine (e.g., blocked by "
                    "security software), and this test serves files locally")
