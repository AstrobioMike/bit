import subprocess
import pytest


def run_cli(cmd, **kwargs):

    res = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if res.returncode != 0:
        pytest.fail(
            f"CLI failed with exit code {res.returncode}\n"
            f"STDOUT:\n{res.stdout}\n"
            f"STDERR:\n{res.stderr}"
        )
