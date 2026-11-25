import os

import pytest


def pytest_addoption(parser):
    parser.addoption("--all", action="store_true",
                     help="also run the tests marked full (e.g., a real assembly, which first builds its conda "
                          "environment if it isn't there yet); the CI runs use this")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--all"):
        return
    skip = pytest.mark.skip(reason="a full run (e.g., builds a conda environment the first time); add --all to run it")
    for item in items:
        if item.get_closest_marker("full"):
            item.add_marker(skip)


def pytest_configure(config):
    """
    The programs the tests start (the bit commands, snakemake) are measured for coverage too, through coverage's
    COVERAGE_PROCESS_START, but only when coverage is being collected: being measured makes each of them slower to
    start (a snakemake dry run takes twice as long), so with --no-cov they run as they are.

    The paths are full ones, for the programs started in another directory to find the settings too, and to leave
    what they measured with the rest rather than in that directory.
    """
    config.addinivalue_line("markers", "full: runs only with --all (e.g., a real assembly, which first builds its "
                                       "conda environment if it isn't there yet)")
    cov = config.pluginmanager.get_plugin("_cov")
    if cov is not None and not getattr(cov, "_disabled", False):
        os.environ["COVERAGE_PROCESS_START"] = os.path.join(str(config.rootpath), "pyproject.toml")
        os.environ.setdefault("COVERAGE_FILE", os.path.join(str(config.rootpath), ".coverage"))


@pytest.fixture(autouse=True)
def sites_answer(monkeypatch):
    """
    The tests don't use the network: the check made before a download (bit.modules.connectivity) finds every site
    answering, unless a test sets up otherwise.
    """
    monkeypatch.setattr("bit.modules.connectivity.probe", lambda url, timeout=None: True)
