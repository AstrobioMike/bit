import importlib
import pytest # type: ignore
from bit.cli.bit import SUBCOMMAND_MAP


@pytest.mark.parametrize("module_path", sorted(set(SUBCOMMAND_MAP.values())))
def test_subcommand_module_imports(module_path):
    importlib.import_module(module_path)
