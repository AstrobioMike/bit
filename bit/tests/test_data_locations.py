import os
from unittest import mock
import pytest

import bit.modules.data_locations as mod
from bit.modules.data_locations import (
    check_and_report_env_variables,
    check_location_var_is_set_and_writable,
    get_variable_path,
    set_variable_path,
    set_env_variables,
    modify_conda_activate_startup_script,
    notify_to_reactivate_conda,
)


# ───────────────────────── get_variable_path ─────────────────────────

def test_get_variable_path_returns_value(monkeypatch):
    monkeypatch.setenv("SOME_BIT_VAR", "/data/here")
    assert get_variable_path("SOME_BIT_VAR") == "/data/here"


def test_get_variable_path_returns_false_when_unset(monkeypatch):
    monkeypatch.delenv("SOME_BIT_VAR", raising=False)
    assert get_variable_path("SOME_BIT_VAR") is False


# ───────────── check_location_var_is_set_and_writable ─────────────

def test_check_location_set_and_writable(tmp_path, monkeypatch):
    monkeypatch.setenv("BIT_TEST_DB", str(tmp_path))
    path, writable = check_location_var_is_set_and_writable("BIT_TEST_DB")
    assert path == str(tmp_path)
    assert writable is True


def test_check_location_unset_exits(monkeypatch, capsys):
    monkeypatch.delenv("BIT_TEST_DB", raising=False)
    with pytest.raises(SystemExit) as e:
        check_location_var_is_set_and_writable("BIT_TEST_DB")
    assert e.value.code == 1
    assert "does not seem to be set" in capsys.readouterr().out


def test_check_location_empty_string_exits(monkeypatch, capsys):
    monkeypatch.setenv("BIT_TEST_DB", "")
    with pytest.raises(SystemExit) as e:
        check_location_var_is_set_and_writable("BIT_TEST_DB")
    assert e.value.code == 1
    assert "does not seem to be set" in capsys.readouterr().out


def test_check_location_not_writable(tmp_path, monkeypatch):
    monkeypatch.setenv("BIT_TEST_DB", str(tmp_path))
    # force the writability probe to report False without touching real perms
    monkeypatch.setattr(mod.os, "access", lambda p, mode: False)
    path, writable = check_location_var_is_set_and_writable("BIT_TEST_DB")
    assert writable is False


# ───────────────────────── set_variable_path ─────────────────────────

def test_set_variable_path_keeps_current_on_no():
    with mock.patch("builtins.input", side_effect=["n"]):
        result = set_variable_path("TAXONKIT_DB", "/existing/path/")
    assert result == "/existing/path/"


def test_set_variable_path_change_to_new(tmp_path):
    target = tmp_path / "newdb"
    with mock.patch("builtins.input", side_effect=["y", str(target)]):
        result = set_variable_path("TAXONKIT_DB", "/old/path/")
    # function appends a trailing separator via os.path.join(path, "")
    assert result == os.path.join(str(target), "")
    assert target.is_dir()          # created since it didn't exist


def test_set_variable_path_no_current(tmp_path):
    target = tmp_path / "freshdb"
    with mock.patch("builtins.input", side_effect=[str(target)]):
        result = set_variable_path("GTDB_DIR", False)
    assert result == os.path.join(str(target), "")


def test_set_variable_path_reprompts_on_bad_yes_no(capsys):
    # first answer invalid, then "n" to keep current
    with mock.patch("builtins.input", side_effect=["maybe", "n"]):
        result = set_variable_path("GO_DB_DIR", "/cur/")
    assert result == "/cur/"
    assert "Must respond with" in capsys.readouterr().out


def test_set_variable_path_reprompts_on_not_writable(tmp_path, monkeypatch, capsys):
    good = tmp_path / "writable_db"
    bad = tmp_path / "unwritable_db"
    good.mkdir()
    bad.mkdir()
    # both dirs exist (so makedirs is skipped); os.access reports the first
    # as not writable, then the second as writable
    monkeypatch.setattr(mod.os, "access", mock.Mock(side_effect=[False, True]))
    with mock.patch("builtins.input", side_effect=[str(bad), str(good)]):
        result = set_variable_path("GO_DB_DIR", False)
    assert result == os.path.join(str(good), "")
    assert "not writable" in capsys.readouterr().out


def test_set_variable_path_reprompts_on_not_absolute(tmp_path, monkeypatch, capsys):
    good = tmp_path / "abs_db"
    good.mkdir()
    # no-op makedirs so the relative first answer doesn't create junk in cwd;
    # access always True so writability passes and isabs is what trips/clears
    monkeypatch.setattr(mod.os, "makedirs", lambda *a, **k: None)
    monkeypatch.setattr(mod.os, "access", lambda p, mode: True)
    with mock.patch("builtins.input", side_effect=["relative_path", str(good)]):
        result = set_variable_path("GTDB_DIR", False)
    assert result == os.path.join(str(good), "")
    assert "absolute path" in capsys.readouterr().out


# ───────────────────────── set_env_variables ─────────────────────────

def test_set_env_variables_iterates_all_keys():
    # patch the two helpers so we test the orchestration, not the prompts
    with mock.patch.object(mod, "get_variable_path", return_value=False) as g, \
         mock.patch.object(mod, "set_variable_path",
                           side_effect=lambda var, cur: f"/{var}/") as s:
        result = set_env_variables()

    assert set(result.keys()) == {"TAXONKIT_DB", "GTDB_DIR",
                                  "GO_DB_DIR", "NCBI_assembly_data_dir"}
    assert result["GTDB_DIR"] == "/GTDB_DIR/"
    assert g.call_count == 4
    assert s.call_count == 4


# ─────────────────── check_and_report_env_variables ───────────────────

def test_check_and_report_all_set(tmp_path, monkeypatch, capsys):
    for var in ("TAXONKIT_DB", "GTDB_DIR", "GO_DB_DIR", "NCBI_assembly_data_dir"):
        d = tmp_path / var
        d.mkdir()
        monkeypatch.setenv(var, str(d))

    check_and_report_env_variables()
    out = capsys.readouterr().out
    assert "bit environment variables are set" in out
    assert "TAXONKIT_DB" in out
    assert "NCBI_assembly_data_dir" in out


def test_check_and_report_warns_when_not_writable(tmp_path, monkeypatch, capsys):
    for var in ("TAXONKIT_DB", "GTDB_DIR", "GO_DB_DIR", "NCBI_assembly_data_dir"):
        d = tmp_path / var
        d.mkdir()
        monkeypatch.setenv(var, str(d))
    monkeypatch.setattr(mod.os, "access", lambda p, mode: False)

    check_and_report_env_variables()
    out = capsys.readouterr().out
    assert "is not writable" in out


# ───────────── modify_conda_activate_startup_script ─────────────

def test_modify_conda_script_rewrites_exports(tmp_path, monkeypatch):
    conda_prefix = tmp_path / "conda"
    activate_dir = conda_prefix / "etc/conda/activate.d"
    activate_dir.mkdir(parents=True)
    startup = activate_dir / "bit.sh"
    startup.write_text(
        "export TAXONKIT_DB=/old/tax\n"   # managed var -> dropped/replaced
        "# a comment line\n"               # unmanaged -> kept
        "export OTHER_VAR=/keep/this\n"    # unmanaged -> kept
    )

    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))

    paths_dict = {
        "TAXONKIT_DB": "/new/tax/",
        "GTDB_DIR": "/new/gtdb/",
        "GO_DB_DIR": "/new/go/",
        "NCBI_assembly_data_dir": "/new/ncbi/",
    }

    modify_conda_activate_startup_script(paths_dict)

    contents = startup.read_text()
    assert "export TAXONKIT_DB=/old/tax" not in contents
    assert "export TAXONKIT_DB=/new/tax/" in contents
    assert "# a comment line" in contents
    assert "export OTHER_VAR=/keep/this" in contents
    for var in paths_dict:
        assert f"export {var}=" in contents


def test_modify_conda_script_no_write_access_writes_to_home(tmp_path, monkeypatch):
    # conda startup exists and is readable, but reported as NOT writable, so the
    # function should write into the user's ~/.config/bit/ instead, and should
    # strip the if/./fi sourcing block while it's at it.
    conda_prefix = tmp_path / "conda"
    activate_dir = conda_prefix / "etc/conda/activate.d"
    activate_dir.mkdir(parents=True)
    startup = activate_dir / "bit.sh"
    original_text = (
        "if [ -f ~/.config/bit/bit.sh ]; then\n"
        ". ~/.config/bit/bit.sh\n"
        "fi\n"
        "export TAXONKIT_DB=/old/tax\n"
        "# a comment line\n"
        "export OTHER_VAR=/keep/this\n"
    )
    startup.write_text(original_text)

    fake_home = tmp_path / "home"
    fake_home.mkdir()

    monkeypatch.setenv("CONDA_PREFIX", str(conda_prefix))
    monkeypatch.setattr(mod.os.path, "expanduser", lambda p: str(fake_home))
    # report the conda startup script as not writable -> home-dir branch
    monkeypatch.setattr(mod.os, "access", lambda p, mode: False)

    paths_dict = {
        "TAXONKIT_DB": "/new/tax/",
        "GTDB_DIR": "/new/gtdb/",
        "GO_DB_DIR": "/new/go/",
        "NCBI_assembly_data_dir": "/new/ncbi/",
    }

    modify_conda_activate_startup_script(paths_dict)

    user_script = fake_home / ".config/bit/bit.sh"
    assert user_script.is_file()
    contents = user_script.read_text()

    # sourcing block (if/./fi) stripped in the home-dir branch
    assert "if [ -f" not in contents
    assert ". ~/.config/bit/bit.sh" not in contents
    # managed var replaced, unmanaged lines kept
    assert "export TAXONKIT_DB=/old/tax" not in contents
    assert "# a comment line" in contents
    assert "export OTHER_VAR=/keep/this" in contents
    for var in paths_dict:
        assert f"export {var}=" in contents

    # the conda-area script is left untouched in this branch
    assert startup.read_text() == original_text


# ───────────────────── notify_to_reactivate_conda ─────────────────────

def test_notify_to_reactivate_conda(monkeypatch, capsys):
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "bit-env")
    notify_to_reactivate_conda()
    out = capsys.readouterr().out
    assert "conda activate bit-env" in out
    assert "reactivate the conda environment" in out
