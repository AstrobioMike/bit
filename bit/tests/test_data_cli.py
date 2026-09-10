"""
Shape tests for `bit data`.

Mirrored by gtotree/tests/cli/test_data.py in the GToTree repo -- if something
changes here, check there too.

These exist because `gtt data get` had drifted into a `source` positional with
`choices` instead of a real subparser level. It parsed the same and its help menu
read fine, but bit.py's _suppress_help_version_on_group_parsers() only recognizes a
group via argparse._SubParsersAction, so TAB after `gtt data get` padded the database
names with -f/-h/-v while `bit data get` offered just the names. bit had the right
shape here already; these lock it in.
"""

import argparse
import io
import contextlib

import pytest # type: ignore

from bit.cli import data as data_cli
from bit.cli.bit import _suppress_help_version_on_group_parsers


HELP_AND_VERSION = {"-h", "--help", "-v", "--version"}

# every `bit data get` source except test-data, which is a group of its own rather
# than a leaf taking the common -q/-f flags
GET_SOURCES = ["go-dbs", "gtdb-data", "ncbi-assembly-data", "ncbi-tax-data"]


def _subparsers_action(parser):
    """The parser's subparsers action, or None if it isn't a group node."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None


def _node(*path):
    """Walk down from the `data` parser by subcommand name."""
    node = data_cli.build_parser()
    for name in path:
        action = _subparsers_action(node)
        assert action is not None, f"'{name}' has no subcommand layer to descend into"
        node = action.choices[name]
    return node


def _flags(parser):
    flags = set()
    for action in parser._actions:
        flags.update(action.option_strings)
    return flags


def _help_text(parser):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        parser.print_help()
    return buf.getvalue()


def test_data_get_is_a_real_subparser_level():
    """
    The load-bearing assertion. A `source` positional with `choices` would satisfy any
    test about *parsing*, so this checks the thing tab completion actually keys off.
    """
    action = _subparsers_action(_node("get"))
    assert action is not None
    assert set(action.choices) == set(GET_SOURCES) | {"test-data"}


def test_group_nodes_hide_help_and_version_from_completion():
    """
    What the shape buys us: the suppressor marks -h/-v SUPPRESS on group parsers, and
    argcomplete skips SUPPRESSed options, so TAB offers only the database names.
    """
    parser = data_cli.build_parser()
    _suppress_help_version_on_group_parsers(parser)

    get_parser = _subparsers_action(parser).choices["get"]
    suppressed = {
        option
        for action in get_parser._actions
        for option in action.option_strings
        if action.help is argparse.SUPPRESS
    }
    assert HELP_AND_VERSION <= suppressed


def test_leaf_nodes_keep_their_flags_visible_to_completion():
    """The flip side -- suppression must stop at the group level, not cascade."""
    parser = data_cli.build_parser()
    _suppress_help_version_on_group_parsers(parser)

    leaf = _subparsers_action(_subparsers_action(parser).choices["get"]).choices["gtdb-data"]
    assert all(action.help is not argparse.SUPPRESS for action in leaf._actions)


@pytest.mark.parametrize("source_name", GET_SOURCES)
def test_every_source_takes_the_same_flags(source_name):
    """
    `-q/--quiet` is here but not in GToTree's mirror of this menu -- bit's getters
    take a quiet argument and GToTree's don't.
    """
    expected = HELP_AND_VERSION | {"-q", "--quiet", "-f", "--force-update"}
    assert _flags(_node("get", source_name)) == expected


def test_test_data_is_a_real_subparser_level_too():
    """
    `test-data` is a group, not a leaf: its datatype used to be a positional with
    `choices`, which is the shape that pads TAB output with -h/-v.
    """
    action = _subparsers_action(_node("get", "test-data"))
    assert action is not None
    assert set(action.choices) == set(data_cli.TEST_DATA_TYPES)


@pytest.mark.parametrize("datatype_name", sorted(data_cli.TEST_DATA_TYPES))
def test_every_test_datatype_is_a_bare_leaf(datatype_name):
    """
    Nothing to configure once a datatype is named, so these carry only -h/-v. The
    chosen name rides on `args.datatype` -- the subparsers action's dest -- which is
    what dl_test_data() reads, so the worker never learned this changed shape.
    """
    assert _flags(_node("get", "test-data", datatype_name)) == HELP_AND_VERSION


@pytest.mark.parametrize("source_name", GET_SOURCES + ["test-data"])
def test_every_source_renders_one_optional_parameters_section(source_name):
    """
    -h and -v have to land in the same explicit argument group. Adding one to the
    group and the other to the parser also "works", but argparse renders the parser's
    default group first and the help menu comes out with two identically-titled
    sections -- which is exactly what these four leaves used to do.
    """
    assert _help_text(_node("get", source_name)).count("Optional Parameters:") == 1
