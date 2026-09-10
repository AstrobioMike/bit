import argparse
import importlib
import pytest # type: ignore
from bit.cli.bit import SUBCOMMAND_MAP, build_parser


@pytest.mark.parametrize("module_path", sorted(set(SUBCOMMAND_MAP.values())))
def test_subcommand_module_imports(module_path):
    importlib.import_module(module_path)


def _subcommand_layers_smuggled_in_as_positionals(parser, path):
    """
    Walk the tree and collect any parser that offers its subcommands as a positional
    with `choices` instead of as a real subparser level.
    """
    subparser_actions = [
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    ]

    if subparser_actions:
        offenders = []
        for name, child in subparser_actions[0].choices.items():
            offenders += _subcommand_layers_smuggled_in_as_positionals(
                child, f"{path} {name}")
        return offenders

    return [
        f"{path} (<{action.dest}>)"
        for action in parser._actions
        if not action.option_strings and getattr(action, "choices", None)
    ]


def test_no_parser_uses_a_positional_with_choices_as_a_subcommand_layer():
    """
    Tab completion is driven by argcomplete walking this tree, and
    _suppress_help_version_on_group_parsers() recognizes a "group" only via
    argparse._SubParsersAction. A parser that instead offers its subcommands as a
    positional with `choices` parses identically and reads fine under `-h`, but TAB
    after it pads the choices with -h/-v and whatever else that parser defines.

    `bit get-workflow` and `bit data get test-data` were both this shape -- TAB after
    either padded the choices with -h/-v (and -l/-w in get-workflow's case).

    Tree-wide on purpose rather than scoped to the parsers that had the problem: the
    shape is an easy one to reach for, and the symptom only ever shows up at a
    terminal, so it goes unnoticed until someone happens to hit TAB there.
    """
    offenders = _subcommand_layers_smuggled_in_as_positionals(build_parser(), "bit")
    assert not offenders, (
        "these parsers use a positional-with-choices as a subcommand layer, so "
        f"argcomplete will pad their choices with flags: {sorted(offenders)}"
    )
