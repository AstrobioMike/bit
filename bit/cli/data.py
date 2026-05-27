import sys
import argparse
import argcomplete # type: ignore
from bit.cli.common import CustomRichHelpFormatter, add_help, add_version_arg


def build_parser(parent_subparsers=None):

    desc = """
        This program manages bit-utilized databases and their location settings. See subcommand-specific
        help menus for more info.
        """

    if parent_subparsers is not None:
        parser = parent_subparsers.add_parser(
            "data",
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False,
        )
    else:
        parser = argparse.ArgumentParser(
            description=desc,
            formatter_class=CustomRichHelpFormatter,
            add_help=False
        )

    add_help(parser)

    add_version_arg(parser)

    subparsers = parser.add_subparsers(dest="subcommand", required=True, metavar='')
    parser.subparsers = subparsers


    ############################################################################
    ### get subcommand ###
    ############################################################################

    get_desc = """
        This subcommand downloads or updates bit's required databases. All sub-subcommands (except test-data) accept
        optional `-q/--quiet` and `-f/--force-update` flags.
        """

    get_parser = subparsers.add_parser(
        "get",
        help="Download/update bit-utilized databases, or get test data",
        description=get_desc,
        epilog="Ex. usage: `bit-data get ncbi-assembly-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(get_parser)

    add_version_arg(get_parser)

    get_subparsers = get_parser.add_subparsers(dest="get_action", required=True, metavar='')
    get_parser.subparsers = get_subparsers


    def add_get_common_args(group):
        group.add_argument(
            "-q", "--quiet",
            help="Exit silently if data already present",
            action="store_true"
        )
        group.add_argument(
            "-f", "--force-update",
            help="Re-download/update the data even if already present",
            action="store_true"
        )


    ### get go-dbs ###

    get_go_dbs_desc = """
        This subcommand downloads or updates the Gene Ontology (GO) databases.
        """

    get_go_dbs_parser = get_subparsers.add_parser(
        "go-dbs",
        help="Download or update GO databases",
        description=get_go_dbs_desc,
        epilog="Ex. usage: `bit-data get go-dbs`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_go_dbs_optional = get_go_dbs_parser.add_argument_group("Optional Parameters")
    add_get_common_args(get_go_dbs_optional)
    add_help(get_go_dbs_optional)

    add_version_arg(get_go_dbs_parser)

    get_go_dbs_parser.set_defaults(func="get_go_dbs")


    ### get gtdb-data ###

    get_gtdb_data_desc = """
        This subcommand downloads or updates the GTDB metadata.
        """

    get_gtdb_data_parser = get_subparsers.add_parser(
        "gtdb-data",
        help="Download or update GTDB metadata",
        description=get_gtdb_data_desc,
        epilog="Ex. usage: `bit-data get gtdb-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_gtdb_data_optional = get_gtdb_data_parser.add_argument_group("Optional Parameters")
    add_get_common_args(get_gtdb_data_optional)
    add_help(get_gtdb_data_optional)

    add_version_arg(get_gtdb_data_parser)

    get_gtdb_data_parser.set_defaults(func="get_gtdb_data")


    ### get ncbi-assembly-data ###

    get_ncbi_assembly_desc = """
        This subcommand downloads or updates NCBI's assembly-summary tables.
        """

    get_ncbi_assembly_parser = get_subparsers.add_parser(
        "ncbi-assembly-data",
        help="Download or update NCBI assembly-summary tables",
        description=get_ncbi_assembly_desc,
        epilog="Ex. usage: `bit-data get ncbi-assembly-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_ncbi_assembly_optional = get_ncbi_assembly_parser.add_argument_group("Optional Parameters")
    add_get_common_args(get_ncbi_assembly_optional)
    add_help(get_ncbi_assembly_optional)

    add_version_arg(get_ncbi_assembly_parser)

    get_ncbi_assembly_parser.set_defaults(func="get_ncbi_assembly_data")


    ### get ncbi-tax-data ###

    get_ncbi_tax_desc = """
        This subcommand downloads or updates NCBI taxonomy data (used by TaxonKit).
        """

    get_ncbi_tax_parser = get_subparsers.add_parser(
        "ncbi-tax-data",
        help="Download or update NCBI taxonomy data",
        description=get_ncbi_tax_desc,
        epilog="Ex. usage: `bit-data get ncbi-tax-data`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_ncbi_tax_optional = get_ncbi_tax_parser.add_argument_group("Optional Parameters")
    add_get_common_args(get_ncbi_tax_optional)
    add_help(get_ncbi_tax_optional)

    add_version_arg(get_ncbi_tax_parser)

    get_ncbi_tax_parser.set_defaults(func="get_ncbi_tax_data")


    ### get test-data ###

    get_test_data_desc = """
        This subcommand downloads test-data files.
        """

    get_test_data_parser = get_subparsers.add_parser(
        "test-data",
        help="Download test data",
        description=get_test_data_desc,
        epilog="Ex. usage: `bit-data get test-data genome`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    get_test_data_required = get_test_data_parser.add_argument_group("Required Parameters")
    get_test_data_optional = get_test_data_parser.add_argument_group("Optional Parameters")

    get_test_data_required.add_argument(
        "datatype",
        choices=["genome", "metagenome"],
        help="What type of test data you'd like to download",
    )

    add_help(get_test_data_optional)

    add_version_arg(get_test_data_optional)

    get_test_data_parser.set_defaults(func="get_test_data")


    ############################################################################
    ### locations subcommand ###
    ############################################################################

    locations_desc = """
        This subcommand checks or sets the location environment variables required by bit.
        """

    locations_parser = subparsers.add_parser(
        "locations",
        help="Check or set data-location environment variables",
        description=locations_desc,
        epilog="Ex. usage: `bit-data locations check`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    add_help(locations_parser)

    add_version_arg(locations_parser)

    locations_subparsers = locations_parser.add_subparsers(dest="locations_action", required=True, metavar='')
    locations_parser.subparsers = locations_subparsers


    ### locations check ###

    locations_check_desc = """
        This subcommand reports the current data-location environment variables and whether the
        paths are writable.
        """

    locations_check_parser = locations_subparsers.add_parser(
        "check",
        help="Report current data-location environment variables",
        description=locations_check_desc,
        epilog="Ex. usage: `bit-data locations check`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    locations_check_optional = locations_check_parser.add_argument_group("Optional Parameters")
    add_help(locations_check_optional)

    add_version_arg(locations_check_optional)
    locations_check_parser.set_defaults(func="locations_check")


    ### locations set ###

    locations_set_desc = """
        This subcommand interactively sets the data-location environment variables.
        """

    locations_set_parser = locations_subparsers.add_parser(
        "set",
        help="Interactively set data-location environment variables",
        description=locations_set_desc,
        epilog="Ex. usage: `bit-data locations set`",
        formatter_class=CustomRichHelpFormatter,
        add_help=False
    )

    locations_set_optional = locations_set_parser.add_argument_group("Optional Parameters")
    add_help(locations_set_optional)

    add_version_arg(locations_set_optional)

    locations_set_parser.set_defaults(func="locations_set")

    return parser


################################################################################


def main():

    parser = build_parser()
    argcomplete.autocomplete(parser)

    if len(sys.argv) == 1: # pragma: no cover
        parser.print_help(sys.stderr)
        sys.exit(0)

    # handling no args for a top-level subcommand so the appropriate help menu is printed
    if len(sys.argv) == 2:
        cmd = sys.argv[1]

        if cmd in ("-h", "--help"):
            parser.print_help(sys.stderr)
            sys.exit(0)

        if cmd in ("-v", "--version"):
            from bit.modules.general import report_version
            report_version()
            sys.exit(0)

        if cmd in parser.subparsers.choices:
            if sys.stdin.isatty():
                parser.subparsers.choices[cmd].print_help(sys.stderr)
                sys.exit(0)
            # else: stdin is being piped, fall through to parse_args()
        else:
            print(f"\n  Invalid subcommand provided: '{cmd}'\n\n  See help below.\n", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)

    # handling no args for a second-level subcommand so the appropriate help menu is printed
    if len(sys.argv) == 3:
        cmd = sys.argv[1]
        sub_cmd = sys.argv[2]

        if cmd in parser.subparsers.choices:
            top_sub_parser = parser.subparsers.choices[cmd]

            if hasattr(top_sub_parser, 'subparsers'):

                if sub_cmd in ("-h", "--help"):
                    top_sub_parser.print_help(sys.stderr)
                    sys.exit(0)

                if sub_cmd in top_sub_parser.subparsers.choices:
                    # only print help if this sub-subcommand itself has required args
                    sub_sub_parser = top_sub_parser.subparsers.choices[sub_cmd]
                    has_required = any(
                        a.required
                        for a in sub_sub_parser._actions
                        if hasattr(a, 'required')
                    )
                    if has_required and sys.stdin.isatty():
                        sub_sub_parser.print_help(sys.stderr)
                        sys.exit(0)
                    # else: command is complete as-is, fall through to parse_args()
                else:
                    print(f"\n  Invalid subcommand provided: '{sub_cmd}'\n\n  See help below.\n", file=sys.stderr)
                    top_sub_parser.print_help(sys.stderr)
                    sys.exit(1)

    args = parser.parse_args()

    from bit.modules.data_locations import (check_and_report_env_variables, set_env_variables,
                                            modify_conda_activate_startup_script, notify_to_reactivate_conda)
    from bit.modules.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data
    from bit.modules.ncbi.get_ncbi_tax_data import get_ncbi_tax_data
    from bit.modules.go.get_go_dbs import get_go_data
    from bit.modules.gtdb.get_gtdb_data import get_gtdb_data
    from bit.modules.get_test_data import dl_test_data

    func_map = {
        "locations_check"      : lambda a: check_and_report_env_variables(),
        "locations_set"        : lambda a: (
                                    modify_conda_activate_startup_script(set_env_variables()),
                                    notify_to_reactivate_conda()
                                 ),
        "get_ncbi_assembly_data" : lambda a: get_ncbi_assembly_data(force_update=a.force_update, quiet=a.quiet),
        "get_ncbi_tax_data"      : lambda a: get_ncbi_tax_data(force_update=a.force_update, quiet=a.quiet),
        "get_go_dbs"             : lambda a: get_go_data(force_update=a.force_update, quiet=a.quiet),
        "get_gtdb_data"          : lambda a: get_gtdb_data(force_update=a.force_update, quiet=a.quiet),
        "get_test_data"          : lambda a: dl_test_data(a),
    }

    func_map[args.func](args)
