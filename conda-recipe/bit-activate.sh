#!/usr/bin/env bash

### handling setting up of argcomplete for bit commands with subcommands ###

ARGCOMPLETE_COMMANDS=(
    bit-extract-seqs
    bit-ez-screen
    bit-genbank
)

# checking interactive shell
case "$-" in
    *i*) ;;
    *) return 0 2>/dev/null || exit 0 ;;
esac

# checking for argcomplete helper
if ! command -v register-python-argcomplete >/dev/null 2>&1; then
    return 0 2>/dev/null || exit 0
fi

for cmd in "${ARGCOMPLETE_COMMANDS[@]}"; do
    if command -v "$cmd" >/dev/null 2>&1; then
        eval "$(register-python-argcomplete "$cmd")"
    fi
done
