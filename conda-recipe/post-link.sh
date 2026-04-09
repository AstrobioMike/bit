printf "will cite" | parallel --citation 2&> /dev/null

# removing coreutils ls so the user-system one is used still
rm -f ${CONDA_PREFIX}/bin/ls

# setting up tab-completion for the bit commands with subcommands
for cmd in bit-extract-seqs bit-ez-screen; do
    eval "$(register-python-argcomplete "$cmd")"
done
