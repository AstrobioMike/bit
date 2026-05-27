from goatools import obo_parser
import os
import pandas as pd
import sys
from bit.modules.go.get_go_dbs import get_go_data


def check_go_dbs():

    get_go_data(quiet=True)
    go_data_dir = os.environ["GO_DB_DIR"]

    return go_data_dir


def resolve_obo_path(obo_arg):

    go_data_dir = check_go_dbs()

    if obo_arg == "goslim_metagenomics":
        return os.path.join(go_data_dir, "goslim_metagenomics.obo")
    elif obo_arg == "go_basic":
        return os.path.join(go_data_dir, "go-basic.obo")
    else:
        return obo_arg


def load_obo(obo_path, load_obsolete=False):

    print("\n\tGO obo file being used:")
    go = obo_parser.GODag(obo_path, load_obsolete=load_obsolete)
    print("")

    return go


def build_go_df(go):

    def get_general_info(go_id):
        go_term = go[go_id]
        return [go_id, go_term.namespace, go_term.depth, go_term.name, go_term.is_obsolete]

    table_list = [get_general_info(go_id) for go_id in go.keys()]
    header = ["GO_term", "namespace", "depth", "name", "is_obsolete"]
    GO_df = pd.DataFrame(table_list, columns=header).sort_values(by=["depth"])

    obsolete_terms = set()

    if True in GO_df["is_obsolete"].values:
        obsolete_terms = set(GO_df.loc[GO_df["is_obsolete"] == True, "GO_term"])
        GO_df = GO_df.loc[GO_df["is_obsolete"] == False, ["GO_term", "namespace", "depth", "name"]]

    return GO_df, obsolete_terms


def count_terms(annots_file, term_counts_dict, obsolete_terms, obo_arg):

    with open(annots_file, "r") as annots:
        for line in annots:

            if len(line.strip().split("\t")) <= 1 or line.startswith("#"):
                continue

            terms = line.strip().split("\t")[1].replace(" ", "")
            for term in terms.split(";"):
                if term in obsolete_terms:
                    continue
                if term in term_counts_dict:
                    term_counts_dict[term] += 1
                else:
                    print(f'\n    So the GO term "{term}" shows up in the provided annotation file')
                    print(f'    "{annots_file}", but it is not present in the obo reference')
                    print(f'    file used: "{obo_arg}". We aren\'t sure how to deal with this. You might')
                    print('    want to pass a different obo file to the `-g` argument.\n')
                    print('        Exiting for now :(\n')
                    sys.exit(1)


def add_percentages(df):

    df = df.copy()
    df["term_perc_of_annotated"] = df["term_counts"] / df["term_counts"].sum() * 100
    return df


_TERM_HEADER = ["GO id", "namespace", "depth", "name"]


def _term_row(go, go_id):
    t = go[go_id]
    return [go_id, t.namespace, t.depth, t.name]


def _print_related_terms(go, term_ids, label, go_id):
    if term_ids:
        rows = [_term_row(go, t) for t in term_ids]
        df = pd.DataFrame(rows)
        print(f"\n{label} info:")
        print(df.to_string(index=False, header=_TERM_HEADER))
    else:
        print(f"\nThere are no {label.lower()} for {go_id}.")


def write_go_tables(GO_df, args):

    out_prefix = args.output_prefix

    def _write_if_nonempty(df, path, label):
        df_filtered = df if args.keep_zeros else df[df["term_counts"] > 0]
        if len(df_filtered.index) == 0:
            if label:
                print(f"\n\tThere were no counts to any {label} terms, so that table wasn't reported.\n")
        else:
            with open(path, "w") as out:
                out.write(df_filtered.to_csv(index=False, sep="\t"))
        return df_filtered

    combined = _write_if_nonempty(GO_df, out_prefix + ".tsv", None)
    if len(combined.index) == 0:
        print("\n\tThere were no counts to any terms :( Sure we're working with the right files here?\n")
        sys.exit(0)

    if args.by_namespace:
        for ns_key, label, suffix in (
            ("biological_process", "Biological Process", "-BP"),
            ("molecular_function", "Molecular Function", "-MF"),
            ("cellular_component", "Cellular Component", "-CC"),
        ):
            ns_df = GO_df[GO_df["namespace"] == ns_key]
            ns_df = add_percentages(ns_df)
            _write_if_nonempty(ns_df, out_prefix + suffix + ".tsv", label)


def get_term_info(args):

    pd.set_option('display.max_colwidth', None)

    go_obo = resolve_obo_path(args.GO_obo_file)
    go = load_obo(go_obo)

    input_go_id = args.GO_term
    if not input_go_id.startswith("GO:"):
        input_go_id = "GO:" + input_go_id

    try:
        input_go_term = go[input_go_id]
    except:
        print(f"{input_go_id} does not seem to be in the GO database :(\n")
        sys.exit()

    input_df = pd.DataFrame([_term_row(go, input_go_id)], columns=_TERM_HEADER)
    print("Input GO term info:")
    print(input_df.to_string(index=False))

    _print_related_terms(go, input_go_term.get_all_parents(), "Parent terms", input_go_id)

    if not args.parents_only:
        _print_related_terms(go, input_go_term.get_all_children(), "Child terms", input_go_id)

    print("")


def summarize_annotations(args):

    go_obo = resolve_obo_path(args.GO_obo_file)
    go = load_obo(go_obo, load_obsolete=True)

    GO_df, obsolete_terms = build_go_df(go)

    term_counts_dict = {key: 0 for key in GO_df["GO_term"]}
    count_terms(args.input_GO_annotations, term_counts_dict, obsolete_terms, args.GO_obo_file)

    GO_df["term_counts"] = term_counts_dict.values()
    GO_df = add_percentages(GO_df)

    write_go_tables(GO_df, args)


def combine_summaries(args):

    # file name is key, and sample name is value
    all_samples = {}

    # setting sample names
    if len(args.sample_names) == 0:
        for file in args.input_files:
            curr_sample = os.path.basename(file).rsplit('.', 1)[0]

            if file in all_samples:
                print('\n    It seems the file "' + file + '" is trying to get in here twice.')
                print("\n    That's not gonna fly :(\n")
                sys.exit(1)

            all_samples[file] = curr_sample

    else:

        # checking if sample names provided the length equals the number of input files
        if len(args.sample_names) != len(args.input_files):
            print("\n    It seems the number of provided sample names doesn't match the number of provided input files :(")
            print("\n    Check usage with `bit go combine-summaries -h`.\n")
            sys.exit(0)

        for i, curr_sample in enumerate(args.sample_names):
            all_samples[args.input_files[i]] = curr_sample

    # keeping a nested dictionary of info for all GO terms that show up in any table
    GO_dict = {}

    # building counts/percents table
    building_tab = pd.DataFrame(columns=["GO_term"])

    ## working on each file
    for sample_key in all_samples:

        # reading current file into pandas dataframe
        curr_tab = pd.read_csv(sample_key, sep="\t")

        # adding to building GO dictionary of all GO terms in the input tables
        for row in curr_tab.itertuples():
            if row.GO_term not in GO_dict:
                GO_dict[row.GO_term] = {'namespace': row.namespace, 'depth': row.depth, 'name': row.name}

        # trimming down to current sample columns
        sample_name = all_samples[sample_key]
        curr_sub_tab = curr_tab[["GO_term", "term_counts", "term_perc_of_annotated"]].rename(columns={
            "term_counts": sample_name + "_counts",
            "term_perc_of_annotated": sample_name + "_perc_of_annotated",
        })

        # merging with master tab on GO_term
        building_tab = building_tab.merge(curr_sub_tab, on="GO_term", how="outer")

    ## replacing NAs with 0s
    building_tab = building_tab.fillna(0)

    ## making GO info dict into dataframe and merging into final table
    go_df = pd.DataFrame.from_dict(GO_dict, orient="index")
    # moving index to column and renaming
    go_df.reset_index(inplace=True)
    go_df.rename(columns = {'index': 'GO_term'}, inplace=True)
    # merging
    final_tab = go_df.merge(building_tab, on="GO_term", how="outer")
    # sorting
    final_tab.sort_values(by=["namespace", "depth"], inplace=True)

    ## writing out
    with open(args.output_file, "w") as out:
        out.write(final_tab.to_csv(index=False, sep="\t"))


def slim_terms(args):

    import subprocess

    initial_obo = resolve_obo_path(args.initial_GO_obo_file)
    slim_obo = resolve_obo_path(args.slimmed_GO_obo_file)


    cmd = ["map_to_slim.py", "--association_file", args.input_GO_annotations, "--slim_out", args.mode, initial_obo, slim_obo]
    with open(args.output_file, "w") as output:
        subprocess.run(cmd, stdout=output)

