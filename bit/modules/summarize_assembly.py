import os


def summarize_assemblies(input_assemblies, output_tsv=False, transpose_output_tsv=False, jobs=1):

    import pandas as pd # type: ignore

    df, use_paths_instead_of_basenames = setup_master_df(input_assemblies)
    df = summarize(df, input_assemblies, use_paths_instead_of_basenames, jobs=jobs)

    if output_tsv:
        if transpose_output_tsv:
            df = df.T
            df.to_csv(output_tsv, sep="\t", index=False, na_rep = "NA")
        else:
            df.to_csv(output_tsv, sep="\t", header=False, na_rep = "NA")

    else:
        display_df = df.copy()
        for idx in display_df.index:
            if idx in ("Assembly", "GC content"):
                continue
            for col in display_df.columns:
                val = display_df.at[idx, col]
                if pd.notna(val):
                    display_df.at[idx, col] = f"{int(val):,}"

        print("")
        print(display_df.to_string(header=False))
        print("")


def setup_master_df(input_assemblies):

    import pandas as pd # type: ignore

    df_colnames = []

    for assembly in input_assemblies:
        assembly_base = os.path.basename(assembly)
        df_colnames.append(assembly_base.rsplit(".", 1)[0])

    # checking for a situation where inputs may have the same basename, due to being from different directories
    # if so, setting a flag and reporting them as the full input paths instead of just basenames
    use_paths_instead_of_basenames = False
    for assembly in df_colnames:
        num_occurences = 0
        for assembly_2 in df_colnames:
            if assembly_2 == assembly:
                num_occurences += 1
        if num_occurences > 1:
            use_paths_instead_of_basenames = True
    if use_paths_instead_of_basenames:
        df_colnames = []
        for assembly in input_assemblies:
            df_colnames.append(assembly.rsplit(".", 1)[0])
            # df_colnames.append(os.path.splitext(assembly)[0])

    df_index = ["Assembly", "Total contigs", "Total length", "Ambiguous characters",
                "GC content", "Maximum contig length", "Minimum contig length", "N50",
                "N75", "N90", "L50", "L75", "L90", "Num. contigs >= 100",
                "Num. contigs >= 500", "Num. contigs >= 1000", "Num. contigs >= 5000",
                "Num. contigs >= 10000", "Num. contigs >= 50000", "Num. contigs >= 100000"]

    return pd.DataFrame(columns=df_colnames, index=df_index), use_paths_instead_of_basenames


def summarize(df, input_assemblies, use_paths_instead_of_basenames, jobs=1):

    all_stats = compute_all_assembly_stats(input_assemblies, jobs=jobs)

    for assembly, stats in zip(input_assemblies, all_stats):
        if use_paths_instead_of_basenames:
            assembly_name = assembly.rsplit(".", 1)[0]
            # assembly_name = os.path.splitext(assembly)[0]
        else:
            assembly_base = os.path.basename(assembly)
            assembly_name = assembly_base.rsplit(".", 1)[0]

        df.at["Assembly", str(assembly_name)] = assembly_name

        # stats is empty for empty input files (which can happen if an assembly produced no contigs)
        # this will leave it in the table, but with NAs (written out as "NA")
        for row, value in stats.items():
            df.at[row, str(assembly_name)] = value

    return df


def compute_all_assembly_stats(input_assemblies, jobs=1):
    """
    Returns a list of per-assembly stats dicts, in the same order as input_assemblies.
    Runs serially if only 1 job is possible, otherwise spreads assemblies across worker processes.
    """

    jobs = max(1, min(jobs, len(input_assemblies)))

    if jobs == 1:
        return [assembly_stats(assembly) for assembly in input_assemblies]

    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    # using 'spawn' rather than the default 'fork' to match gen-reads and stay safe if
    # this is ever called from a multi-threaded context (fork can deadlock there and
    # warns on Python 3.12+); executor.map keeps results in input order
    mp_ctx = multiprocessing.get_context("spawn")
    ex = ProcessPoolExecutor(max_workers=jobs, mp_context=mp_ctx)
    try:
        return list(ex.map(assembly_stats, input_assemblies))
    except KeyboardInterrupt:
        ex.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        ex.shutdown(wait=True)


def assembly_stats(assembly):
    """
    Returns a dict of summary stats for one assembly, keyed by the master df row names.
    Returns an empty dict if the file is empty. Kept top-level so it can be pickled for worker processes.
    """

    import pyfastx # type: ignore
    import tempfile

    if os.stat(assembly).st_size == 0:
        return {}

    stats = {}

    # pyfastx writes a .fxi index file next to the input, this can cause problems if it is from a read-only area
    # using a symlink in a temp dir avoid this
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_link = os.path.join(tmpdir, os.path.basename(assembly))
        os.symlink(os.path.abspath(assembly), tmp_link)
        fasta = pyfastx.Fasta(tmp_link)

        stats["Total contigs"] = len(fasta)
        stats["Total length"] = fasta.size

        num_ambiguous_chars = 0
        for key in fasta.composition:
            if key not in ["A","T","G","C"]:
                num_ambiguous_chars += fasta.composition[key]

        stats["Ambiguous characters"] = num_ambiguous_chars
        stats["GC content"] = round(fasta.gc_content, 2)
        stats["Maximum contig length"] = len(fasta.longest)
        stats["Minimum contig length"] = len(fasta.shortest)

        info_at_50 = fasta.nl(50)
        info_at_75 = fasta.nl(75)
        info_at_90 = fasta.nl(90)
        stats["N50"] = info_at_50[0]
        stats["N75"] = info_at_75[0]
        stats["N90"] = info_at_90[0]
        stats["L50"] = info_at_50[1]
        stats["L75"] = info_at_75[1]
        stats["L90"] = info_at_90[1]

        for min_len in (100, 500, 1000, 5000, 10000, 50000, 100000):
            stats[f"Num. contigs >= {min_len}"] = fasta.count(min_len)

    return stats
