from pathlib import Path
from itertools import zip_longest


def print_two_columns(items, indent="            "):
    names = sorted(str(p if isinstance(p, str) else p.name) for p in items)

    col_width = max(len(name) for name in names) + 4

    half = (len(names) + 1) // 2
    left = names[:half]
    right = names[half:]

    for a, b in zip(left, right + [""] * (len(left) - len(right))):
        print(f"{indent}{a:<{col_width}}{b}")


def pytest_sessionstart(session):

    scripts_dir = Path(__file__).resolve().parents[1] / "scripts"

    script_files = [p for p in scripts_dir.iterdir() if p.is_file()]

    bash_scripts = []
    helper_scripts = []
    python_scripts = []
    modularized = []
    not_modularized = []

    for file in script_files:
        # pseudo-ignoring bit-version for this purpose
        if file.name == "bit-version":
            python_scripts.append(file)
            modularized.append(file)
            continue

        # checking if helper script
        if file.name.startswith("helper-bit-"):
            helper_scripts.append(file)
            continue  # moving on if helper script

        # getting first line to check if bash script
        with file.open() as fh:
            first_line = fh.readline().strip()

        if first_line.startswith("#!/usr/bin/env bash"):
            bash_scripts.append(file)
            continue  # moving on if bash script

        # the rest are all python if not bash
        python_scripts.append(file)

        # now reading more to check if i've modularized it yet (based on if __name__ == "__main__" is in first 10 lines)
        text = file.read_text().splitlines()
        first_ten = "\n".join(text[:10])

        if 'if __name__ == "__main__"' in first_ten:
            modularized.append(file)
        else:
            not_modularized.append(file)

    print("\n  ====================================================================================================")
    print("  ================================= Progress on modularizing scripts =================================")
    print("  ====================================================================================================\n")

    print(f"    Total script files: {len(script_files)}\n")
    print(f"    Total helper scripts: {len(helper_scripts)}")
    print(f"    Total bash scripts: {len(bash_scripts)}")
    print(f"    Total python scripts: {len(python_scripts)}")
    print(f"    Python scripts modularized: {len(modularized)}")
    print(f"    Python scripts not yet modularized: {len(not_modularized)}\n")
    print(f"    Percent of the way done: {round(len(modularized) / len(python_scripts) * 100, 2)}%\n")

    if not_modularized:
        print("        Files not yet modularized:\n")
        print_two_columns(not_modularized)

    print("\n  ====================================================================================================\n")
