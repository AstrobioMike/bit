from setuptools import setup
import glob

# scripts ported to pyproject.toml [project.scripts] entry points
# exclude them from the glob so they don't conflict
ported = {
    "bit/scripts/bit-dl-ncbi-assemblies",
}

setup(
    scripts=[s for s in glob.glob("bit/scripts/*") if s not in ported]
)
