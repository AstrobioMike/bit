# dependencies are handled by conda installation

from setuptools import setup, find_packages
import glob

scripts = glob.glob("bit/scripts/*")

setup(
    name="bit",
    version="1.12.4",
    description="A set of bioinformatics scripts and workflows",
    license="GPLv3",
    author="Mike Lee",
    author_email="Michael.Lee0517@gmail.com",
    url="https://github.com/AstrobioMike/bit",
    packages=find_packages(),
    scripts=scripts,
    include_package_data=True,
    package_data={"bit": ["tests/data/*"]},
)
