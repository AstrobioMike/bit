from setuptools import setup
import glob

setup(
    scripts=glob.glob("bit/scripts/*")
)
