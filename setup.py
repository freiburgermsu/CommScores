# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.rst") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

from os import system
system("pip install git+https://github.com/freiburgermsu/ModelSEEDpy.git")
print("ModelSEEDpy is installed")

setup(
    name="CommScores",
    version="0.0.1",
    description="A Python package for quantifying microbial interactions",
    long_description_content_type="text/x-rst",
    long_description=readme,
    author="Andrew Freiburger",
    author_email="andrewfreiburger@gmail.com",
    url="https://github.com/freiburgermsu/CommScores",
    license=license,
    packages=find_packages(exclude=("docs")),
    # package_data={
    #     "commscores": ["config.cfg", "community/*.html", "core/*.html", "data/*"],
    # },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
    ],
    include_package_data =True,
    install_requires=[
        "modelseedpy", "optlang", "numpy", "deepdiff", "pprint", "sigfig"
    ],
    project_urls={
        "Documentation": "https://modelseedpy.readthedocs.io/en/latest/",
        "Issues": "https://github.com/freiburgermsu/CommScores/issues",
    },
)
