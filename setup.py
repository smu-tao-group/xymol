#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


# get version info
def _get_version():
    with open("xymol/__init__.py", encoding="utf-8") as init_file:
        for line in init_file:
            if line.startswith("__version__"):
                version_info = {}
                exec(line, version_info)
                return version_info["__version__"]
    raise ValueError("version number is missing!")


with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='xymol',
    description=(
        "eXplain Your MOLecule (XYMOL): A Python package to understand"
        " and explain atom/bond contributions of small molecules"
        " in machine learning models."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=_get_version(),
    url='https://github.com/smu-tao-group/xymol',
    author='Hao Tian',
    author_email='htian1997@gmail.com',
    license='Apache License 2.0',
    packages=setuptools.find_packages(),
    install_requires=[
        "rdkit",
        "numpy",
        "matplotlib",
        "deepchem"
    ],
    python_requires='>=3.7'
)
