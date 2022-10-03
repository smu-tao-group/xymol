#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


# get version info
def _get_version():
    with open('xymol/__init__.py') as fp:
        for line in fp:
            if line.startswith('__version__'):
                g = {}
                exec(line, g)
                return g['__version__']


with open("README.md", "r") as fh:
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
    ],
    python_requires='>=3.7'
)
