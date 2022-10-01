#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='AtomContr',
    description=(
        "A Python package to understand and explain"
        "atom contribution for small molecules."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    version='0.0.1',
    url='https://github.com/smu-tao-group/AtomContr',
    author='Hao Tian',
    author_email='htian1997@gmail.com',
    license='Apache License 2.0',
    packages=setuptools.find_packages(),
    package_dir={'AtomContr': 'AtomContr'},
    install_requires=[
        "rdkit"
    ],
    zip_safe=False
)
