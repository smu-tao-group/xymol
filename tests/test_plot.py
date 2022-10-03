#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from xymol import XYMOL


SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"


def test_plot_mol():
    xymol = XYMOL(SMILES)
    file_name = "molecule.png"
    xymol.plot_mol(file_name)

    # check file exist
    assert os.path.exists(file_name)
    # remove file
    os.system(f"rm {file_name}")


def test_plot_similarity_map():
    xymol = XYMOL(SMILES)
    weights = [0.1] * len(xymol.mol.GetAtoms())
    file_name = "molecule.png"
    xymol.plot_similarity_map(weights=weights, file_name=file_name)

    # check file exist
    assert os.path.exists(file_name)
    # remove file
    os.system(f"rm {file_name}")
