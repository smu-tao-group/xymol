#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test create_map function
"""

import os
from rdkit import Chem
import numpy as np
from xymol import XYMOL

SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"
TEST_WEIGHT = 4.0

class TestModel:
    def predict(self, data) -> float:

        # return 0.0 if predicting for the original smiles
        if data[0][0] == 1:
            return [0.0]

        # return the test weight if predicting for drop_atoms
        return [TEST_WEIGHT] * len(data)

class TestFeaturizer:
    def featurize(self, smiles) -> list:

        #return list of zeros if drop_atoms is passed
        if isinstance(smiles, np.ndarray):
            feat = []
            for _ in range(len(smiles)):
                feat.append([0, 0, 0])
            return feat

        #return list of ones if the original smiles is passed
        return [[1, 1, 1]]


def test_create_map():
    xymol = XYMOL(SMILES)
    file_name = "molecule.png"

    weights = xymol.create_map(TestFeaturizer(), TestModel(), file_name)

    mol = Chem.MolFromSmiles(SMILES)
    assert len(weights) == len(mol.GetAtoms())

    # check that all values in weights equal correct test difference (drop_atom pred minus parent_pred)
    assert weights[0] == TEST_WEIGHT
    assert weights.count(weights[0]) == len(weights)

    # check file exists
    assert os.path.exists(file_name)
    # remove file
    os.remove(file_name)
