#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test create_map function
"""

import os
import pytest
import numpy as np
from xymol import XYMOL

SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"

class TestModel:
    def predict(self, data) -> float:

        # return 0.0 if predicting for the original smiles
        if data[0][0] == 1:
            return [0.0]

        # return 1.0 if predicting for drop_atoms
        return [1.0] * len(data)

class TestFeaturizer:
    def featurize(self, smiles) -> list:

        # return list of zeros if drop_atoms is passed
        if isinstance(smiles, np.ndarray):
            feat = []
            for _ in range(len(smiles)):
                feat.append([0, 0, 0])
            return feat

        # return list of ones if the original smiles is passed
        return [[1, 1, 1]]

xymol = XYMOL(SMILES)

def test_create_map():
    file_name = "molecule.png"

    xymol.create_map(TestFeaturizer(), TestModel(), file_name)

    # check file exists
    assert os.path.exists(file_name)
    # remove file
    os.remove(file_name)

def test_create_map_without_featurizer():
    with pytest.raises(ValueError):
        xymol.create_map("", TestModel)

def test_create_map_without_model():
    with pytest.raises(ValueError):
        xymol.create_map(TestFeaturizer, "")
