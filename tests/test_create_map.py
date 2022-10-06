#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test create_map function
"""

import os
from rdkit import Chem
from deepchem import deepchem as dc
from xymol import XYMOL


SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"

metric = dc.metrics.mean_absolute_error
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='GraphConv', reload=False)
train_dataset, valid_dataset, test_dataset = datasets

def test_create_map():
    n_tasks = len(tasks)
    model = dc.models.GraphConvModel(n_tasks=n_tasks, mode="regression")
    model.fit(train_dataset, nb_epoch=5)

    xymol = XYMOL(SMILES)
    file_name = "molecule.png"

    weights = xymol.create_map("GraphConv", model, file_name)

    mol = Chem.MolFromSmiles(SMILES)
    assert len(weights) == len(mol.GetAtoms())

    # check file exists
    assert os.path.exists(file_name)
    # remove file
    os.remove(file_name)
