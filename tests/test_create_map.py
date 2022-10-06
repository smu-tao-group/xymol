#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test create_map function
"""

import os
from xymol import XYMOL
from deepchem import deepchem as dc


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
    xymol.create_map(dc.feat.ConvMolFeaturizer, model, file_name)

    # check file exist
    assert os.path.exists(file_name)
    # remove file
    os.system(f"rm {file_name}")

