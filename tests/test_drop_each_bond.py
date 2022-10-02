#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rdkit import Chem
from xymol import XYMOL


smiles = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"


def test_drop_each_bond():
    xymol = XYMOL(smiles)
    _, bond_smiles = xymol.drop_each_bond()

    # create a rdkit mol
    mol = Chem.MolFromSmiles(smiles)

    assert len(bond_smiles) == len(mol.GetBonds())
