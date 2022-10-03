#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test drop_each_bond function.
"""

from rdkit import Chem
from xymol import XYMOL


SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"
SMILES_AROMATIC = "Fc1ccccc1F"


def test_drop_each_bond():
    xymol = XYMOL(SMILES)
    _, bond_smiles = xymol.drop_each_bond()

    # create a rdkit mol
    mol = Chem.MolFromSmiles(SMILES)

    assert len(bond_smiles) == len(mol.GetBonds())


def test_drop_each_bond_sanitize():
    xymol = XYMOL(SMILES_AROMATIC)
    _, bond_smiles = xymol.drop_each_bond()
    assert len(bond_smiles) == len(xymol.mol.GetBonds())
