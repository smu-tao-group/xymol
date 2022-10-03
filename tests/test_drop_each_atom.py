#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test drop_each_atom function.
"""

from rdkit import Chem
from xymol import XYMOL


SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"
SMILES_AROMATIC = "Fc1ccccc1F"


def test_drop_each_atom():
    xymol = XYMOL(SMILES)
    atom_smiles = xymol.drop_each_atom()

    # create a rdkit mol
    mol = Chem.MolFromSmiles(SMILES)

    assert len(atom_smiles) == len(mol.GetAtoms())


def test_drop_each_atom_sanitize():
    xymol = XYMOL(SMILES_AROMATIC)
    atom_smiles = xymol.drop_each_atom()
    assert len(atom_smiles) == len(xymol.mol.GetAtoms())
