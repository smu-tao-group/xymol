#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rdkit import Chem
from xymol import XYMOL


smiles = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"


def test_drop_each_atom():
    xymol = XYMOL(smiles)
    atom_smiles = xymol.drop_each_atom()

    # create a rdkit mol
    mol = Chem.MolFromSmiles(smiles)

    assert len(atom_smiles) == len(mol.GetAtoms())
