#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test show and hide atom number.
"""

from xymol import XYMOL


SMILES = "CCC1(CCC(C)C)C(=O)NC(=O)NC1=O"


def test_show_atom_number():
    xymol = XYMOL(SMILES)
    xymol.show_atom_number()
    for atom in xymol.mol.GetAtoms():
        assert atom.GetAtomMapNum() == atom.GetIdx()

    xymol.show_atom_number(zero_index=False)
    for atom in xymol.mol.GetAtoms():
        assert atom.GetAtomMapNum() == atom.GetIdx() + 1


def test_hide_atom_number():
    xymol = XYMOL(SMILES)
    xymol.show_atom_number()
    for atom in xymol.mol.GetAtoms():
        assert atom.GetAtomMapNum() == atom.GetIdx()

    xymol.hide_atom_number()
    for atom in xymol.mol.GetAtoms():
        assert atom.GetAtomMapNum() == 0
