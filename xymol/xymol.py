#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""XYMOL class
"""

from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops


class XYMOL:
    def __init__(self, smiles: str) -> None:
        self.smiles = smiles
        self.mol = self._build_mol()

    def _build_mol(self) -> Chem.rdchem.Mol:
        mol = Chem.MolFromSmiles(self.smiles)
        # reorder atoms
        order = rdmolfiles.CanonicalRankAtoms(mol)
        mol = rdmolops.RenumberAtoms(mol, order)
        return mol

    def _sanitize(self, mol) -> Chem.rdchem.Mol:
        while Chem.SanitizeMol(mol, catchErrors=True).real != 0:
            try:
                Chem.SanitizeMol(mol)
            except Exception as exc:
                atom_idx = int(exc.args[0].split(" ")[2])
                mol.GetAtomWithIdx(atom_idx).SetIsAromatic(False)
        return mol

    def drop_each_atom(self) -> list:
        mol = deepcopy(self.mol)
        smiles_list = []

        for drop_atom in range(len(mol.GetAtoms())):
            mol2 = deepcopy(mol)
            mol2.GetAtomWithIdx(drop_atom).SetAtomicNum(0)
            mol2 = Chem.DeleteSubstructs(mol2, Chem.MolFromSmarts('[#0]'))
            mol2 = self._sanitize(mol2)
            smiles_list.append(Chem.MolToSmiles(mol2))

        return smiles_list
