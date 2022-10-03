#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""XYMOL class
"""

from typing import Tuple
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import Draw, rdmolfiles, rdmolops
from rdkit.Chem.Draw import SimilarityMaps


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

    def drop_each_bond(self) -> Tuple[list, list]:
        mol = deepcopy(self.mol)
        rwmol = Chem.RWMol(mol)
        smiles = []
        bonds = []

        for bond in mol.GetBonds():
            rwmol_copy = deepcopy(rwmol)
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            rwmol_copy.RemoveBond(begin_atom, end_atom)
            rwmol_copy = self._sanitize(rwmol_copy)
            bonds.append([begin_atom, end_atom])
            smiles.append(Chem.MolToSmiles(rwmol_copy))

        return bonds, smiles

    def show_atom_number(
        self, prop_label='molAtomMapNumber', zero_index=True
    ) -> None:
        for atom in self.mol.GetAtoms():
            index = atom.GetIdx() + (zero_index is False)
            atom.SetProp(prop_label, str(index))

    def hide_atom_number(self) -> None:
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(0)

    def plot_mol(self, file_name="mol.png", **kwargs) -> None:
        img = Draw.MolToImage(self.mol, **kwargs)
        img.save(file_name)

    def plot_similarity_map(
        self, weights, file_name="similarity_map.png", **kwargs
    ) -> None:
        img = SimilarityMaps.GetSimilarityMapFromWeights(
            self.mol, weights, **kwargs
        )
        img.savefig(file_name, bbox_inches='tight', dpi=600)
