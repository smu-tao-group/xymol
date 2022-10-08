#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""XYMOL class
"""

from typing import Tuple
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import Draw, rdmolfiles, rdmolops
from rdkit.Chem.Draw import SimilarityMaps
import numpy as np


class XYMOL:
    """XYMOL class
    """
    def __init__(self, smiles: str) -> None:
        """Init class

        Args:
            smiles (str): SMILES representation.
        """
        self.smiles = smiles
        self.mol = self._build_mol()

    def _build_mol(self) -> Chem.rdchem.Mol:
        """Build molecule

        Returns:
            Chem.rdchem.Mol: rdkit molecule.
        """
        mol = Chem.MolFromSmiles(self.smiles)
        # reorder atoms
        order = rdmolfiles.CanonicalRankAtoms(mol)
        mol = rdmolops.RenumberAtoms(mol, order)
        return mol

    def _sanitize(self, mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
        """Sanitize molecule

        Args:
            mol (Chem.rdchem.Mol): rdkit molecule.

        Returns:
            Chem.rdchem.Mol: sanitized molecule.
        """
        while Chem.SanitizeMol(mol, catchErrors=True).real != 0:
            try:
                Chem.SanitizeMol(mol)
            except Exception as exc:
                atom_idx = int(exc.args[0].split(" ")[2])
                mol.GetAtomWithIdx(atom_idx).SetIsAromatic(False)
        return mol

    def drop_each_atom(self) -> list:
        """Drop each atom, one at a time

        Returns:
            list: a list of SMILES.
        """
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
        """Drop each bond, one at a time

        Returns:
            Tuple[list, list]: a list of bonds, a list of SMILES.
        """
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
        """Show atom number

        Args:
            prop_label (str, optional): Defaults to 'molAtomMapNumber'.
            zero_index (bool, optional): Atom index start from 0.
        """
        for atom in self.mol.GetAtoms():
            index = atom.GetIdx() + (zero_index is False)
            atom.SetProp(prop_label, str(index))

    def hide_atom_number(self) -> None:
        """Hide atom number
        """
        for atom in self.mol.GetAtoms():
            atom.SetAtomMapNum(0)

    def plot_mol(self, file_name="mol.png", **kwargs) -> None:
        """Plot molecule structure

        Args:
            file_name (str, optional): img file name. Defaults to "mol.png".
        """
        img = Draw.MolToImage(self.mol, **kwargs)
        img.save(file_name)

    def plot_similarity_map(
        self, weights: list, file_name="similarity_map.png", **kwargs
    ) -> None:
        """Plot similarity map

        Args:
            weights (list): weights by each atom.
            file_name (str, optional): Defaults to "similarity_map.png".
        """
        img = SimilarityMaps.GetSimilarityMapFromWeights(
            self.mol, weights, **kwargs
        )
        img.savefig(file_name, bbox_inches='tight', dpi=600)

    def create_map(
        self, featurizer, model, filename="similarity_map.png"
        ) -> None:
        """Create similarity map

        Args:
            featurizer (object): must have a featurize() function.
            model (object): must have a predict() function.
        """
        if not(
            hasattr(featurizer, "featurize")
            and callable( getattr(featurizer, "featurize") )
        ):
            raise ValueError("Featurizer must have featurize() function.")

        if not(
            hasattr(model, "predict")
            and callable( getattr(model, "predict") )
        ):
            raise ValueError("Model must have predict() function.")


        parent_feat = np.array(featurizer.featurize(self.smiles))
        parent_prediction = model.predict(parent_feat)[0]

        drop_atoms = np.array(self.drop_each_atom())
        drop_atoms_feat = np.array(featurizer.featurize(drop_atoms))

        predictions = model.predict(drop_atoms_feat)

        weights = [pred - parent_prediction for pred in predictions]
        self.plot_similarity_map(weights, filename)
