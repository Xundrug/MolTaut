#!/usr/bin/env python
# coding: utf-8
import numpy as np
from scipy.spatial.distance import euclidean, cdist

from openbabel import pybel
from openbabel.openbabel import OBMolAtomIter, OBMolBondIter
from moltaut_src.aevs import calc_descriptors
from moltaut_src.sasa import get_sasa
import pickle
import torch
from torch_geometric.data import Data


def get_atom_coords(at):
    acoords = np.array([at.GetX(), at.GetY(), at.GetZ()])
    return acoords


def get_mol_coords(obmol):
    mol_coords = []
    for at in OBMolAtomIter(obmol):
        coords = get_atom_coords( at )
        mol_coords.append(coords)
    return np.array(mol_coords)


def get_cdist(obmol):
    mol_coords = get_mol_coords(obmol)
    dmatrix = cdist(mol_coords, mol_coords)
    return dmatrix


def one_hot(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(
                    x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))


def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def get_atom_features(obmol):
    dmatrix = get_cdist(obmol)
    pmol = pybel.Molecule(obmol)
    sasas = get_sasa(obmol)
    m = []
    for at in OBMolAtomIter(obmol):
        idx = at.GetIdx()
        o = []
        o += one_hot(at.GetAtomicNum(), [1, 6, 7, 8, 16, 9, 17])
        o += [at.IsAromatic()]
        o += [at.IsInRingSize(3),
              at.IsInRingSize(4),
              at.IsInRingSize(5),
              at.IsInRingSize(6),
              at.IsInRingSize(7),
              at.IsInRingSize(8)]
        o += [at.GetImplicitHCount()]
        o += [at.GetFormalCharge()]
        o += [at.GetHyb()]
        o += calc_descriptors(pmol, idx, dmatrix)
        o += [sasas[int(idx)]]
        m.append(o)
    return m

def get_bond_pair(obmol):
    res = [[],[]]
    for bond in OBMolBondIter(obmol):
        res[0] += [bond.GetBeginAtomIdx()-1, bond.GetEndAtomIdx()-1]
        res[1] += [bond.GetEndAtomIdx()-1, bond.GetBeginAtomIdx()-1]
    return res

def mol2vec(obmol):
    node_f= get_atom_features(obmol)
    edge_index = get_bond_pair(obmol)
    
    batch = np.zeros(len(node_f), )
    data = Data(x=torch.tensor(node_f, dtype=torch.float32),
                edge_index=torch.tensor(edge_index, dtype=torch.long),
                batch=torch.tensor(batch, dtype=torch.long))
    return data



