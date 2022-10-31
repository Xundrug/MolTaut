#!/usr/bin/env python
# coding: utf-8

from openbabel import pybel
import numpy as np
from scipy.spatial.distance import euclidean, cdist
from openbabel.openbabel import OBMolAtomIter, OBAtomAtomIter
from itertools import combinations
import json
import pickle

# idx : start from 1


def get_atom_coords(at):
    acoords = np.array([at.GetX(), at.GetY(), at.GetZ()])
    return acoords


def get_atominc_num(at):
    if at.GetFormalCharge() < 0:
        return -2
    elif at.GetFormalCharge() > 0:
        return 2
    else:
        return at.GetAtomicNum()


def cutoff_func(Rc, Rij):
    if Rij <= Rc:
        return 0.5 * (np.cos(np.pi * Rij / Rc) + 1)
    else:
        return 0


#def get_distance(at1, at2):
#    acoord1 = get_atom_coords(at1)
#    acoord2 = get_atom_coords(at2)
#    d = euclidean(acoord1, acoord2)
#    return d

def get_distance(acoord1, acoord2):
    d = euclidean(acoord1, acoord2)
    return d


def get_cdist(obmol):
    mol_coords = get_mol_coords(obmol)
    dmatrix = cdist(mol_coords, mol_coords)
    return dmatrix


def get_mol_coords(obmol):
    mol_coords = [] 
    for at in OBMolAtomIter(obmol):
        coords = get_atom_coords( at )
        mol_coords.append(mol_coords)
    return np.array(mol_coords)
        

class CalcRadialSymmetry(object):
    def __init__(self, pmol, idx, dmatrix):
        self.obmol = pmol.clone.OBMol
        self.center_atom = self.obmol.GetAtom(idx) # idx start from 1
        self.cidx = idx
        self.dmatrix = dmatrix
        #self.cutoff_radiis = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 10.0]
        self.cutoff_radiis = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.5]
        self.atom_types = [1, 6, 7, 8, 16, 9, 17, -2, 2]  # -2->neg +2->pos
        self.eta = 0.0001
        self.Rs = 0.0

    def calc_single_symmetry(self, atom_type, Rc):
        G = 0
        for at in OBMolAtomIter(self.obmol):
            atomicnum = get_atominc_num(at)
            if atomicnum == atom_type:
                #Rij = get_distance(at, self.center_atom)
                #acoords = get_atom_coords( at )
                #Rij = get_distance( acoords, self.center_coords )
                Rij = self.dmatrix[at.GetIdx()-1, self.cidx-1]
                Gj = np.exp(-1.0 * self.eta * np.power(Rij -
                                                       self.Rs, 2)) * cutoff_func(Rc, Rij)
                G += Gj
        return G

    def calc_symmetry(self):
        symmetry = []
        for atom_type in self.atom_types:
            for Rc in self.cutoff_radiis:
                G = self.calc_single_symmetry(atom_type, Rc)
                symmetry.append(G)
        return symmetry


class CalcAngularSymmetry(object):
    def __init__(self, pmol, idx, dmatrix):
        self.obmol = pmol.clone.OBMol
        self.catom = self.obmol.GetAtom(idx)
        self.cidx = idx
        #self.ccoords = get_atom_coords(self.catom)
        self.dmatrix = dmatrix

        self.atom_types = [1, 6, 7, 8, 16, 9, 17, -2, 2]  # -2->neg +2->pos
        #self.cutoff_radiis = [0.8, 1.5, 2.0, 2.5, 3.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 10.0]
        #self.cutoff_radiis = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 10.0]
        self.cutoff_radiis =  [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.5]
        self.atoms = self.get_atoms()
        self.peta = 0.0001
        self.pzeta = 0.5
        self.plambda = 0.5
        self.atomtype_combinations = self.atom_types_combination()

    def get_atoms(self):
        atoms = []
        for at in OBMolAtomIter(self.obmol):
            if at.GetIdx() == self.catom.GetIdx():
                continue
            atoms.append(at)
        return atoms

    def atom_types_combination(self):
        atomtype_combinations = list(combinations(self.atom_types, 2))
        for atom_type in self.atom_types:
            atomtype_combinations.append((atom_type, atom_type))
        return atomtype_combinations

    def calc_single_symmetry(self, atom_type_j, atom_type_k, Rc):
        G = 0
        for atj in self.atoms:
            if atj.GetAtomicNum() != atom_type_j:
                continue
            for atk in self.atoms:
                if atk.GetAtomicNum() != atom_type_k:
                    continue
                if atj.GetIdx() == atk.GetIdx():
                    continue
                #Rij = get_distance(self.catom, atj)
                #Rik = get_distance(self.catom, atk)
                #Rjk = get_distance(atj, atk)
                #atj_coords = get_atom_coords( atj )
                #atk_coords = get_atom_coords( atk )
                #Rij = get_distance(self.ccoords, atj_coords) 
                #Rik = get_distance(self.ccoords, atk_coords)
                #Rjk = get_distance(atj_coords, atk_coords)
                Rij = self.dmatrix[self.cidx-1, atj.GetIdx()-1]
                Rik = self.dmatrix[self.cidx-1, atk.GetIdx()-1]
                Rjk = self.dmatrix[atj.GetIdx()-1, atk.GetIdx()-1]
                theta = self.obmol.GetAngle(atk, self.catom, atj)
                fij = cutoff_func(Rc, Rij)
                fik = cutoff_func(Rc, Rik)
                fjk = cutoff_func(Rc, Rjk)
                Gijk = np.power(1 + self.plambda * np.cos(theta), self.pzeta) * \
                    np.exp(-1.0 * self.peta * (Rij**2 + Rik**2 + Rjk**2)) * fij * fik * fjk
                G += Gijk
        return np.power(2, 1 - self.pzeta)* G

    def calc_symmetry(self):
        symmetry = []
        for atom_type_j, atom_type_k in self.atomtype_combinations:
            for Rc in self.cutoff_radiis:
                G = self.calc_single_symmetry(atom_type_j, atom_type_k, Rc)
                symmetry.append(G)
        return symmetry

def calc_descriptors( pmol, idx, dmatrix):
    obmol = pmol.OBMol
    #dmatrix = get_cdist(obmol)
    radial_symmetry = CalcRadialSymmetry(pmol, idx, dmatrix)
    angular_symmetry = CalcAngularSymmetry(pmol, idx, dmatrix)
    radial_fp = radial_symmetry.calc_symmetry()
    angular_fp = angular_symmetry.calc_symmetry()
    
    ani_descriptors = radial_fp + angular_fp
    return ani_descriptors


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='calculate ani descriptor for torsion energy')
    parser.add_argument('--input', type=str, help='the sdf file')
    parser.add_argument('--output', type=str, help='the output pickle file')
    args = parser.parse_args()

    infile = args.input
    pickle_out = args.output
    func(infile, pickle_out)
