#!/usr/bin/env python
# coding: utf-8


from ase.optimize import  BFGS
from ase import units
from ase.constraints import FixInternals, FixBondLengths, Hookean
from ase.atoms import Atoms
import torchani
from openbabel import pybel
from openbabel.openbabel import OBMolAtomIter, OBAtomAtomIter
from torchani.ase import Calculator
import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
ani_2x_model = torchani.models.ANI2x()
ani_2x_model = ani_2x_model.to(device)
calculator = ani_2x_model.ase()


def pmol_to_asemol(pmol):
    mol_block = pmol.write("mol")
    lines = mol_block.split("\n")
    L1 = lines[3]

    # The V2000 dialect uses a fixed field length of 3, which means there
    # won't be space between the numbers if there are 100+ atoms, and
    # the format doesn't support 1000+ atoms at all.
    if L1.rstrip().endswith('V2000'):
        natoms = int(L1[:3].strip())
    else:
        natoms = int(L1.split()[0])
    positions = []
    symbols = []
    for line in lines[4:4 + natoms]:
        x, y, z, symbol = line.split()[:4]
        symbols.append(symbol)
        positions.append([float(x), float(y), float(z)])
    return Atoms(symbols=symbols, positions=positions)


def set_coords(obmol, new_coords):
    i = 0
    for at in OBMolAtomIter(obmol):
        at_coords = new_coords[i, :]
        x = at_coords[0]
        y = at_coords[1]
        z = at_coords[2]
        at.SetVector(x, y,  z)
        i += 1
    return


def get_bond_length_dict(obmol):
    num_bonds = obmol.NumBonds()
    bonds = [obmol.GetBond(i) for i in range(num_bonds)]
    
    bond_length_dict = {}
    for i in range(num_bonds):
        bond_length_dict[i] = bonds[i].GetLength()
    return bond_length_dict

def get_atypical_bond(now_obmol, origin_obmol):
    atypical_ids = []
    origin_length_dict = get_bond_length_dict(origin_obmol)
    now_length_dict = get_bond_length_dict(now_obmol)
    for key, origin_length in origin_length_dict.items():
        now_length = now_length_dict[key]
        if abs(now_length - origin_length) > 0.25:
            #print(now_length)
            #print(origin_length_dict)
            #print(now_length_dict)
            erbd = now_obmol.GetBond(key)
            atypical_ids.append([erbd.GetBeginAtomIdx()-1, erbd.GetEndAtomIdx()-1, origin_length])
    return atypical_ids

def set_constrain(asemol, atypical_bonds):
    for id1, id2, bd_length in atypical_bonds:
        c = Hookean(a1=id1, a2=id2, rt=bd_length, k=30.0)
        asemol.set_constraint(c)
    #bds = [[bdinf[2], [bdinf[0], bdinf[1]]] for bdinf in atypical_bonds]
    #c = FixInternals(bonds=bds)
    #asemol.set_constraint(c)
    return


def optimize_step(asemol, obmol, fmax, steps):
    asemol.set_calculator(calculator)
    opt = BFGS(asemol, logfile="/dev/null")
    opt.run(fmax=fmax, steps=steps)
    coords = asemol.get_positions()
    set_coords(obmol, coords)
    dE = asemol.get_potential_energy()
    return obmol, dE


def optimize(oobmol, fmax=0.05, steps=24):
    """
    strict: fmax=0.001
    """
    opmol = pybel.Molecule(oobmol)


    pmol1 = opmol.clone
    obmol1 = pmol1.OBMol
    pmol2 = opmol.clone
    obmol2 = pmol2.OBMol

    asemol = pmol_to_asemol(pmol1)
    obmol1, dE = optimize_step(asemol, obmol1, fmax, steps=steps)
  
    atypical_bonds = get_atypical_bond( obmol1, oobmol )
    if len(atypical_bonds) == 0:
        return obmol1, dE
    else:
        #print("-------------:", atypical_bonds)
        asemol = pmol_to_asemol(pmol2)
        #set_constrain(asemol, atypical_bonds)
        obmol2, dE = optimize_step(asemol, obmol2, fmax, steps=steps)
        #print(dE, atypical_bonds)
        #print(get_bond_length_dict(obmol2))
        dE = 10000.0
        return obmol2, dE


