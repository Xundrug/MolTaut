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


def get_neighbor(obat):  
    nobats = [nobat for nobat in OBAtomAtomIter(obat)]
    if len(nobats) != 1:
        raise RuntimeError("Hydrogen neighbors is not one")
    obnat = nobats[0]
    return obnat

def get_hydrogen_dist(hat):
    obhat = hat.OBAtom
    obnat = get_neighbor(obhat)
    elem = obnat.GetAtomicNum()
    dist = obhat.GetDistance(obnat)
    hid = obhat.GetId()
    nid = obnat.GetId()
    return {"d": dist, "hid": hid, "nid": nid, "elem": elem}

def get_atypical_bond(pmol):
    atypical_ids = []
    for at in pmol.atoms:
        if at.atomicnum != 1:
            continue
        res = get_hydrogen_dist(at)
        if res["d"] > 1.2 and res["elem"] != 16:
            atypical_ids.append([res["hid"], res["nid"], res["elem"]])
        elif res["d"] > 1.4 and res["elem"] == 16:
            atypical_ids.append([res["hid"], res["nid"], res["elem"]])
    return atypical_ids


def set_constrain(asemol, atypical_bonds):
    """
    dihedral_indices: start from 0
    """
    #c = FixBondLengths(atypical_bonds)
    for id1, id2, elem in atypical_bonds:
        if elem == 16:
            c = Hookean(a1=id1, a2=id2, rt=1.4, k=7.0)
        else:
            c = Hookean(a1=id1, a2=id2, rt=1.2, k=5.0)
        asemol.set_constraint(c)
    return


def optimize_step(asemol, obmol, fmax):
    asemol.set_calculator(calculator)
    opt = BFGS(asemol, logfile="/dev/null")
    opt.run(fmax=fmax, steps=1)
    coords = asemol.get_positions()
    set_coords(obmol, coords)
    dE = asemol.get_potential_energy()
    return obmol, dE


def optimize(oobmol, fmax=0.05):
    """
    strict: fmax=0.001
    """
    opmol = pybel.Molecule(oobmol)
    pmol = opmol.clone
    obmol = pmol.OBMol

    asemol = pmol_to_asemol(pmol)
    obmol, dE = optimize_step(asemol, obmol, fmax)

    atypical_bonds = get_atypical_bond( pybel.Molecule(obmol) )
    if len(atypical_bonds) == 0:
        return obmol, dE
    else:
        print("-------------:", atypical_bonds)
        asemol = pmol_to_asemol(opmol)
        set_constrain(asemol, atypical_bonds)
        obmol, dE = optimize_step(asemol, opmol.OBMol, fmax)
        return obmol, dE


