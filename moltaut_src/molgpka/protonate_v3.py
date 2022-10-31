from src.molgpka.predict_pka import predict
from copy import deepcopy
from rdkit import Chem

from rdkit import Chem
from rdkit.Chem import AllChem,Draw
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

import json
import numpy as np
import random
import os
import copy

def modify_mol(mol, acid_dict, base_dict):
    for at in mol.GetAtoms():
        idx = at.GetIdx()
        if idx in set(acid_dict.keys()):
            value = acid_dict[idx]
            nat = at.GetNeighbors()[0]
            nat.SetProp("ionization", "A")
            nat.SetProp("pKa", str(value))
        elif idx in set(base_dict.keys()):
            value = base_dict[idx]
            at.SetProp("ionization", "B")
            at.SetProp("pKa", str(value))
        else:
            at.SetProp("ionization", "O")
    nmol = AllChem.RemoveHs(mol)
    return nmol

def refined_idx(mol, oacid_dict, obase_dict):
    mc = modify_mol(mol, oacid_dict, obase_dict)
    acids, bases = [], []
    for at in mc.GetAtoms():
        props = at.GetPropsAsDict()
        acid_or_basic = props.get('ionization', False)
        pKa = float(props.get('pKa', False))
        idx = at.GetIdx()
        if acid_or_basic == "A":
            acids.append( [idx, pKa] )
        elif acid_or_basic == "B":
            bases.append( [idx, pKa] )
        else:
            continue
    return bases, acids, mc

def protonate_mol(smi, ph, tph):
    omol = Chem.MolFromSmiles(smi)
    obase_dict, oacid_dict, omol = predict(omol)
    mc = modify_mol(omol, oacid_dict, obase_dict)
    for at in mc.GetAtoms():
        props = at.GetPropsAsDict()
        acid_or_basic = props.get('ionization', False)
        pKa = float(props.get('pKa', False))
        idx = at.GetIdx()
        if acid_or_basic == "A" and pKa <= ph-tph:
            hnum = at.GetNumExplicitHs()
            at.SetFormalCharge(-1)
            at.SetNumExplicitHs(hnum-1)
        elif acid_or_basic == "B" and pKa >= ph+tph:
            hnum = at.GetNumExplicitHs()
            at.SetFormalCharge(1)
            at.SetNumExplicitHs(hnum+1)
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mc)))
    return smi

if __name__=="__main__":
    mol = Chem.MolFromSmiles("CN(C)CCCN1C2=CC=CC=C2SC2=C1C=C(C=C2)C(C)=O")
    pt_smis = protonate_mol(mol, ph=7.0, tph=0.5)
    print(pt_smis)


   

