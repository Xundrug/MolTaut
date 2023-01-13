from rdkit import Chem
import copy


def match_bonds(mm):
    tsmarts = ["[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!X1;!$([NH,NH2,OH,SH]-[*;r]);!$(*=,#[*;!R])]"]
    tpatterns = [Chem.MolFromSmarts(tsm) for tsm in tsmarts]
    matches = []
    for tpat in tpatterns:
        tms = mm.GetSubstructMatches(tpat)
        matches.extend(list(tms))
    return matches


def match_atoms(mm):
    fsmarts = ["[$([#6]([F,Cl])-[*;r])]"]
    fpatterns = [Chem.MolFromSmarts(fsm) for fsm in fsmarts]
    fatom_idxs = []
    for fpat in fpatterns:
        fms = mm.GetSubstructMatches(fpat)
        fatom_idxs.extend(list(fms))
    fatom_idxs = sum(fatom_idxs, ())
    return fatom_idxs


def get_bond_idxs(mm):
    bonds_idxs = match_bonds(mm)
    atom_idxs = match_atoms(mm)
    
    filter_bond_idxs = []
    for bond_idx in bonds_idxs:
        begin_idx = bond_idx[0]
        end_idx = bond_idx[1]
        if (begin_idx in atom_idxs) or (end_idx in atom_idxs):
            continue
        filter_bond_idxs.append(bond_idx)
    return filter_bond_idxs


def get_bonds(mol):
    bond_matches = get_bond_idxs(mol)

    bonds = []
    for begin_idx, end_idx in bond_matches:
        bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
        bonds.append(bond.GetIdx())
    return bonds

def get_frags(mol):
    bonds = get_bonds(mol)
    if len(bonds) == 0:
        return []
    m = copy.deepcopy(mol)
    labels = [[i,i] for i in range(1,len(bonds)+1)]
    nm = Chem.FragmentOnBonds(m, bonds, addDummies=True, dummyLabels=labels)
    frag_mols = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=True)

    frag_smiles = []
    for fmol in frag_mols:
        fsmi = Chem.MolToSmiles(fmol)
        for i in range(1, len(bonds)+1):
            fsmi = fsmi.replace(f'[{i}*]',f'[*:{i}]')
        frag_smiles.append(fsmi)
    return frag_smiles


if __name__ == "__main__":
    smi = "CN(C)CCCN1C2=CC=CC=C2SC2=C1C=C(C=C2)C(C)=O"
    smi = "Brc1cnn2c1nc(cc2NCc1cccnc1)c1ccccc1"
    smi = "Cc1n[nH]c(c12)OC(N)=C(C#N)C2(C(C)C)c(cc3C(F)(F)F)cc(c3)N4CCCC4"
    smi = "Nc1nc2c([nH]1)cccn2"
    smi = "c1ncccc1-c(n2)[nH]c(c23)CCCc4c3cc(F)cc4"
    smi = "Cc1c2c([nH]n1)OC(=C([C@@]2(c3cc(cc(c3)N4CCCC4)C(F)(F)F)C(C)C)C#N)N"
    mol = Chem.MolFromSmiles(smi)
    frags = get_frags(mol)
    print(frags)
