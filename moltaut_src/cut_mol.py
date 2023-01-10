from rdkit import Chem
import copy

bond_smarts = ["[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([NH,NH2,OH,SH]-[*;r]);!$(*=,#[*;!R])]", "[#6+0;!$(*=,#[!#6])]!@!=!#[#16;H0]"]
bond_patterns = [Chem.MolFromSmarts(bs) for bs in bond_smarts]

def get_bonds(mol):
    bond_matches = []
    for bp in bond_patterns:
        bms = mol.GetSubstructMatches(bp)
        bond_matches.extend(list(bms))
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
    mol = Chem.MolFromSmiles(smi)
    frags = get_frags(mol)
    print(frags)
