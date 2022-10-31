from rdkit.Chem import AllChem as Chem

def gen_confs_set(smi, num_confs):
    mol = Chem.MolFromSmiles(smi)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    cids = Chem.EmbedMultipleConfs(mol, num_confs, Chem.ETKDG())
    for conf in cids:
        converged =  Chem.MMFFOptimizeMolecule(mol,confId=conf)
        Chem.UFFOptimizeMolecule(mol,confId=conf)
    
    blocks = []
    for conf in cids:
        block = Chem.MolToMolBlock(mol, confId=conf)
        blocks.append(block)
    return blocks
