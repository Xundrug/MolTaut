import os
#os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=4
#os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=4
#os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=6
#os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=4
#os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=6

from openbabel import pybel
import numpy as np
import pandas as pd
from rdkit import Chem

from multiprocessing import Pool
from moltaut_src.utils import filter_mol
import torch
from moltaut_src.descriptor import mol2vec
from moltaut_src.models import load_model
from moltaut_src.optimize_mol import optimize
from moltaut_src.gen_confs import gen_confs_set
from moltaut_src.get_vmrs import enumerate_vmrs
import warnings
warnings.filterwarnings("ignore")

nmodel, imodel = load_model()

def linker_to_hydrogen(smi):
    mol = Chem.MolFromSmiles(smi)

    linker_aids = []
    for at in mol.GetAtoms():
        if at.GetSymbol() == "*":
            idx = at.GetIdx()
            linker_aids.append(idx)

    emol = Chem.RWMol(mol)
    for idx in linker_aids:
        emol.ReplaceAtom(idx, Chem.Atom(6))
    nmol = emol.GetMol()
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(nmol)))
    return smi

def predict_single(pmol, model, fmax):
    device = "cpu"
    obmol = pmol.OBMol
    if fmax:
        obmol, dE = optimize(obmol, fmax)
    else:
        dE = 0.0
    data = mol2vec(obmol)
    npmol = pybel.Molecule(obmol)
    with torch.no_grad():
        data = data.to(device)
        solv = model(data).cpu().numpy()[0][0]
    #npmol.write("mol", str(-1.0 * solv * 100000) + ".mol")
    return solv, dE

def predict_multicore_wrapper(param):
    idx, smi, block, fmax = param
    nmodel, imodel = load_model()
    pmol = pybel.readstring("mol", block)
    solv, dE = predict_single(pmol, nmodel, fmax)
    dE = dE / 27.2114 * 627.5094
    return [idx, smi, solv, dE, solv+dE]

def predict_by_smi(smi, fmax,  num_confs):
    pmol = pybel.readstring("smi", smi)
    blocks = gen_confs_set(smi, num_confs)
    params = zip(blocks, [fmax for i in range(len(blocks))])

    pool = Pool()
    score = pool.map(predict_multicore_wrapper, params)
    pool.close()
    
    score_sort = sorted(score, key=lambda x: x[2])
    solv, dE, dG = score_sort[0]
    return solv, dE 


def predict_by_smis(smis, fmax, num_confs):
    params = []
    for idx, smi in enumerate(smis):
        blocks = gen_confs_set(smi, num_confs)
        for block in blocks:
            params.append([idx, smi, block, fmax])

    pool = Pool()
    score = pool.map(predict_multicore_wrapper, params)
    pool.close()

    output = []
    df_score = pd.DataFrame(score)
    for smi_idx, dfsg in df_score.groupby(0):
        #print(smi_idx)
        dfsg_sorted = dfsg.sort_values(4)
        smi = dfsg_sorted.iloc[0, 1]
        solv = dfsg_sorted.iloc[0, 2]
        dE = dfsg_sorted.iloc[0, 3]
        output.append([smi_idx, smi, solv, dE])
    #print("output:", len(output))
    return output


def predict_by_mol(pmol, fmax, model):
    if not filter_mol( pmol ):
        print("#### Warning filter molecule")
        return
    solv, dE = predict_single(pmol, model, fmax)
    dE = dE / 27.2114 * 627.5094
    return solv, dE

def calc_solv(tauts, fmax, num_confs, is_fragment):
    if is_fragment:
        tauts_smis_include_linker = [Chem.MolToSmiles(taut.mol) for taut in tauts]
        tauts_smis_exclude_linker = [linker_to_hydrogen(smi) for smi in tauts_smis_include_linker]
        
        output = predict_by_smis(tauts_smis_exclude_linker, fmax, num_confs) 
        res = []
        for smi_idx, tsmi, solv, dE in output:
            lsmi = tauts_smis_include_linker[smi_idx]
            res.append([lsmi, solv, dE])
    else:
        tauts_smis = [taut.smi for taut in tauts]
        output = predict_by_smis(tauts_smis, fmax, num_confs)
        res = []
        for smi_idx, tsmi, solv, dE in output:
            res.append([tsmi, solv, dE])
    df = pd.DataFrame(res)
    if len(df) == 0:
        return df 
    df[3] = df[1] + df[2]
    df[3] = df[3] - df[3].min()
    df.columns = ["smi", "solv", "internal", "dG"] 
    return df

def rank_tauts(tauts, num_confs, fmax=0.01, is_fragment=True):
    df = calc_solv(tauts, fmax, num_confs, is_fragment)
    smirks_rules = [taut.smirks for taut in tauts]
    #print(len(tauts))
    #print(len(df))
    #print(len(smirks_rules))
    df["smirks"] = smirks_rules
    df["dG"] = df["dG"] * 0.72
    df = df.sort_values("dG")
    return df

        
if __name__=="__main__":
    #smi = "Brc1cnn2c1nc(cc2NCc1cccnc1)c1ccccc1"
    smi = "Cc1n[nH]c(c12)OC(N)=C(C#N)C2(C(C)C)c(cc3C(F)(F)F)cc(c3)N4CCCC4"
    #smi = "CS(=O)(=O)c1ccc(cc1)c1cccn2c1nc(n2)Nc1ccc(cc1)N1CCOCC1"
    vmrs = enumerate_vmrs(smi)
    print(vmrs)
    import sys
    index = int(sys.argv[1])
    tauts = vmrs[index].tauts
    df = rank_tauts(tauts)
    print(df)



