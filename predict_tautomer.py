import warnings
warnings.filterwarnings("ignore")
from itertools import product
import random
import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from rdkit.Chem import AllChem
from moltaut_src.tautomer import enumerate_tauts
from moltaut_src.combine_frag import link_fragment
from moltaut_src.rank_tautomer import rank_tauts
from moltaut_src.molgpka.protonate import protonate_mol
from moltaut_src.get_vmrs import enumerate_vmrs
from collections import namedtuple
import os
import argparse
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

un = rdMolStandardize.Uncharger()


def is_need_mol(mol, element_list=[1, 6, 7, 8, 9, 15, 16, 17]):
    if mol is not None:
        elements = all(
            [at.GetAtomicNum() in element_list for at in mol.GetAtoms()])
        if elements:
            return True
        else:
            return False


def get_lower_energy_tauts(smi, energy_range, num_confs):
    vmrs = enumerate_vmrs(smi)

    data = namedtuple("lowerEnergyTauts", "smi smirks_index energy lower")

    lower_energy_tauts = []
    for vmr in vmrs:
        tauts = vmr.tauts
        if len(tauts) == 1:
            lower_energy_tauts.append([
                data(
                    smi=vmr.smi,
                    smirks_index=-1,
                    energy=0.0,
                    lower=True
                )])
        else:
            df_score = rank_tauts(tauts, num_confs)
            conts = []
            for idx, row in df_score.iterrows():
                smirks_index = row[4][0]
                taut_smi = row[0]
                energy = row[3]
                if energy <= energy_range:
                    conts.append(
                        data(
                            smi=taut_smi,
                            smirks_index=smirks_index,
                            energy=energy,
                            lower=True
                        ))
                else:
                    conts.append(
                        data(
                            smi=taut_smi,
                            smirks_index=smirks_index,
                            energy=energy,
                            lower=False
                        ))
            lower_energy_tauts.append(conts)
    return lower_energy_tauts


def combine_lower_energy_tauts(lower_energy_tauts):
    tauts_product = list(product(*lower_energy_tauts))
    lower_energy_mols, upper_energy_mols = [], []
    for tauts in tauts_product:
        smis, energies, labels = [], [], []
        for taut in tauts:
            smis.append(taut.smi)
            energies.append(taut.energy)
            labels.append(taut.lower)
        dG = sum(energies)
        m = link_fragment(smis)
        if all(labels):
            lower_energy_mols.append([Chem.MolToSmiles(m), dG])
        else:
            upper_energy_mols.append([Chem.MolToSmiles(m), dG])
    return lower_energy_mols, upper_energy_mols


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


def is_cut_mol(mm):
    bonds_idxs = match_bonds(mm)
    atom_idxs = match_atoms(mm)

    filter_bond_idxs = []
    for bond_idx in bonds_idxs:
        begin_idx = bond_idx[0]
        end_idx = bond_idx[1]
        if (begin_idx in atom_idxs) or (end_idx in atom_idxs):
            continue
        filter_bond_idxs.append(bond_idx)
    if len(filter_bond_idxs) == 0:
        return False
    else:
        return True


def generate_tautomer_cutmol(smi, num_confs, energy_range, ph, tph):
    lower_energy_tauts = get_lower_energy_tauts(
        smi,
        energy_range,
        num_confs)
    lower_energy_mols, upper_energy_mols = combine_lower_energy_tauts(
        lower_energy_tauts)
    df_res_lower = pd.DataFrame(lower_energy_mols)
    dfs_res_lower = df_res_lower.sort_values(1)

    if len(upper_energy_mols) == 0:
        dfs_res_upper = pd.DataFrame({0: [], 1: [], 2: []})
    else:
        dfs_res_upper = pd.DataFrame(upper_energy_mols)
        dfs_res_upper = dfs_res_upper.sort_values(1)
        dfs_res_upper[2] = dfs_res_upper[0]
    return dfs_res_lower, dfs_res_upper


def generate_tautomer_non_cutmol(mm, num_confs, energy_range, ph, tph):
    tauts = enumerate_tauts(mm)
    df_res = rank_tauts(tauts, num_confs, is_fragment=False)
    df_res = df_res.iloc[:, [0, 3]]
    df_res.columns = [0, 1]

    dfs_res_lower = df_res[df_res[1] <= energy_range].copy()
    dfs_res_lower = dfs_res_lower.sort_values(1)
    dfs_res_upper = df_res[df_res[1] > energy_range].copy()
    if len(dfs_res_upper) == 0:
        dfs_res_upper = pd.DataFrame({0: [], 1: [], 2: []})
    else:
        dfs_res_upper = dfs_res_upper.sort_values(1)
        dfs_res_upper[2] = dfs_res_upper[0]
    return dfs_res_lower, dfs_res_upper


def func(smi, cutmol, energy_range=2.8, ph=7.0, tph=1.0, num_confs=3):
    mm = Chem.MolFromSmiles(smi)
    mm = un.uncharge(mm)
    mm = Chem.MolFromSmiles(Chem.MolToSmiles(mm))
    if cutmol:
        if is_cut_mol(mm):
            dfs_res_lower, dfs_res_upper = generate_tautomer_cutmol(
                smi, energy_range=energy_range, num_confs=num_confs, ph=ph, tph=tph)
        else:
            dfs_res_lower, dfs_res_upper = generate_tautomer_non_cutmol(
                mm, energy_range=energy_range, num_confs=num_confs, ph=ph, tph=tph)
    else:
        dfs_res_lower, dfs_res_upper = generate_tautomer_non_cutmol(
            mm, energy_range=energy_range, num_confs=num_confs, ph=ph, tph=tph)
    dfs_res_lower[2] = dfs_res_lower[0].map(
        lambda x: protonate_mol(x, ph, tph))
    return dfs_res_lower, dfs_res_upper


def generate_conf(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    cids = AllChem.EmbedMultipleConfs(mol, 1, AllChem.ETKDG())
    for conf in cids:
        converged = AllChem.MMFFOptimizeMolecule(mol, confId=conf)
        AllChem.UFFOptimizeMolecule(mol, confId=conf)
    return mol, cids

def write_file(datas, sdf_path):
    conf_data = []
    for data in datas:
        tsmi = data['tsmi']
        psmis = data['psmis']
        score = data['score']
        label = data['label']
        if label == "high_energy":
            continue
        for smi in psmis:
            mol, cids = generate_conf(smi)
            mol.SetProp("tautomer smiles", tsmi)
            mol.SetProp("protonation smiles", smi)
            mol.SetProp("Score", score)
            mol.SetProp("Label", label)
            conf_data.append([mol, cids])

    sdw = Chem.SDWriter(sdf_path)
    for mol, cids in conf_data:
        for cid in cids:
            sdw.write(mol, confId=cid)
    sdw.close()
    return


def construct_data(dfs, label):
    datas = []
    for idx, row in dfs.iterrows():
        tsmi = row[0]
        score = row[1]
        psmis = row[2]

        data = {}
        data['tsmi'] = tsmi
        data['psmis'] = psmis
        data['score'] = str(round(score, 2))
        data['label'] = label
        datas.append(data)
    return datas


def get_taut_data(smi, cutmol, num_confs, energy_cutoff, ph, tph):
    dfs_res_lower, dfs_res_upper = func(
        smi,
        cutmol=cutmol,
        energy_range=energy_cutoff,
        num_confs=num_confs,
        ph=ph,
        tph=tph)
    datas_lower = construct_data(
        dfs_res_lower,
        label="low_energy")
    datas_upper = construct_data(
        dfs_res_upper,
        label="high_energy")
    fdatas = datas_lower + datas_upper
    return fdatas


def run():
    parser = argparse.ArgumentParser(
        description='calculate low-energy tautomer for small molecules')
    parser.add_argument(
        '--smi',
        type=str,
        default='O=c1ccnc[nH]1',
        help='the molecular smiles')
    parser.add_argument(
        '--low_energy_tautomer_cutoff',
        type=float,
        default=2.8,
        help='the energy cutoff for low energy')
    parser.add_argument('--cutmol', type=int, default=True,
                        help='determine to frag the molecule')
    parser.add_argument(
        '--num_confs',
        type=int,
        default=50,
        help='the number of conformation for solvation energy prediction')
    parser.add_argument(
        '--ph',
        type=float,
        default=7.0,
        help='the target pH for protonation states generation')
    parser.add_argument(
        '--tph',
        type=float,
        default=1.0,
        help='pH tolerance for protonation states generation')
    parser.add_argument(
        '--output',
        type=str,
        default="moltaut_output.sdf",
        help='the output SDF file name')

    args = parser.parse_args()

    smi = args.smi
    energy_cutoff = args.low_energy_tautomer_cutoff
    cutmol = args.cutmol
    num_confs = args.num_confs
    ph = args.ph
    tph = args.tph
    output = args.output

    data = get_taut_data(smi, cutmol, num_confs, energy_cutoff, ph, tph)
    print(data)
    write_file(data, output)
    return


if __name__ == "__main__":
    run()
