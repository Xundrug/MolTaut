#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import sys
import glob
import pickle
from rdkit import Chem

from multiprocessing import Pool
import argparse


def parse_log(filename):
    molfile = filename.replace("log", "mol")
    mol = Chem.MolFromMolFile(molfile, removeHs=False)

    cmd = "grep 'DeltaG (solv)' {}".format(filename)
    r = os.popen(cmd)
    conts = r.read()
    r.close()
    
    solvdG = float(conts.strip().strip("\n").split()[4])
    conf_data = {
        "mol": None,
        "solvdG": None,
        "name": None}

    conf_data["solvdG"] = solvdG
    conf_data["mol"] = mol
    conf_data["name"] = os.path.join(
        os.path.basename(
            os.path.dirname(filename)),
        os.path.basename(filename))
    return conf_data


def collect_logs(index):
    logs_file = logs_file_chunks[index]

    smd_datas = []
    for frag_id, filename in enumerate(logs_file):
        try:
            data = parse_log(filename)
            smd_datas.append(data)
        except:
            print(filename)
    with open("tmp/scan_log_{}.pickle".format(index), "wb") as f:
        f.write(pickle.dumps(smd_datas))
    return


def split_list(a, n):
    k, m = divmod(len(a), n)
    ko = list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)]
              for i in range(n))
    return ko


def merge_pickle(cores, pickle_out="logs_all.pickle"):
    pickle_files = [
        "tmp/scan_log_{}.pickle".format(idx) for idx in range(cores)]

    data = []
    for pfile in pickle_files:
        with open(pfile, "rb") as f:
            conts = pickle.load(f)
        data.extend(conts)

    with open(pickle_out, "wb") as f:
        f.write(pickle.dumps(data))
    return


if __name__ == "__main__":
    os.system("rm -r tmp")

    if not os.path.exists("tmp"):
        os.makedirs("tmp")

    parser = argparse.ArgumentParser(description='parse gaussian scan logs')
    parser.add_argument('--cores', type=int, help='number of cpu')
    parser.add_argument(
        '--logs',
        type=str,
        help='the dirname of folder for scan logs')
    parser.add_argument(
        '--out',
        type=str,
        help='the file name of pickle for output')
    args = parser.parse_args()

    cores = args.cores
    pickle_out = args.out
    scan_dirname = args.logs

    all_logs_file = glob.glob(os.path.join(scan_dirname, "*", "*.log"))
    print("log", all_logs_file[0])
    logs_file_chunks = split_list(all_logs_file, cores)

    pool = Pool(cores)
    pool.map(collect_logs, list(range(cores)))
    pool.close()

    merge_pickle(cores, pickle_out)
