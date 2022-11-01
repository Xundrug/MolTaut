# MolTaut: A Tool for Rapid Generation of Favorable Tautomer in Aqueous Solution

Fast and proper treatment of tautomeric state for drug-like molecules is critical in computer-aided drug discovery since the major tautomer of a molecule determines its pharmacophore features and physical properties. We present MolTaut, a tool for the rapid generation of favorable states of drug-like molecules in water. MolTaut works by enumerating possible tautomer states with tautomeric transformation rules, ranking tautomers with their relative internal energies and solvation energies calculated by AI-based models, and generating preferred ionization states according to predicted microscopic pKa. We found that the ranking ability of the AI-based tautomer scoring approach is comparable to the DFT method from which the AI models try to learn. To facilitate the usage of MolTaut, we made a web server, which is available at http://moltaut.xundrug.cn. MolTaut is a handy tool for generating low energy tatuomers when preparing chemical structures in computer-aided drug discovery. Besides, MolTaut can be used for analyzing substitution effect on the tautomeric equilibrium in lead optimization.


![image](https://github.com/Xundrug/MolTaut/blob/master/images/Picture1.png)


## Requirements

* Python 3.6
* openbabel >= 3.0
* numpy 1.18.1
* RDKit 2020.09.1.0
* scipy
* pandas 0.25.3
* freesasa
* pytorch 1.10.2
* pytorch geometric 2.0.3
* torch-scatter 2.0.9 
* torch-sparse 0.6.12
* torchani 2.2   
* ase 3.22.1

You also can create the python environment by conda configure file:
```
conda env create -f environment.yaml
```
If you run torch-sparse with error, please uninstall the package `torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric`:
```
pip uninstall torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric
```
and then reinstall them:
```
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.10.0+cpu.html
```

Or you can use the environment created by conda-pack, and activate the python env by conda. Some errors may occur when you do `from openbabel import pybel`, you just need to reinstall openbabel by conda. The environment file download URL is as fellow:  
```
https://drive.google.com/file/d/1xhJRTJa49Qdj1R00PISWGVKHK2WeQKnJ/view?usp=share_link

mkdir moltaut_env
mv solv_rdkit_2020_env.tar.gz moltaut_env
tar -zxvf solv_rdkit_2020_env.tar.gz
source activate moltaut_env/bin/active
conda remove openbabel
conda install openbabel -c conda-forge
```

## Usage

```
python predict_tautomer.py --help

usage: predict_tautomer.py [-h] [--smi SMI] [--cutoff CUTOFF]
                           [--cutmol CUTMOL] [--num_confs NUM_CONFS] [--ph PH]
                           [--tph TPH] [--output OUTPUT]

calculate low-energy tautomer for small molecules

optional arguments:
  -h, --help            show this help message and exit
  --smi SMI             the molecular smiles
  --cutoff CUTOFF       the energy cutoff for low energy
  --cutmol CUTMOL       determine to frag the molecule
  --num_confs NUM_CONFS
                        the number of conformation for solvation energy
                        prediction
  --ph PH               the target pH for protonation states generation
  --tph TPH             pH tolerance for protonation states generation
  --output OUTPUT       the output SDF file name

```

Or you can access the web server of MolTaut--[URL](http://moltaut.xundrug.cn/)

