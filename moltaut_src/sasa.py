import freesasa
from openbabel.openbabel import OBMolAtomIter, OBAtomAtomIter
from openbabel import pybel
from collections import namedtuple
import numpy as np
freesasa.setVerbosity(1)


element_dict = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S", 9: "F", 17: "Cl"}

def get_atom_coords(at):
    acoords = np.array([at.GetX(), at.GetY(), at.GetZ()])
    return acoords


class DerivedClassifier(freesasa.Classifier):
    purePython = True

    def radius(self, residueName, atomName):
        radius = {"N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1, "F": 1.5, "Cl": 1.8,
                  "Br": 2.0, "I": 2.2, "C": 1.9, "H": 0.0}
        return radius[atomName]

def create_freesasa_structure(mol):
    new_s = freesasa.Structure()

    for atom in OBMolAtomIter(mol):
        x, y, z = get_atom_coords(atom)
        
        #resnr, restype, atom_name, atom_symbol = get_residue_info(atom)
        resnr = 0
        restype = 'UNL'
        #atom_name = atom.GetType()
        atom_name = element_dict[ atom.GetAtomicNum() ]
        atom_symbol = element_dict[ atom.GetAtomicNum() ]
        new_s.addAtom(atom_symbol, restype, str(resnr), "A", x, y, z)
    classifier = DerivedClassifier()
    new_s.setRadiiWithClassifier(classifier)
    return new_s

def calc_atoms_sasa(mol):
    data = namedtuple('sasa', 'value symbol idx')
    datas = []

    mol_freesasa = create_freesasa_structure(mol)
    sasa = freesasa.calc(mol_freesasa)

    for atom in OBMolAtomIter(mol):
        atom_idx = atom.GetIdx()-1
        sa = sasa.atomArea(atom_idx)
        datas.append(data(value=sa, symbol=element_dict[ atom.GetAtomicNum() ], idx=atom_idx+1))
    return datas

def get_sasa(mol):
    datas = calc_atoms_sasa(mol)
    sasa_dict = {}
    for data in datas:
        idx = data.idx
        sasa_dict[idx] = data.value
    return sasa_dict

if __name__=="__main__":
    filename = "t.sdf"
    pmol = next(pybel.readfile("sdf", filename))
    mol = pmol.OBMol
    datas = calc_atoms_sasa(mol)
    sasa_info = get_sasa(mol)
    print(sasa_info)
