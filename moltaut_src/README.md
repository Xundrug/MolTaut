# Scripts associated with atomic feature calculating, molecular conformation optimizing, energies calculating, splitting molecules into fragments, and tautomer generation.

## Manifest

* `aevs.py`: This script is to calculate atomic environment vectors using openbabel (AEVs).
* `sasa.py`: It's a script to calculate atomic SASA by FreeSASA.
* `descriptor.py`: Used to combine all atomic features and calculate atom type information.
* `models.py`: The model of MolSolv, which is built by Pytorch and PyG.
* `gen_confs.py`: To generate molecular conformation by RDKit.
* `optimize_mol.py`: To optimize molecular conformation by ANI-2x model.
* `rank_tautomer.py`: Used to calculate internal energies and solvation energies to get the scores, then rank the tautomers by the scores.
* `smirks_tansform_all.txt`: The SMIRTKS string for tautomer generation using RDKit reaction module.
* `tautomer.py`: This script is used to generate tautomers for each molecule.
* `cut_mol.py`: This script is used to split the molecule into the different fragments.
* `combine_frag.py`: Used to combine different fragment molecules by the linker information to get the full molecule.
* `molgpka`: The package for micro-pka prediction and protonation generation, details are shown in https://github.com/Xundrug/MolGpKa.

