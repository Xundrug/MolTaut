#! /bin/bash

if [ ! -f solv_rdkit_2020_env.tar.gz ]; then
	echo "please download solv_rdkit_2020_env.tar.gz from google drive"
	exit
fi

mkdir solv_env
mv solv_rdkit_2020_env.tar.gz solv_env
cd solv_env
tar -zxvf solv_rdkit_2020_env.tar.gz
source activate bin/active

conda remove openbabel -y
conda install openbabel -c conda-forge -y


