#! /bin/bash

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1xhJRTJa49Qdj1R00PISWGVKHK2WeQKnJ' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1xhJRTJa49Qdj1R00PISWGVKHK2WeQKnJ" -O solv_rdkit_2020_env.tar.gz && rm -rf /tmp/cookies.txt

if [ ! -f solv_rdkit_2020_env.tar.gz ]; then
	echo "Auto download conda-pack env file failed. Please download solv_rdkit_2020_env.tar.gz from google drive manually."
	exit
fi

mkdir solv_env
mv solv_rdkit_2020_env.tar.gz solv_env
cd solv_env
tar -zxvf solv_rdkit_2020_env.tar.gz
source activate bin/active

conda remove openbabel -y
conda install openbabel -c conda-forge -y
pip install mols2grid

