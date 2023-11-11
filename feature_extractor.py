import numpy as np
import argparse
import yaml
import time
from tqdm import tqdm
import h5py
import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# Define a function to calculate Morgan fingerprints for a given SMILES string
def calculate_morgan_fingerprint(smiles_string, radius, nBits, useChirality):
    #print("I am process", os.getpid())
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is not None:
        morgan_fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useChirality=useChirality)
        morgan_fingerprint_array = np.zeros((1,))
        AllChem.DataStructs.ConvertToNumpyArray(morgan_fingerprint, morgan_fingerprint_array)
        return morgan_fingerprint_array
    else:
        return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="SMILES string")
    parser.add_argument("--save_path", type=str, required=True, help="path to store features")
    parser.add_argument("--experiment", type=str, required=True, help="experiment name")
    args = parser.parse_args()

    with open('./config.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    print("Loading file...")
    smiles_df = pd.read_csv(args.input_file, usecols=['SMILES'])
    smiles_strings = smiles_df['SMILES'].tolist()

    calculate_morgan = partial(calculate_morgan_fingerprint, radius=config['radius'], 
                               nBits=config['nBits'], useChirality=config['useChirality'])
    # Create a ProcessPoolExecutor to parallelize fingerprint calculation
    print("Calculate fingerprints...")
    t_start = time.time()
    morgan_fingerprints = []
    # for smile in tqdm(smiles_strings):
    #     morgan_fingerprints.append(calculate_morgan(smile))
    chunk_size = config['chunk_size']
    for i in range(0, len(smiles_strings), config['chunk_size']):
        print(f"Processing chunk {i} to {i+chunk_size}")
        with ProcessPoolExecutor(max_workers=config['num_processes']) as executor:
            tmp = list(tqdm(executor.map(calculate_morgan, smiles_strings[i: i+chunk_size]), total=len(smiles_strings[i: i+chunk_size])))
            morgan_fingerprints.extend(tmp)

    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")

    # morgan_fingerprints is a list of numpy arrays containing the fingerprints
    smile_valid = []
    fingerprint_valid = []
    smile_invalid = []
    for i, smile in enumerate(smiles_strings):
        if morgan_fingerprints[i] is None:
            smile_invalid.append(smile)
        else:
            smile_valid.append(smile)
            fingerprint_valid.append(morgan_fingerprints[i])
    
    print("Saving features...")
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    with h5py.File(os.path.join(args.save_path, f'{args.experiment}.h5'), "w") as hdf5_file:
        # Create datasets for SMILES and Morgan fingerprints
        hdf5_file.create_dataset("SMILES", data=np.array(smile_valid, dtype="S"), compression="gzip")
        hdf5_file.create_dataset("feature", data=fingerprint_valid, compression="gzip")
        hdf5_file.create_dataset("invalid_SMILES", data=np.array(smile_invalid, dtype="S"))
            


