#! /usr/bin/env python3

import torch
import mol2
import numpy as np
import os 
from sklearn.model_selection import train_test_split
import Bio

files = "../Data/final_data/"
structures={}

for i in os.listdir(files):
    structures[i] = ['cavity6.mol2', 'protein.mol2', 'site.mol2', 'ligand.mol2']

features = {}

for i in structures:
    cur_feat = {}
    protein = f"{files}{i}/{structures[i][0]}"
    cavity = f"{files}{i}/{structures[i][-1]}"
    site = f"{files}{i}/{structures[i][1]}"
    ligand = f"{files}{i}/{structures[i][2]}"

    cur_prot = mol2.Protein(i, protein, cavity, site, ligand)
    # Transform the matrix into tensors
    cur_feat["protein_TENSOR"] = torch.from_numpy(cur_prot.get_protein())
    cur_feat["cavity_TENSOR"] = torch.from_numpy(cur_prot.get_cavity())
    cur_feat["site_TENSOR"] = torch.from_numpy(cur_prot.get_site())
    cur_feat["ligand_TENSOR"] = torch.from_numpy(cur_prot.get_ligand())
    # Save each feature in the dictionary
    features[i] = cur_feat
