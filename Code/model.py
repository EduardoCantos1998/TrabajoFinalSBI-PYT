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
    structures[i] = ['cavityALL.mol2', 'protein.mol2', 'site.mol2', 'ligand.mol2']

features = {}

def get_distances(X):
    distances = np.zeros((len(X),len(X)))
    for i, coord1 in enumerate(X):
        for j, coord2 in enumerate(X):
            distances[i,j] = np.linalg.norm(coord1-coord2)
    return distances

def get_distances2(X,Y):
    distances = np.zeros((len(X),len(Y)))
    for i, coord1 in enumerate(X):
        for j, coord2 in enumerate(Y):
            distances[i,j] = np.linalg.norm(coord1-coord2)
    return distances

for i in structures:
    cur_feat = {}
    protein = f"{files}{i}/{structures[i][0]}"
    cavity = f"{files}{i}/{structures[i][-1]}"
    site = f"{files}{i}/{structures[i][1]}"
    ligand = f"{files}{i}/{structures[i][2]}"

    cur_prot = mol2.Protein(i, protein, cavity, site, ligand)

    # Calculate the distances
    protein_distance = get_distances(cur_prot.get_protein())   
    cavity_distance = get_distances(cur_prot.get_cavity())   
    site_distance = get_distances(cur_prot.get_site())   
    site_ligand_dist = get_distances2(cur_prot.get_ligand(), cur_prot.get_siteCB())

    # Transform the matrix into tensors
    cur_feat["protein_TENSOR"] = torch.from_numpy(cur_prot.get_protein())
    cur_feat["cavity_TENSOR"] = torch.from_numpy(cur_prot.get_cavity())
    cur_feat["site_TENSOR"] = torch.from_numpy(cur_prot.get_site())
    cur_feat["ligand_TENSOR"] = torch.from_numpy(cur_prot.get_ligand())
    cur_feat["protein_DISTANCES"] = torch.from_numpy(protein_distance)
    cur_feat["cavity_DISTANCES"] = torch.from_numpy(cavity_distance)
    cur_feat["site_DISTANCES"] = torch.from_numpy(site_distance)
    cur_feat["ligand_DISTANCES"] = torch.from_numpy(ligand_distance)
    cur_feat["site_ligand_DISTANCES"] = torch.from_numpy(site_ligand_dist)
    # Save each feature in the dictionary
    features[i] = cur_feat

proteins = features.keys()
print(proteins)
