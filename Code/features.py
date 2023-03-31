#! /usr/bin/env python3

import torch
import mol2
import numpy as np
import pickle
import os 
from random import sample

files = "../Data/final_data/"
structures={}

for i in os.listdir(files):
    structures[i] = ['cavityALL.mol2', 'protein.mol2', 'site.mol2', 'ligand.mol2']

random_files = sample(list(structures.keys()), 1000)

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

def get_angles(coords):
    """
    Calcula los ángulos phi y psi de una proteína a partir de sus coordenadas.
    El angulo se calcula en radianes. Utiliza la formula de la tangente inversa.

    Args:
        coords (list): Lista de listas que contienen las coordenadas de los
        átomos de la proteína.

    Returns:
        phi (list): Lista de los ángulos phi de la proteína.
        psi (list): Lista de los ángulos psi de la proteína.
    """
    coords=coords.tolist()
    phi = []
    psi = []
    for i in range(1, len(coords) - 1):
        x1, y1, z1 = coords[i - 1]
        x2, y2, z2 = coords[i]
        x3, y3, z3 = coords[i + 1]

        # Calcula el ángulo phi
        phi_i = np.arctan2(y2 - y1, x2 - x1) - np.arctan2(z2 - z1, np.sqrt((x2 - x1)**2 + (y2 - y1)**2))
        phi.append(phi_i)

        # Calcula el ángulo psi
        psi_i = np.arctan2(y2 - y1, x2 - x1) - np.arctan2(z3 - z2, np.sqrt((x3 - x2)**2 + (y3 - y2)**2))
        psi.append(psi_i)

    return phi, psi

def extract_features(data):
    cur = 1
    print(f"################################# {cur}/{len(data)}")
    for i in data:
        print(f"Now processing: {i}")
        cur_feat = {}
        print(f"Loading protein.")
        protein = f"{files}{i}/{structures[i][0]}"
        print(f"Loading cavity.")
        cavity = f"{files}{i}/{structures[i][-1]}"
        print(f"Loading site.")
        site = f"{files}{i}/{structures[i][1]}"
        print(f"Loading ligand.")
        ligand = f"{files}{i}/{structures[i][2]}"

        print("Creating Protein instance.")
        cur_prot = mol2.Protein(i, protein, cavity, site, ligand)

        # Calculate the distances
        print("Calculating distances.")
        # protein_distance = get_distances(cur_prot.get_protein())   
        cavity_distance = get_distances(cur_prot.get_cavity())   
        # site_distance = get_distances(cur_prot.get_site())   
        site_ligand_dist = get_distances2(cur_prot.get_ligand(), cur_prot.get_siteCB())

        # Calculate the angles
        print("Calculating angles.")
        protein_angles_phi, protein_angles_psi = get_angles(cur_prot.get_proteinCA())  
        site_angles_phi, site_angles_psi = get_angles(cur_prot.get_siteCA())  

        # Transform the matrix into tensors
        print("Transforming data to tensors.")
        cur_feat["protein_TENSOR"] = torch.from_numpy(cur_prot.get_protein())
        cur_feat["cavity_TENSOR"] = torch.from_numpy(cur_prot.get_cavity())
        cur_feat["site_TENSOR"] = torch.from_numpy(cur_prot.get_site())
        cur_feat["ligand_TENSOR"] = torch.from_numpy(cur_prot.get_ligand())
        # cur_feat["protein_DISTANCES"] = torch.from_numpy(protein_distance)
        cur_feat["cavity_DISTANCES"] = torch.from_numpy(cavity_distance)
        # cur_feat["site_DISTANCES"] = torch.from_numpy(site_distance)
        cur_feat["site_ligand_DISTANCES"] = torch.from_numpy(site_ligand_dist)
        cur_feat["protein_PHI"] = torch.tensor(protein_angles_phi)
        cur_feat["protein_PSI"] = torch.tensor(protein_angles_psi)
        cur_feat["site_PHI"] = torch.tensor(site_angles_phi)
        cur_feat["site_PSI"] = torch.tensor(site_angles_psi)

        # Save each feature in the dictionary
        print(f"Saving features: {i}")
        features[i] = cur_feat
        cur += 1 
        print(f"################################# {cur}/{len(data)}")
    return features

data = extract_features(random_files)

print("Saving into file.")
out_fd = open("features.p", "wb")

pickle.dump(data, out_fd)

out_fd.close()
print("Done.")