#! /usr/bin/env python3

import torch
import mol2
import numpy as np
import os 

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

def get_angles(coords):
    """
    Calcula los ángulos phi y psi de una proteína a partir de sus coordenadas.
    El angulo se calcula en radianes.

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

    # Calculate the angles
    protein_angles_phi, protein_angles_psi = get_angles(cur_prot.get_proteinCA())  
    site_angles_phi, site_angles_psi = get_angles(cur_prot.get_siteCA())  

    # Transform the matrix into tensors
    cur_feat["protein_TENSOR"] = torch.from_numpy(cur_prot.get_protein())
    cur_feat["cavity_TENSOR"] = torch.from_numpy(cur_prot.get_cavity())
    cur_feat["site_TENSOR"] = torch.from_numpy(cur_prot.get_site())
    cur_feat["ligand_TENSOR"] = torch.from_numpy(cur_prot.get_ligand())
    cur_feat["protein_DISTANCES"] = torch.from_numpy(protein_distance)
    cur_feat["cavity_DISTANCES"] = torch.from_numpy(cavity_distance)
    cur_feat["site_DISTANCES"] = torch.from_numpy(site_distance)
    cur_feat["site_ligand_DISTANCES"] = torch.from_numpy(site_ligand_dist)
    cur_feat["protein_PHI"] = torch.from_numpy(protein_angles_phi)
    cur_feat["protein_PSI"] = torch.from_numpy(protein_angles_psi)
    cur_feat["site_PHI"] = torch.from_numpy(site_angles_phi)
    cur_feat["site_PSI"] = torch.from_numpy(site_angles_psi)

    # Save each feature in the dictionary
    features[i] = cur_feat
