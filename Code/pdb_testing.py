#! /usr/bin/env python3

from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import sys
import pickle
import pandas as pd
import df_maker
import warnings
from scipy.signal import find_peaks

warnings.filterwarnings("ignore")

print("Reading file.")
try:
    protein = sys.argv[1]
except IndexError:
    raise IndexError("Please introduce a PDB file.")

try:
    with open(sys.argv[2], "rb") as file:
        USER_model = pickle.load(file)
except IndexError:
    print("No user model provided.")

print("Loading model.")
with open("model_6.pckl", "rb") as file:
    model_6 = pickle.load(file)
with open("model_7.pckl", "rb") as file:
    model_7 = pickle.load(file)
with open("model_8.pckl", "rb") as file:
    model_8 = pickle.load(file)

print("Creating dataframe.")
pdb_df = df_maker.NewDataFrame(pdb_file = protein, type="pdb")

print("Extracting sequnce.")
amino_keys = {
    0: 'ALA', 1: 'ARG', 2: 'ASN', 3: 'ASP', 4: 'CYS', 5: 'GLN', 
    6: 'GLU', 7: 'GLY', 8: 'HIS', 9: 'ILE', 10: 'LEU', 11: 'LYS', 
    12: 'MET', 13: 'PHE', 14: 'PRO', 15: 'SER', 16: 'THR', 
    17: 'TRP', 18: 'TYR', 19: 'VAL'
}
#
#SS_keys = {
#    0: '-', 1: 'B', 2: 'T', 3: 'S', 4: 'G', 5: 'E', 6: 'H', 7: 'I'
#}
#
#parser = PDBParser()
#structure = parser.get_structure("protein", protein)
#
## initialize a PPBuilder to extract the sequence
#ppb = PPBuilder()
#
## initialize an empty sequence string
#sequence = ""
#
## iterate over all chains and append their sequences to the overall sequence string
#for chain in structure.get_chains():
#    sequence += "".join([residue.get_resname() for residue in chain.get_residues() if residue.get_resname() not in ["HOH", "H2O"]])
#
#sequence = list(sequence)

print("Predicting file.")
try:
    prediction = USER_model.predict(pdb_df)
    print("Trying with user model.")
except NameError:
    print("Trying with model 6.")
    prediction = model_6.predict(pdb_df)
    if 1 not in prediction:
        print("Trying with model 7.")
        prediction = model_7.predict(pdb_df)
        if 1 not in prediction:
            print("Trying with model 8.")
            prediction = model_8.predict(pdb_df)

#predicted_sequence = ""
#for pred, seq in zip(prediction, sequence):
#    if pred == 1:
#        predicted_sequence += seq

#print(prediction)

name = protein[:-4]

print("")
atoms = []
all_atoms = []
with open(protein, "r") as file:
    for line in file:
        if line.startswith("ATOM") and line[13:16].rstrip() == "CA":
            atoms.append(line)
        if line.startswith("ATOM"):
            all_atoms.append(line)
        
print("Finding the binding residues.")
res_pos_list = []
for pred, ato in zip(prediction, atoms):
    if pred == 1:
        res_pos_list.append((ato[17:21], ato[23:26]))

# To make sure the file is empty
with open(f"{name}_prediction.pdb", "w") as file:
    file.write("")

print("Saving data to new file.")
with open(f"{name}_prediction.pdb", "a") as file:
    for res, pos in res_pos_list:
        for ato in all_atoms:
            if res == ato[17:21] and pos == ato[23:26]:
                file.writelines(ato)

print(f"Predicted Binding site to '{name}_prediction.pdb'.")

aa_all_pos = []
for pred, ato in zip(prediction, atoms):
        aa_all_pos.append((ato[17:21], ato[23:26]))

Coded_AA = pdb_df['AA'].values
Decoded_AA = pd.Series(Coded_AA).map(amino_keys).values
pdb_df['AA'] = Decoded_AA
for i, pos in enumerate(aa_all_pos):
    pdb_df.at[i, 'AA'] = pdb_df.at[i, 'AA'] + '-' + str(pos[1])

potential_aa = []
for i in range(len(prediction)):
    if prediction[i] == 1:
        potential_aa.append(pdb_df.loc[i, 'AA'])

peaks, _ = find_peaks(prediction, height=0.5)
sites = []
for i in range(len(peaks)):
    start = peaks[i]
    try:
        end = peaks[i+1]
    except IndexError:
        end = len(prediction)
    site = pdb_df.loc[start:end-1, 'AA'].values
    if any(aa in site for aa in potential_aa):
        sites.append(site)

sites_sorted = sorted(sites, key=len)

with open(f"{name}_binding_site_predictions.txt", "w") as f:
    f.write(f"Amino_b acids that can form part of a binding site: {', '.join(potential_aa)}\n\n")
    bs_count = 1
    for site in sites_sorted:
        aa_string = " ".join(site)
        f.write(f"BS predicted {bs_count}: {aa_string}\n\n")
        bs_count += 1
print("Saved predicted binding sites to file.")