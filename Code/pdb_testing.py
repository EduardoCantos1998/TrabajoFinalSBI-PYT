#! /usr/bin/env python3

from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import sys
import pickle
import pandas as pd
import df_maker

print("Reading file.")
try:
    protein = sys.argv[1]
except IndexError:
    raise "Please introduce a PDB file."

print("Loading model.")
with open("model.pckl", "rb") as file:
    model = pickle.load(file)

print("Creating dataframe.")
pdb_df = df_maker.NewDataFrame(pdb_file = protein, type="pdb")

print("Extracting sequnce.")

amino_keys = {
    0: 'ALA', 1: 'ARG', 2: 'ASN', 3: 'ASP', 4: 'CYS', 5: 'GLN', 
    6: 'GLU', 7: 'GLY', 8: 'HIS', 9: 'ILE', 10: 'LEU', 11: 'LYS', 
    12: 'MET', 13: 'PHE', 14: 'PRO', 15: 'SER', 16: 'THR', 
    17: 'TRP', 18: 'TYR', 19: 'VAL'
}

SS_keys = {
    0: '-', 1: 'B', 2: 'T', 3: 'S', 4: 'G', 5: 'E', 6: 'H', 7: 'I'
}

parser = PDBParser()
structure = parser.get_structure("protein", protein)

# initialize a PPBuilder to extract the sequence
ppb = PPBuilder()

# initialize an empty sequence string
sequence = ""

# iterate over all chains and append their sequences to the overall sequence string
for chain in structure.get_chains():
    sequence += "".join([residue.get_resname() for residue in chain.get_residues() if residue.get_resname() not in ["HOH", "H2O"]])

sequence = list(sequence)

print("Predicting file.")
prediction = model.predict(pdb_df)
prediction = list(prediction)
predicted_sequence = ""
for pred, seq in zip(prediction, sequence):
    if pred == 1:
        predicted_sequence += seq

pdb_df["PREDICTED_ATOM"] = prediction

try:
    name = sys.argv[2]
except IndexError:
    name = protein[-8:-4]
print("Saving data to new file.")
with open(f"{name}_prediction.pdb", "w") as file:
    pass

print(f"Predicted Binding site to '{name}_prediction.pdb'.")