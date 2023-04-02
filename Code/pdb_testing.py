#! /usr/bin/env python3

from Bio.PDB import PDBParser
import numpy as np
import sys
import pickle
import model as df_gen # The NewDataFrame generator function

print("Reading file.")
try:
    protein = sys.argv()[1]
except IndexError:
    raise "Please introduce a PDB file."

print("Loading model.")
with open("model.pckl", "rb") as file:
    model = pickle.load(file)

path = sys.argv()[1]

print("Creating dataframe.")
pdb_df = df_gen.NewDataFrame(pdb_file=path, type="pdb")

print(pdb_df)

print("Predicting file.")
prediction = model.predict(pdb_df)

pdb_df["PREDICTED_ATOM"] = prediction

print("Saving data to new file.")
with open("prediction.pdb", "w") as file:
    pass