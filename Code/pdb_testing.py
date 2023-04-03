#! /usr/bin/env python3

from Bio.PDB import PDBParser
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

path = sys.argv[1]

print("Creating dataframe.")
pdb_df = df_maker.NewDataFrame(pdb_file = protein, type="pdb")

print(pdb_df)

print("Predicting file.")
prediction = model.predict(pdb_df)
print(prediction)

pdb_df["PREDICTED_ATOM"] = prediction

print("Saving data to new file.")
with open("prediction.pdb", "w") as file:
    pass

print("Predicted Binding site to 'prediction.pdb'.")