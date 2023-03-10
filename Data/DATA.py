#! /usr/bin/env python

import pandas as pd
import urllib.request
import sys
import os
import warnings

try:
    sys.argv[1]
except IndexError:
    raise Exception("Please introduce a valid .tsv file.")

warnings.filterwarnings("ignore")

print("Loading data...")
df = pd.read_csv(sys.argv[1], sep = "\t", on_bad_lines="skip")

pdb_codes = []

try:
    list1 = df["PDB ID(s) of Target Chain"].tolist()
except:
    raise Exception("Please introduce a valid file from BindingDB.")

def extract_seqs(input, output):
    try:
        for i in input:
            if i != "nan":
                for j in i.split(","):
                    if j not in output:
                        output.append(j)
    except AttributeError:
        pass

print("Extracting PDB files")
extract_seqs(list1,pdb_codes)

print(f"{len(pdb_codes)} files found.")

print("Starting download")
for code in set(pdb_codes):
    pdb_url = f"https://files.rcsb.org/download/{code}.pdb"
    pdb_path = f"PDB/{code}.pdb"
    if f"{code}.pdb" not in os.listdir("PDB/"):
        urllib.request.urlretrieve(pdb_url, pdb_path)
        print(f"PDB file {code} downloaded and saved to {pdb_path}")
    else:
        print(f"PDB file {code} already found in {pdb_path}")











