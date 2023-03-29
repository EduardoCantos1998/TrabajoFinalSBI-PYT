#! /usr/bin/env python3

from Bio.PDB import PDBParser
import numpy as np
import torch
import sys

try:
    protein = sys.argv()[1]
except IndexError:
    raise "Please introduce a PDB file."

parser = PDBParser()
structure = parser.get_structure("proteina", protein)
coords = []

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                coords.append(atom.get_coord())

coords = np.matrix(coords)
coords_TENSOR = torch.from_numpy(coords)
print(coords_TENSOR)