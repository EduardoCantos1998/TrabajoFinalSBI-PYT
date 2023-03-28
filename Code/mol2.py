#! /usr/bin/env python3

import os

files = "../Data/final_data/"

structures = {}

class Protein:
    def __init__(self, name, protein, cavity, site):
        self.name = name
        self.protein = protein
        self.cavity = cavity
        self.site = site

    def get_protein(self):
        with open(self.protein, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write:
                    atoms = line.split()
                    struct.append((atoms[2],atoms[3],atoms[4]))
        return struct

    def get_cavity(self):
        with open(self.cavity, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write:
                    atoms = line.split()
                    struct.append((atoms[2],atoms[3],atoms[4]))
        return struct

    def get_site(self):
        with open(self.site, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write:
                    atoms = line.split()
                    struct.append((atoms[2],atoms[3],atoms[4]))
        return struct

    def __hash__(self):
        return hash(self.name)


for i in os.listdir(files):
    structures[i] = ['cavityALL.mol2', 'protein.mol2', 'site.mol2']

for i in structures:
    protein = f"{files}{i}/{structures[i][1]}"
    cavity = f"{files}{i}/{structures[i][0]}"
    site = f"{files}{i}/{structures[i][2]}"

    cur_prot = Protein(i, protein, cavity, site)
    

