#! /usr/bin/env python3

import numpy as np

class Protein():
    def __init__(self, name, protein, cavity, site, ligand):
        self.name = name
        self.protein = protein
        self.cavity = cavity
        self.site = site
        self.ligand = ligand

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
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_protein = np.matrix(struct)
        return np_protein

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
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_cavity = np.matrix(struct)
        return np_cavity

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
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_site = np.matrix(struct)
        return np_site

    def get_ligand(self):
        with open(self.ligand, "r") as file:
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
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_ligand = np.matrix(struct)
        return np_ligand

    def __hash__(self):
        return hash(self.name)
