#! /usr/bin/env python3

"""
The way this code works is the following:

        Extract the coordinates of the protein/ligand/cavity/site from a mol2 file. We've
        also included the alpha carbons for calculating the distances within the proteins,
        and the beta carbons to calculate the distance between the known binding 
        site/cavity and the ligand.

"""

import numpy as np

class Protein():
    def __init__(self, name, protein, cavity, site, ligand):
        self.name = name
        self.protein = protein
        self.cavity = cavity
        self.site = site
        self.ligand = ligand

    ### Define the protein coordinates
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

    def get_proteinCA(self):
        with open(self.protein, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write and line.split()[1] == "CA":
                    atoms = line.split()
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_protein = np.matrix(struct)
        return np_protein

    def get_proteinCB(self):
        with open(self.protein, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write and line.split()[1] == "CB":
                    atoms = line.split()
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4]), atoms[7]])
        np_protein = np.matrix(struct)
        return np_protein

    ### Define the cavity coordinates
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

    ### Define the site coordinates
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

    def get_siteCA(self):
        with open(self.site, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write and line.split()[1] == "CA":
                    atoms = line.split()
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4])])
        np_site = np.matrix(struct)
        return np_site

    def get_siteCB(self):
        with open(self.site, "r") as file:
            struct = []
            write = False
            for line in file:
                if line.startswith("@<TRIPOS>ATOM"):
                    write = True
                    continue
                elif line.startswith("@"):
                    write = False
                if write and line.split()[1] == "CB":
                    atoms = line.split()
                    struct.append([float(atoms[2]),float(atoms[3]),float(atoms[4]), atoms[7]])
        np_site = np.matrix(struct)
        return np_site

    ### Define the ligand coordinates
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
