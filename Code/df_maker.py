import mol2
import freesasa
import pandas as pd
import numpy as np
import math
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.SASA import *
from Bio.PDB import *
import Bio.PDB
from Bio.PDB import PDBParser, Selection
from Bio.PDB.DSSP import DSSP
import tempfile
import os

def NewDataFrame(protein = None, pdb_file = None, type = "mol2"):
    """
    This function creates a new pandas dataframe each time it's called. It takes
    three parameters:

        - protein: this is the Protein object that we've generated using the mol2 file.
        - pdb_file: this is the path to the pdb file we're interested in.
        - type: this is the type of file. It takes two options:
            -'mol2': this comes by default
            -'pdb'

    The function will return the format in a pandas data frame. It will also
    calculate the protein's average distance between alpha carbons. Only the
    alpha carbons will be saved for the sake of simplicity. It will also take
    the distance between the protein and the cavity, the protein and the ligand.
    """

    def get_avg_distance(X):
        avg_distances = []
        for i, coord1 in enumerate(X):
            distances = []
            for j, coord2 in enumerate(X):
                distances.append(np.linalg.norm(coord1-coord2))
            avg_distances.append(sum(distances)/len(distances))
        return avg_distances

    def get_avg_distance2(X, Y):
        avg_distances = []
        for i, coord1 in enumerate(X):
            distances = []
            for j, coord2 in enumerate(Y):
                distances.append(np.linalg.norm(coord1-coord2))
            avg_distances.append(sum(distances)/len(distances))
        return avg_distances

    def get_angles(coords):
        """
        Calculates the phi and psi angles of a protein from its coordinates. The
        angle is calculated in radians. It uses the inverse tangent formula.

        Args:
            coords (list): List of lists containing the coordinates of the
            protein atoms.

        Returns:
            phi (list): List of the phi angles of the protein.
            psi (list): List of the psi angles of the protein.
        """
        try:
            coords=coords.tolist()
        except AttributeError:
            pass
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
    
    # Define the pdb protein
    def b_fact_calculator(structure):
        # Calculate the average B-factor for each residue
        b_fact = []
        for model in structure:
            for chain in model:
                ppb = PPBuilder()
                for pp in ppb.build_peptides(chain):
                    #residues = pp.get_sequence()
                    for residue in pp:
                        b_factor_sum = sum(atom.bfactor for atom in residue)
                        b_factor_avg = b_factor_sum / len(residue)
                        b_fact.append(b_factor_avg)
                        #print(f"{residue.get_full_id()[3][1]} {residue.get_resname()} {b_factor_avg}")
        return b_fact
    
    # Define SASA
    def getSASA(alpha_carbons):
        """
        Calculates the SASA of each amino acid based on the coordinates of its alpha carbon.

        Args:
        alpha_carbons (list): A list of the coordinates of the alpha carbons.

        Returns:
        sasa_values (list): A list of the SASA values for each amino acid.
        """

        # Convert the alpha_carbons list to a numpy array
        alpha_carbons = np.array(alpha_carbons)

        # Define the radius of a carbon atom (in Angstroms)
        carbon_radius = 1.7

        # Calculate the pairwise distance matrix between all alpha carbons
        pairwise_distances = np.linalg.norm(alpha_carbons[:, np.newaxis, :] - alpha_carbons, axis=-1)

        # Calculate the SASA for each amino acid
        sasa_values = []
        for i, d_i in enumerate(pairwise_distances):
            # Calculate the surface area of a sphere with radius equal to the distance
            # between the alpha carbon of the current amino acid and all other alpha carbons
            surface_area = 4 * np.pi * np.power(d_i, 2)

            # Subtract the surface area of the spheres corresponding to neighboring amino acids
            neighboring_surface_area = 0
            for j, d_j in enumerate(pairwise_distances[i+1:], i+1):
                if np.less(d_j, (d_i + carbon_radius)).any():
                    neighboring_surface_area += 4 * np.pi * np.power(d_j, 2)

            sasa = surface_area - neighboring_surface_area
            sasa_values.append(sasa)

        return sasa_values

    def calculate_secondary_structure(structure):
        ss_list = []
        aa_list = ['G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W', 'N', 'Q', 'S', 'T', 'Y', 'C', 'D', 'E', 'K', 'R', 'H']
        # Create a DSSP object
        for i in structure:
            model = i
            dssp = DSSP(model, pdb_file)
            # Iterate over each residue and print its secondary structure
            for residue in dssp:
                if residue[1] in aa_list:
                    #res_id = residue[0]
                    ss = residue[2]
                    ss_list.append(ss)
        return ss_list
        
    # Charges
    def get_residue_pI(amino_list):
        # Crear un objeto PPBuilder para obtener la secuencia de aminoácidos
        #ppb = PPBuilder()
        #seq = ppb.build_peptides(structure)[0].get_sequence()

        # Crear un diccionario de cargas por aminoácido
        residue_pI = {
            'ARG': 10.76,
            'HIS': 7.59,
            'LYS': 9.74,
            'ASP': 2.77,
            'GLU': 3.22,
            'SER': 5.68,
            'THR': 5.60,
            'ASN': 5.41,
            'GLN': 5.65,
            'CYS': 5.07,
            'SEC': 5.07,  # assuming selenocysteine has the same pI as cysteine
            'GLY': 5.97,
            'PRO': 6.30,
            'ALA': 6.00,
            'VAL': 5.96,
            'ILE': 6.02,
            'LEU': 5.98,
            'MET': 5.74,
            'PHE': 5.48,
            'TYR': 5.66,
            'TRP': 5.89,
        }

        # Calcular la carga neta de cada residuo
        pI_list = []
            
        for res in amino_list:
            pI_list.append(residue_pI[res])
        
        return pI_list
    
    # Hydrophobicity
    def get_hydrophobicity(structure):
        # Crear un objeto PPBuilder para obtener la secuencia de aminoácidos
        ppb = PPBuilder()
        seq = ppb.build_peptides(structure)
        sequences = [peptide.get_sequence() for peptide in seq]

        # Escala de hidrofobicidad Kyte-Doolittle
        kd_hydrophobicity = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
                        'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
                        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
                        'W': -0.9, 'Y': -1.3}

        # Calcular la hidrofobicidad de cada residuo
        hydrophobicity_list = []
        for chain in sequences:
            for residue in chain:
                aa = ProteinAnalysis(str(residue)).get_amino_acids_percent()[str(residue)]
                hydrophobicity = aa * kd_hydrophobicity[str(residue)]
                hydrophobicity_list.append(hydrophobicity)

        return hydrophobicity_list

    # Entropy
    def get_entropies(structure):
        """
        Calculates the entropy of each residue in a protein.

        Args:
            pdb_id (str): A string identifier for the protein.
            pdb_file (str): The path to the PDB file.

        Returns:
            entropies (list): A list of the entropy values for each residue in the protein.
        """
        # Create an empty list to store the entropies
        entropies = []
        aa_list = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'CYS', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS']

        # Iterate over each residue in the protein
        for residue in structure.get_residues():
            if residue.get_resname().upper() in aa_list:
                # Get the one-letter code for the residue
                residue_letter = residue.get_resname().upper()
                # Calculate the probability of the residue occurring in the protein
                residue_count = 0
                total_residues = 0
                for other_residue in structure.get_residues():
                    if is_aa(other_residue):
                        total_residues += 1
                        if other_residue.get_resname().upper() == residue_letter:
                            residue_count += 1
                residue_probability = residue_count / total_residues
                # Calculate the entropy of the residue
                residue_entropy = -1 * residue_probability * math.log2(residue_probability)
                entropies.append(residue_entropy)

        return entropies

    # Sequence
    def get_sequence(pdb_file):
        """
        Parses a PDB file and returns a list of the residue letters.

        Args:
        pdb_file (str): The path to the PDB file.

        Returns:
        sequence (list): A list of the one-letter codes for each residue in the protein.
        """
        structure = Bio.PDB.PDBParser().get_structure('pdb', pdb_file)
        ppb = Bio.PDB.PPBuilder()
        sequence = []
        for pp in ppb.build_peptides(structure):
            sequence += [residue.resname for residue in pp]
        return sequence

    # Define the dataframe
    if type == "mol2":
        #new_df = pd.DataFrame(protein.get_proteinCA())
        print("Finding Site Atoms")
        site_atoms = protein.get_siteCA().tolist()
        #protein_atoms = protein.get_proteinCA().tolist()
    #elif type == "pdb":
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    # create an empty list to store the alpha carbons
    alpha_carbons = []

    print("Finding AlphaCarbons")
    # iterate over all the atoms in the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # check if the residue is an amino acid
                if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    # get the alpha carbon atom
                    alpha_carbon = residue['CA']
            
                    # add the alpha carbon atom to the list
                    alpha_carbons.append(alpha_carbon)

    print("Create Dataframe")
    # create a DataFrame with the coordinates of all the alpha carbons
    new_df = pd.DataFrame([list(atom.get_coord()) for atom in alpha_carbons], columns=['X_COORD', 'Y_COORD', 'Z_COORD'])
    alpha_carbons = [list(atom.get_coord()) for atom in alpha_carbons]

    if type == "mol2":
        print("Calculating angles.")
        proteinCA_angles_PHI, proteinCA_angles_PSI = get_angles(alpha_carbons)
        print("Saving binding site.")
        binding = []
        for atom in protein.get_proteinCA().tolist():
            if atom in site_atoms:
                binding.append(1)
            else:
                binding.append(0)

    elif type == "pdb":
        print("Calculating angles.")
        proteinCA_angles_PHI, proteinCA_angles_PSI = get_angles(alpha_carbons)
    
    proteinCA_angles_PHI.append(0)
    proteinCA_angles_PHI.insert(0,0)
    proteinCA_angles_PSI.append(0)
    proteinCA_angles_PSI.insert(0,0)
    #if type == "mol2":
        #protein_distances = get_avg_distance(np.matrix(alpha_carbons))
        #protein_ligand = get_avg_distance2(protein.get_proteinCA(), protein.get_ligand())
        #protein_cavity = get_avg_distance2(protein.get_proteinCA(), protein.get_cavity())
    #elif type == "pdb":
    print("Finding protein's inner distance.")
    protein_distances = get_avg_distance(np.matrix(alpha_carbons))
    # Add the values to the dataframe.
    #if type == "mol2":
        #new_df.columns = ["X_COORD", "Y_COORD", "Z_COORD"]
    print("Saving features in DataFrame.")
    new_df["PROTEIN_AVG_LENGTH"] = protein_distances
    if type == "mol2":
        #new_df["PROTEIN_LIGAND_LENGTH"] = protein_ligand
        #new_df["PROTEIN_CAVITY_LENGTH"] = protein_cavity
        new_df["BINDING_ATOM"] = binding 
    #elif type == "pdb":
        #new_df["PROTEIN_LIGAND_LENGTH"] = [12] * len(protein_distances)
        #new_df["PROTEIN_CAVITY_LENGTH"] = [12] * len(protein_distances)
    new_df["PROTEIN_PSI"] = proteinCA_angles_PSI
    new_df["PROTEIN_PHI"] = proteinCA_angles_PHI
    new_df["AA"] = get_sequence(pdb_file)
    new_df["pI"] = get_residue_pI(get_sequence(pdb_file))
    #if type == "mol2":
    #    new_df["SASA"] = getSASA(alpha_carbons)
    #if type == "pdb": 
    SASA = getSASA(alpha_carbons)
    mean_SASA = []
    print("Calculating average SASA.")
    mean_SASA = []

    for i in SASA:
        mean_i = sum(i) / len(i)
        mean_SASA.append(mean_i)
    new_df["SASA"] = mean_SASA
    new_df["SECONDARY_STRUCTURE"] = calculate_secondary_structure(structure)
    new_df["B-FACTOR"] = b_fact_calculator(structure)
    new_df["HIDROPHOBICITY"] = get_hydrophobicity(structure)
    new_df["ENTROPY"] = get_entropies(structure)

    amino_keys = {
    'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLN': 5, 'GLU': 6, 'GLY': 7, 
    'HIS': 8, 'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14, 
    'SER': 15, 'THR': 16, 'TRP': 17, 'TYR': 18, 'VAL': 19
    }
    SS_dict = {'-': 0, 'B': 1, 'T': 2, 'S': 3, 'G': 4, 'E': 5, 'H': 6, 'P': 7, 'I': 8}

    new_df.AA = new_df.AA.map(amino_keys)
    new_df.SECONDARY_STRUCTURE = new_df.SECONDARY_STRUCTURE.map(SS_dict)
    print("Finished generating DataFrame.")
    return new_df