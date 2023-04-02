import mol2
import pandas as pd
import numpy as np
import os
from random import sample
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.SASA import *
from Bio.PDB import *
import pickle
from Bio.PDB.DSSP import DSSP
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

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
        coords=coords.tolist()
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
    def b_fact_calculator(protein_name, file):
        # Load the PDB file
        parser = PDBParser()
        structure = parser.get_structure(protein_name, file)
        # Calculate the average B-factor for each residue
        b_fact = []
        for model in structure:
            for chain in model:
                ppb = PPBuilder()
                for pp in ppb.build_peptides(chain):
                    residues = pp.get_sequence()
                    for residue in pp:
                        b_factor_sum = sum(atom.bfactor for atom in residue)
                        b_factor_avg = b_factor_sum / len(residue)
                        b_fact.append(b_factor_avg)
                        #print(f"{residue.get_full_id()[3][1]} {residue.get_resname()} {b_factor_avg}")
        return b_fact
    
    ## Define SASA
    #def get_SASA(protein_name, file):
    #    # Cargar la estructura proteica desde un archivo PDB
    #    parser = PDBParser()
    #    structure = parser.get_structure(protein_name, file)
#
    #    # Crear un objeto PPBuilder para obtener la secuencia de aminoácidos
    #    ppb = PPBuilder()
    #    seq = ppb.build_peptides(structure)[0].get_sequence()
#
    #    # Crear un objeto NeighborSearch para buscar residuos cercanos
    #    ns = NeighborSearch(list(structure.get_atoms()))
#
    #    sasa_list = []
    #    for residue in structure.get_residues():
    #        if is_aa(residue):
    #            if 'CA' not in residue:
    #                sasa_list.append(0.0)
    #            else:
    #                center = residue['CA'].get_coord()
    #                neighbors = ns.search(center, 1.4)
    #                if len(neighbors) < 3:
    #                    sasa_list.append(0.0)
    #                else:
    #                    A = neighbors[0]['CA'].get_coord() - center
    #                    B = neighbors[1]['CA'].get_coord() - center
    #                    C = neighbors[2]['CA'].get_coord() - center
    #                    area = 0.5 * np.linalg.norm(np.cross(B - A, C - A))
    #                    sasa_list.append(area)
    #    return sasa_list


    #def calculate_secondary_structure(pdb_file):
    #            # Load the PDB file
    #    parser = PDBParser()
    #    structure = parser.get_structure("protein", pdb_file)

    #    # Extract the secondary structure information
    #    ss_list = []
    #    for model in structure:
    #        for chain in model:
    #            for residue in chain:
    #                #res_id = residue.get_id()[1]
    #                #res_name = residue.get_resname()
    #                if "CA" not in residue:
    #                    continue
    #                sec_struc = residue.get_full_id()[3][0]
    #                ss_list.append(sec_struc)
    #    return ss_list

    #ss = calculate_secondary_structure(pdb_file)

    # Charges
    def get_residue_charges(protein_name, file):
        # Cargar la estructura proteica desde un archivo PDB
        parser = PDBParser()
        structure = parser.get_structure(protein_name, file)

        # Crear un objeto PPBuilder para obtener la secuencia de aminoácidos
        ppb = PPBuilder()
        seq = ppb.build_peptides(structure)[0].get_sequence()

        # Crear un diccionario de cargas por aminoácido
        residue_charges = {
            'ARG': 1,
            'HIS': 0.5,
            'LYS': 1,
            'ASP': -1,
            'GLU': -1,
            'SER': 0,
            'THR': 0,
            'ASN': 0,
            'GLN': 0,
            'CYS': 0,
            'SEC': 0,
            'GLY': 0,
            'PRO': 0,
            'ALA': 0,
            'VAL': 0,
            'ILE': 0,
            'LEU': 0,
            'MET': 0,
            'PHE': 0,
            'TYR': -0.5,
            'TRP': -0.5,
        }

        # Calcular la carga neta de cada residuo
        charges_list = []
            
        for residue in structure.get_residues():
            if is_aa(residue):
                res_name = residue.get_resname()
                res_charge = residue_charges.get(res_name, 0)
                charges_list.append(res_charge)
        return charges_list
    

    def get_hydrophobicity(protein_name, file):
        # Cargar la estructura proteica desde un archivo PDB
        parser = PDBParser()
        structure = parser.get_structure(protein_name, file)

        # Crear un objeto PPBuilder para obtener la secuencia de aminoácidos
        ppb = PPBuilder()
        seq = ppb.build_peptides(structure)[0].get_sequence()

        # Escala de hidrofobicidad Kyte-Doolittle
        kd_hydrophobicity = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
                        'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
                        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
                        'W': -0.9, 'Y': -1.3}

        # Calcular la hidrofobicidad de cada residuo
        hydrophobicity_list = []
        for residue in seq:
                aa = ProteinAnalysis(str(residue)).get_amino_acids_percent()[str(residue)]
                hydrophobicity = aa * kd_hydrophobicity[str(residue)]
                hydrophobicity_list.append(hydrophobicity)

        return hydrophobicity_list

    # Define the dataframe
    if type == "mol2":
        new_df = pd.DataFrame(protein.get_proteinCA())
        site_atoms = protein.get_siteCA().tolist()
        protein_atoms = protein.get_proteinCA().tolist()
    elif type == "pdb":
        parser = PDBParser()
        structure = parser.get_structure("proteina", protein)
        
        # create an empty list to store the alpha carbons
        alpha_carbons = []

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
        new_df = pd.DataFrame(np.matrix(alpha_carbon))
    else: 
        raise("Please introduce a valid data type. Options are 'pdb' or 'mol2'.")
    if type == "mol2":
        proteinCA_angles_PHI, proteinCA_angles_PSI = get_angles(protein.get_proteinCA())
        binding = []
        for atom in protein_atoms:
            if atom in site_atoms:
                binding.append(1)
            else:
                binding.append(0)
    elif type == "pdb":
        proteinCA_angles_PHI, proteinCA_angles_PSI = get_angles(alpha_carbons)

    proteinCA_angles_PHI.append(0)
    proteinCA_angles_PHI.insert(0,0)
    proteinCA_angles_PSI.append(0)
    proteinCA_angles_PSI.insert(0,0)
    if type == "mol2":
        protein_distances = get_avg_distance(protein.get_proteinCA())
        protein_ligand = get_avg_distance2(protein.get_proteinCA(), protein.get_ligand())
        protein_cavity = get_avg_distance2(protein.get_proteinCA(), protein.get_cavity())
    elif type == "pdb":
        protein_distances = get_avg_distance(alpha_carbons)
    # Add the values to the dataframe.
    new_df.columns = ["X_COORD", "Y_COORD", "Z_COORD"]
    new_df["PROTEIN_AVG_LENGTH"] = protein_distances
    if type == "mol2":
        new_df["PROTEIN_LIGAND_LENGTH"] = protein_ligand
        new_df["PROTEIN_CAVITY_LENGTH"] = protein_cavity
        new_df["BINDING_ATOM"] = binding 
    new_df["PROTEIN_PSI"] = proteinCA_angles_PSI
    new_df["PROTEIN_PHI"] = proteinCA_angles_PHI
    new_df["CHARGES"] = get_residue_charges("protein", pdb_file)
    #new_df["SASA"] = get_SASA(protein.name, pdb_file)
    #new_df["SECONDARY_STRUCTURE"] = ss
    new_df["B-FACTOR"] = b_fact_calculator("protein", pdb_file)
    new_df["HIDROPHOBICITY"] = get_hydrophobicity("protein", pdb_file)

    return new_df

if __name__ == "__main__":
    # Generate the files
    files = "../Data/final_data/"
    structures={}

    for i in os.listdir(files):
        structures[i] = ['cavityALL.mol2', 'protein.mol2', 'site.mol2', 'ligand.mol2']

    # Generate the random files sample.
    random_files = sample(list(structures.keys()), 5020)

    print("Generating DataFrame dictionary.")
    df_dict = {}
    time = 0
    for i in random_files:
        print(f"#################### {time}/{len(random_files)}")
        time += 1
        try:
            path = f"../Data/final_data/{i}/"
            my_protein = mol2.Protein(name=f"i", protein=f"{path}protein.mol2", cavity=f"{path}cavity6.mol2", ligand=f"{path}ligand.mol2", site=f"{path}site.mol2")
            df_dict[i] = NewDataFrame(my_protein, f"{path}protein.pdb")
        except ValueError:
            pass
        except FileNotFoundError:
            pass

    print("Saving dictionary to pickle.")
    with open("dictionary.pckl", "wb") as file:
        pickle.dump(df_dict, file)
    print("File saved to 'dictionary.pckl'.")

    # Split the data into training, validation and test sets
    model = RandomForestClassifier(max_depth=2,random_state=0)

    print("Spliting the keys in train/test/validation sets.")
    train_keys, test_keys = train_test_split(list(df_dict.keys()), test_size=0.2, random_state=0)
    train_keys, val_keys = train_test_split(train_keys, test_size=0.2, random_state=0)

    # Fit the model on the training data and evaluate on the validation set
    best_accuracy = 0
    print("Starting the training.")
    for key in train_keys:
        print(f"Current protein: {key}")
        matrix = df_dict[key]
        X_train = matrix.drop("BINDING_ATOM", axis=1)
        y_train = matrix.BINDING_ATOM

        for val_key in val_keys:
            val_matrix = df_dict[val_key]
            X_val = val_matrix.drop("BINDING_ATOM", axis=1)
            y_val = val_matrix.BINDING_ATOM

            # Train the model on the training set and evaluate on the validation set
            print("Fitting data.")
            model.fit(X_train, y_train)
            y_val_pred = model.predict(X_val)
            val_accuracy = accuracy_score(y_val_pred, y_val)

            # Keep track of the best hyperparameters
            if val_accuracy > best_accuracy:
                print("Saving accuracy.")
                best_accuracy = val_accuracy
                best_params = model.get_params()

    # Train the final model on the combined training and validation sets
    train_val_keys = train_keys + val_keys
    X_train_val = pd.concat([df_dict[key].drop("BINDING_ATOM", axis=1) for key in train_val_keys])
    y_train_val = pd.concat([df_dict[key].BINDING_ATOM for key in train_val_keys])
    model.set_params(**best_params)
    print("Fitting with validation data.")
    model.fit(X_train_val, y_train_val)

    # Evaluate the final model on the test set
    accuracy = []
    print ("Calculating accuracy.")
    for key in test_keys:
        matrix = df_dict[key]
        X_test = matrix.drop("BINDING_ATOM", axis=1)
        y_test = matrix.BINDING_ATOM
        y_pred = model.predict(X_test)
        accuracy.append(accuracy_score(y_pred, y_test))
        print(f"{key} accuracy: {accuracy_score(y_pred, y_test)}")

    print(f"Average accuracy: {sum(accuracy)/len(accuracy)}")

    print("Saving model to pickle.")
    with open("model.pckl", "wb") as file:
        pickle.dump(model, file)

    print("Model saved in 'model.pckl'.")