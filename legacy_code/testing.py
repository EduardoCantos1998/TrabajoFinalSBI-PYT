import mol2
import features as ft
import numpy as np

path = "../Data/final_data/1w5a_1/"

protein = mol2.Protein(name="1w5a_1", protein = f"{path}protein.mol2",cavity=f"{path}cavity6.mol2",site=f"{path}site.mol2",ligand=f"{path}ligand.mol2")

print(protein.get_site().shape)
