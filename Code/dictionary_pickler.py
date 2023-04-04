import mol2
import os
import df_maker
from random import sample
import pickle
import warnings

warnings.filterwarnings("ignore")

# Generate the files
files = "../Data/final_data/"
structures={}

for i in os.listdir(files):
    structures[i] = ['cavityALL.mol2', 'protein.mol2', 'site.mol2', 'ligand.mol2']

# Generate the random files sample.
random_files = sample(list(structures.keys()), 2000)

print("Generating DataFrame dictionary.")
df_dict = {}
time = 1
errors = 0
successes = 0
for i in random_files:
    print(f"############################################ Successes = {successes} Errors = {errors}")
    print(f"Protein {i}:{time:>5}/{len(random_files)}")
    time += 1
    try:
        path = f"../Data/final_data/{i}/"
        my_protein = mol2.Protein(name=f"i", protein=f"{path}protein.mol2", cavity=f"{path}cavity6.mol2", ligand=f"{path}ligand.mol2", site=f"{path}site.mol2")
        df_dict[i] = df_maker.NewDataFrame(my_protein, f"{path}protein.pdb")
        successes += 1
    except ValueError:
        print("Value error.")
        errors += 1
    except FileNotFoundError:
        print("File not found.")
        errors += 1
    except KeyError:
        print("CA not found.")
        errors += 1
    except Exception:
        print("DSSP failed")
        errors += 1

print("##########################################")
print(f"Errors: {errors}/{len(random_files)}")
print(f"Successes: {successes}/{len(random_files)}")
print(f"Error/Success ratio: {errors/successes:.3f}")
print("##########################################")

print("Saving dictionary to pickle.")
with open("dictionary.pckl", "wb") as file:
    pickle.dump(df_dict, file)
print("File saved to 'dictionary.pckl'.")