
import os
import urllib.request
import warnings

warnings.filterwarnings("ignore")

print("Reading files.")
files = os.listdir("final_data/")

print("Starting download.")
cur_prot = 0
print(f"############################## {cur_prot}/{len(files)}")
try:
    for i in files:
        code = i.split(sep="_")[0].upper()
        pdb_url = f"https://files.rcsb.org/download/{code}.pdb"
        pdb_path = f"final_data/{i}/protein.pdb"
        fa_url = f"https://www.rcsb.org/fasta/entry/{code}"
        fa_path = f"final_data/{i}/protein.fa"
        if f"protein.pdb" not in os.listdir(f"final_data/{i}/"):
            urllib.request.urlretrieve(pdb_url, pdb_path)
            print(f"PDB file {code} downloaded and saved to {pdb_path}")
            urllib.request.urlretrieve(fa_url, fa_path)
            print(f"Fasta file {code} downloaded and saved to {fa_path}")
        else:
            print(f"PDB file {code} already found in {pdb_path}")
        print(f"############################## {cur_prot}/{len(files)}")
        cur_prot += 1
except urllib.error.HTTPError:
    print(f"File {code} not found.")
