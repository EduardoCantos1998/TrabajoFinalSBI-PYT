from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import sys

protein = sys.argv[1]
# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure('1XYZ', protein)

# Create a DSSP object
for i in structure:
    model = i
    dssp = DSSP(model, protein)
    # Iterate over each residue and print its secondary structure
    for residue in dssp:
        res_id = residue[0]
        ss = residue[2]
        print(f"Residue {res_id} is in {ss} structure.")