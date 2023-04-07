# SBIXPYT: RF approach for PPI

**Authors**: Espitia S., Gary; Marco D., Alejandro; Cantos G., Eduardo

![Results from 3IMX](3IMX_result.png?raw=true "Results from 3IMX")

## Table of Contents
- [SBIXPYT: RF approach for PLI](#sbixpyt-rf-approach-for-pli)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Training](#training)
  - [Script description](#script-description)
  - [Requirements](#requirements)
  - [Command line Installation](#command-line-installation)
  - [Usage (Tutorial)](#usage-tutorial)
    - [Running the code](#running-the-code)
  - [Output](#output)
  - [Theory](#theory)
  - [Result Benchmark (Analysis)](#result-benchmark-analysis)
  - [License](#license)
  - [References](#references)

## Introduction 
This project is designed to determine the binding site of proteins using Random Forest (RF).
The program takes a PDB file and generates an output called output.pdb that has information that can be visualized in [Jmol](https://jmol.sourceforge.net/), [PyMOL](https://pymol.org/2/), [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/index.html), etc.

## Training
We will use the subset of the scPDB dataset generated by the [PUResNet](#references) team as an starting point. Particularly we used the following data from their [repository](https://github.com/jivankandel/PUResNet): `protein.mol2`, `site.mol2`. We added the .pdb files using `mol2pdb.py` and used them for training.
We will also use the BindingDB as a simple visual validation set.
We extracted the PDB files from the datasets with only PDBs included in articles, the subset of files drawn from the ChEMBL, and also the subset of the files from patents and those published in papers. 

## Script description
- `Code/mol2.py`: `get_protein`, `get_proteinCA`, and `get_proteinCB` methods extract the coordinates of all atoms, alpha-carbons (CA), and beta-carbons (CB), `get_cavity`, `get_site`, `get_siteCA`, `get_siteCB` in binding cavity, binding site, and ligand of the protein from the PDB file, and return them as NumPy matrices.
which downloads PDB (Protein Data Bank) and FASTA files for a list of protein codes, which are present in the `"final_data/"` directory.
- `mol2pdb.py`: Extracts the IDs from each folder and retrieves the pdb file from [RCSB](https://files.rcsb.org/).
- `df_maker.py`: extracts features (coordinates, aminoacid, binding atom, entropy, charge, hidrophobicity, secundary structure, solvent accessible surface area(SASA), b-factor, phi and psi angles and alpha-carbon distance). betacarbons; and generates a pandas dataframe and converts them to integers.
- `dictionary_pickler.py`: it iterates the proteins to generate a dataframe using `df_maker.py`; and keeps them in a dictionary. It could be considered the first selection step because we create a sample from the database (5020).
- `DATA.py`: this is used to extract the PDB codes from the BindingDB tsv and download into a folder in the same directory called "PDB". A zip of this can be found in the DATA folder.

## Requirements
This is an Python script that particularly uses the following dependencies to take into account: biopython, df_maker, freesasa, mol2, numpy, pandas, scikit-learn and DSSP. 

:warning: Particularly DSSP is meant to be ran in a **Ubuntu** or **Mac OSX**; at the moment of this release DSSP may not work in other distributions or operating systems.


## Command line Installation

```bash
git clone https://github.com/EduardoCantos1998/TrabajoFinalSBI-PYT
cd TrabajoFinalSBI-PYT
```

and proceed to create a python env to run the scripts.

```bash
# Create a virtual environment:
python -m venv venv # Or name it as desired 
source venv/bin/activate
pip install -r requirements.txt
```

**To install DSSP follow these [instructions](https://ssbio.readthedocs.io/en/latest/instructions/dssp.html#installation-instructions-ubuntu).**


## Usage (Tutorial)
This is the workflow for the general use of the tool:

```mermaid
graph  LR
A1[PUResNet data] --> B1[mol2pdb.py] 
B1 --> G
subgraph OPTIONALLY: creating your own pickle fa:fa-jar with your own data;
  subgraph invoked;
    direction LR;
    B[model.py]; 
    C[mol2.py]-.-o E;
    D[df_maker.py]-.-o E;
    E[dictionary_pickler.py]-.- dictionary.pckl -.-o B;
  end
end
G(((custom data))) --> E;
invoked -.-o B
  A(((input.pdb))) --> F;
subgraph testing pdb with our data;
  F[pdb_testing.py];
end;
B -.-model.pckl-.-o F;
F ==> Z(((output.pdb)))
D --> F
```
It takes as an input a PDB file which is evaluated using `model.py` then; the output will be a list of the aminoacids and sites belonging to a binding site. 

### Running the code
We go to the folder
```bash
cd Code
python3 pdb_testing.py [PATH_TO_PROTEIN]
```
In case the user wants to use its own model; as an option is also possible to:
```bash
python3 pdb_testing.py [PATH_TO_PROTEIN] [PATH_TO_MODEL]
```

## Output
This is an example of the output, as a concept:

```
binding_site_prediction = [0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1]
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N']

# Obtener los aminoácidos que corresponden con el binding site
binding_site_amino_acids = [amino_acids[i] for i, val in enumerate(binding_site_prediction) if val == 1]

print(binding_site_amino_acids)
# Output: ['D', 'E', 'G', 'H', 'M', 'N']
```
And the following output files:
- `{name}_binding_site_predictions.txt`
- `{name}_prediction.pdb`

Displayed in UCSFChimera:
![Results from 5T2W](5T2W_result.png?raw=true "prediction in blue")

_**Fig. 1**: This is an visualization of the results, being beige the local sequence; being blue the prediction and red the XFC ligand for this interaction._

## [Theory](theory.md)
## Result Analysis
The models were trained with different weights for the positive binding site value. This was necessary since the data for the binding sites was unbalanced in favor of the negative value for binding sites. We obtained 4 main models. The main models are the ones with a weight of 6 and 7, and an accuracy of 86.11% and 76.41% respectively. They obtained the highest accuracy while training, but also gave the best predictions during our own visual testing. The other two models have a weight of 8 and 10, with an accuracy of 66.65% and 60% respectively. These last two models had a really low accuracy, and they didn’t prove to be useful when predicting the binding site. We can see that the optimal weight is around 6. A weight of 5 would have been too low, since 6 is already giving really few positive binding sites results. We encourage the user to train the model with a lower weight if predictions are unfitting. There might be some cases where the prediction might be too poor or too generous. For situations like this one it would then be best to develop a new model and use it in their prediction. When the weight is not specified, it predicts all the atoms are not binding sites. This also resulted in a 95% accuracy, which came as a surprise since we know that all the proteins have a binding site. This is what leads us to generate different models with different weights (Fig. 2).

![Model accuracy v. weight](results_graph.jpeg?raw=true "Model accuracy in weight")

_**Fig. 2**: We can observe that as the weight increases, the accuracy also does as well. For weights lower than 6 we did observe that the accuracy was higher, but that didn’t correspond to better results. We observed the optimum weight to be around 6._

## [License](LICENSE)

## References
1. Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, et al. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research. 2011;12(85):2825–30. 
2. Desaphy J, Bret G, Rognan D, Kellenberger E. sc-PDB: a 3D-database of ligandable binding sites—10 years on. Nucleic Acids Research. 2015 Jan 28;43(D1):D399–404. 
3. Ho TK. Random decision forests. In: Proceedings of 3rd International Conference on Document Analysis and Recognition. 1995. p. 278–82 vol.1. 
4. Kandel J, Tayara H, Chong KT. PUResNet: prediction of protein-ligand binding sites using deep residual neural network. Journal of Cheminformatics. 2021 Sep 8;13(1):65. 
5. Raza K. Protein features identification for machine learning-based prediction of protein-protein interactions [Internet]. Bioinformatics; 2017 May [cited 2023 Apr 5]. Available from: http://biorxiv.org/lookup/doi/10.1101/137257
6. Šikić M, Tomić S, Vlahoviček K. Prediction of Protein–Protein Interaction Sites in Sequences and 3D Structures by Random Forests. Stormo GD, editor. PLoS Comput Biol. 2009 Jan 30;5(1):e1000278. 
7. Ma W, Bao W, Cao Y, Yang B, Chen Y. Prediction of Protein-Protein Interaction Based on Deep Learning Feature Representation and Random Forest. In: Huang DS, Jo KH, Li J, Gribova V, Premaratne P, editors. Intelligent Computing Theories and Application. Cham: Springer International Publishing; 2021. p. 654–62. (Lecture Notes in Computer Science). 
8. pickle — Python object serialization [Internet]. Python documentation. [cited 2023 Apr 6]. Available from: https://docs.python.org/3/library/pickle.html
9. Casadio R, Martelli PL, Savojardo C. Machine learning solutions for predicting protein–protein interactions. WIREs Comput Mol Sci [Internet]. 2022 Nov [cited 2023 Apr 5];12(6). Available from: https://onlinelibrary.wiley.com/doi/10.1002/wcms.1618
10. Mahesh B. Machine Learning Algorithms -A Review. 2019. 
11. Kabsch W, Sander C. Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983 Dec;22(12):2577–637. 
12. Jamasb AR, Day B, Cangea C, Liò P, Blundell TL. Deep Learning for Protein–Protein Interaction Site Prediction. In: Cecconi D, editor. Proteomics Data Analysis [Internet]. New York, NY: Springer US; 2021 [cited 2023 Apr 5]. p. 263–88. (Methods in Molecular Biology). Available from: https://doi.org/10.1007/978-1-0716-1641-3_16
13. Das S, Chakrabarti S. Classification and prediction of protein–protein interaction interface using machine learning algorithm. Sci Rep. 2021 Jan 19;11(1):1761. 
14. ChEMBL Database [Internet]. [cited 2023 Apr 7]. Available from: https://www.ebi.ac.uk/chembl/
15. 3.3. Metrics and scoring: quantifying the quality of predictions [Internet]. scikit-learn. [cited 2023 Apr 6]. Available from: https://scikit-learn/stable/modules/model_evaluation.html
