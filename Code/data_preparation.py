import torch
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
import pickle

with open("features.pckl", "rb") as f:
    features = pickle.load(f)

print(features["2ijg_2"]["site_ligand_DISTANCES"])

# for i in features:
#    for j in features[i]:
#        print(f"Protein '{i}', Feature '{j}': ",features[i][j].numel())

class ProteinDataset(Dataset):
    def __init__(self, features):
        self.features = features
        self.proteins = list(self.features.keys())

    def __len__(self):
        return len(self.proteins)

    def __getitem__(self, idx):
        protein = self.proteins[idx]
        sample = {
            "protein_TENSOR": self.features[protein]["protein_TENSOR"],
            "cavity_TENSOR": self.features[protein]["cavity_TENSOR"],
            "site_TENSOR": self.features[protein]["site_TENSOR"],
            "ligand_TENSOR": self.features[protein]["ligand_TENSOR"],
            #"proteinCA_DISTANCES": self.features[protein]["proteinCA_DISTANCES"],
            #"cavity_DISTANCES": self.features[protein]["cavity_DISTANCES"],
            #"siteCA_DISTANCES": self.features[protein]["siteCA_DISTANCES"],
            #"site_ligand_DISTANCES": self.features[protein]["site_ligand_DISTANCES"],
            #"protein_PHI": self.features[protein]["protein_PHI"],
            #"protein_PSI": self.features[protein]["protein_PSI"],
            #"site_PHI": self.features[protein]["site_PHI"],
            #"site_PSI": self.features[protein]["site_PSI"]
        }
        return sample

# Dividimos los datos en conjuntos de entrenamiento, validación y prueba
train_data, test_data = train_test_split(list(features.keys()), test_size=0.2, random_state=42)
train_data, val_data = train_test_split(train_data, test_size=0.2, random_state=42)

# Creamos las instancias de ProteinDataset para cada conjunto
train_set = ProteinDataset({k: v for k, v in features.items() if k in train_data})
val_set = ProteinDataset({k: v for k, v in features.items() if k in val_data})
test_set = ProteinDataset({k: v for k, v in features.items() if k in test_data})

# Creamos las instancias de DataLoader para cada conjunto
train_loader = DataLoader(train_set, batch_size=32, shuffle=True)
val_loader = DataLoader(val_set, batch_size=32, shuffle=False)
test_loader = DataLoader(test_set, batch_size=32, shuffle=False)
