import data_preparation as dp
import torch
import pickle
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm

class CNN1D(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, kernel_size, stride):
        super(CNN1D, self).__init__()
        self.conv1d = nn.Conv1d(input_dim, hidden_dim, kernel_size=kernel_size, stride=stride)
        self.maxpool1d = nn.MaxPool1d(kernel_size=2, stride=2)
        self.linear1 = nn.Linear(hidden_dim * (sequence_length // 2), output_dim)

    def forward(self, x):
        x = self.conv1d(x)
        x = nn.functional.relu(x)
        x = self.maxpool1d(x)
        x = x.view(x.size(0), -1)
        x = self.linear1(x)
        return x
        
# Definir los parámetros
input_dim = 12
hidden_dim = 32
output_dim = 1
kernel_size = 3
stride = 1
learning_rate = 0.001
num_epochs = 100

# Crear el modelo y definir la función de pérdida y el optimizador
model = CNN1D(input_dim, hidden_dim, output_dim, kernel_size, stride)
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Entrenamiento
for epoch in range(num_epochs):
    train_loss = 0
    for i, data in enumerate(dp.train_loader, 0):
        inputs, labels = data
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
    train_loss /= len(dp.train_loader)
    print(f"Epoch {epoch+1}/{num_epochs}, Training Loss: {train_loss:.4f}")

#with open("model.pckl", "wb") as fl:
#    pickle.dump(model, fl)