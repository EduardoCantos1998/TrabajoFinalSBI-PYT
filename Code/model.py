import data_preparation as dp
import torch
from tqdm import tqdm
import pickle

model = MyModel()

num_epochs = 10
criterion = torch.nn.MSELoss()  # Función de pérdida
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)  # Optimizador

for epoch in range(num_epochs):
    # Entrenamiento
    model.train()
    for batch in dp.train_loader:
        inputs = batch
        labels = ... # Calculamos las etiquetas para el conjunto de entrenamiento
        outputs = model(inputs)
        loss = criterion(outputs, labels)  # Calculamos la función de pérdida
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()

    # Validación
    model.eval()
    with torch.no_grad():
        for batch in dp.val_loader:
            inputs = batch
            labels = ... # Calculamos las etiquetas para el conjunto de validación
            outputs = model(inputs)
            val_loss = criterion(outputs, labels)  # Calculamos la función de pérdida en el conjunto de validación

    # Imprimimos la tasa de pérdida en cada época
    print(f"Epoch {epoch}, train loss: {train_loss/len(dp.train_loader)}, val loss: {val_loss/len(dp.val_loader)}")

with open("model.p", "wb") as fl:
    pickle.dump(model, fl)