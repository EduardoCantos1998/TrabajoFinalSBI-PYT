import pickle
import pandas as pd
import sys 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

"""
This script generates a model using the dictionary.pckl we generated with
dictionary_pickler.py. This process tends to take around 30 minutes, so do only
when the current models don't suffice. It accepts two arguments:

INPUT:
    - Weight of the value of positive cases: this affects how likely the model
      is to assign a 1 to the atom.
    - Name [OPTIONAL]: it is the name of the model pickle to be generated. 

OUTPUT:
    - Pickle of the model. It can be generated using your name as defined above.
      If left undefined, it will name it 'model_WEIGHT.pckl'.
"""

print("Loading data.")
with open("dictionary.pckl", "rb") as file:
    df_dict = pickle.load(file)

try:
    if sys.argv[1].isnumeric():
        pass
    else:
        raise ValueError("Please introduce a weight in numerical form.")
except IndexError:
    raise IndexError("Please introduce a weight in numerical form.")

print("Defining Model.")
class_weights = {0: 1, 1: sys.argv[1]}
# Split the data into training, validation and test sets
model = RandomForestClassifier(max_depth=2,random_state=0, class_weight=class_weights)

print("Spliting the keys in train/test/validation sets.")
train_keys, test_keys = train_test_split(list(df_dict.keys()), test_size=0.3, random_state=0)
train_keys, val_keys = train_test_split(train_keys, test_size=0.2, random_state=0)

# Fit the model on the training data and evaluate on the validation set
best_accuracy = 0
cur_prot = 1
print("Starting the hyperparameter tuning.")
for key in train_keys:
    print(f"Current protein: {key} ### Current protein {cur_prot:>3}/{len(df_dict.keys())}")
    cur_prot += 1
    matrix = df_dict[key]
    X_train = matrix.drop("BINDING_ATOM", axis=1)
    y_train = matrix.BINDING_ATOM

    for val_key in val_keys:
        print(f"Current validation key: {val_key}")
        val_matrix = df_dict[val_key]
        X_val = val_matrix.drop("BINDING_ATOM", axis=1)
        y_val = val_matrix.BINDING_ATOM

        # Train the model on the training set and evaluate on the validation set
        print("Fitting data.")
        model.fit(X_train, y_train)
        print("Predicting against validation.")
        y_val_pred = model.predict(X_val)
        print("Calculating accuracy.")
        val_accuracy = accuracy_score(y_val_pred, y_val)

        # Keep track of the best hyperparameters
        if val_accuracy > best_accuracy:
            print("Saving best parameters.")
            best_accuracy = val_accuracy
            best_params = model.get_params()

# Train the final model on the combined training and validation sets
print("Creating ")
train_val_keys = train_keys + val_keys
X_train_val = pd.concat([df_dict[key].drop("BINDING_ATOM", axis=1) for key in train_val_keys])
y_train_val = pd.concat([df_dict[key].BINDING_ATOM for key in train_val_keys])
model.set_params(**best_params)
print("Fitting with validation data.")
model.fit(X_train_val, y_train_val)

# Evaluate the final model on the test set
accuracy = []
print ("Calculating accuracy for evaluation set.")
for key in test_keys:
    matrix = df_dict[key]
    X_test = matrix.drop("BINDING_ATOM", axis=1)
    y_test = matrix.BINDING_ATOM
    y_pred = model.predict(X_test)
    accuracy.append(accuracy_score(y_pred, y_test))
    print(f"{key} accuracy: {accuracy_score(y_pred, y_test)}")

print(f"Average accuracy: {sum(accuracy)/len(accuracy)}")

print("Saving model to pickle.")

try:
    name = sys.argv[2]
except IndexError:
    name = f"model_{sys.argv[1]}"

with open(f"{name}.pckl", "wb") as file:
    pickle.dump(model, file)

print(f"Model saved in '{name}.pckl'.")