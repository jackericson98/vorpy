import pandas as pd
import os
import matplotlib.pyplot as plt
from tkinter import Tk, Button, Label, filedialog
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score


# Function to calculate Mean Absolute Percentage Error
def mean_absolute_percentage_error(y_true, y_pred):
    y_true, y_pred = pd.Series(y_true), pd.Series(y_pred)
    return (abs((y_true - y_pred) / y_true).replace([float('inf'), -float('inf')], pd.NA).dropna().mean()) * 100


# Function to load multiple files, combine them, and train two models
def load_data_and_train():
    global model, model_pow_vol, X_test, y_test, X_test_pow_vol, y_pred, y_pred_pow_vol

    combined_data = pd.DataFrame()
    features = ['pow vol', 'pow sa', 'vdw vol', 'rad', 'pow nbors']

    while True:
        filepath = filedialog.askopenfilename(title="Select file for training",
                                              filetypes=(("CSV files", "*.csv"), ("All files", "*.*")))
        if not filepath:  # Cancelled
            break
        if os.path.exists(filepath):
            try:
                data = pd.read_csv(filepath)
                combined_data = pd.concat([combined_data, data], ignore_index=True)
                label.config(text=f'Loaded {combined_data.shape[0]} rows')
            except Exception as e:
                label.config(text=f'Error loading: {e}')

    if not combined_data.empty:
        X = combined_data[features]
        y = combined_data['aw vol']
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        X_train_pow_vol = X_train[['pow vol']]
        X_test_pow_vol = X_test[['pow vol']]

        # Train full-feature model
        model = LinearRegression()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        # Train baseline model (pow vol only)
        model_pow_vol = LinearRegression()
        model_pow_vol.fit(X_train_pow_vol, y_train)
        y_pred_pow_vol = model_pow_vol.predict(X_test_pow_vol)

        # Calculate metrics
        r2 = r2_score(y_test, y_pred)
        r2_pow = r2_score(y_test, y_pred_pow_vol)

        mape = mean_absolute_percentage_error(y_test, y_pred)
        mape_pow = mean_absolute_percentage_error(y_test, y_pred_pow_vol)

        label.config(text=f'R²: {r2:.3f} ({r2_pow:.3f}), MAPE: {mape:.2f}% ({mape_pow:.2f}%)')
    else:
        label.config(text="No data loaded.")
        model = model_pow_vol = None


def predict_from_loaded_data():
    if model:
        try:
            y_pred = model.predict(X_test)
            mse = mean_squared_error(y_test, y_pred)
            label.config(text=f'Model MSE: {mse:.2f}')
        except Exception as e:
            label.config(text=f'Prediction error: {e}')
    else:
        label.config(text='Train the model first.')


def print_results():
    if model and y_pred is not None:
        print("Actual:", y_test.to_numpy())
        print("Predicted:", y_pred)
    else:
        print("Model not trained or predictions not made.")


def plot_results():
    if model and y_pred is not None:
        plt.figure(figsize=(6, 6))
        plt.scatter(y_test, y_pred, alpha=0.7)
        plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--')
        plt.xlabel('Actual aw vol')
        plt.ylabel('Predicted aw vol')
        plt.title('Actual vs Predicted')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    else:
        print("Model not trained or predictions not made.")


# Initialize GUI
root = Tk()
root.title("aw vol Predictor")
root.geometry("400x350")

# Global variables
model = model_pow_vol = None
X_test = y_test = X_test_pow_vol = y_pred = y_pred_pow_vol = None

# GUI elements
label = Label(root, text="Load CSV files to train model", font=("Helvetica", 12))
label.pack(pady=20)

train_button = Button(root, text="Load Data & Train Model", command=load_data_and_train)
train_button.pack(pady=5)

predict_button = Button(root, text="Predict from Latest Data", command=predict_from_loaded_data)
predict_button.pack(pady=5)

print_button = Button(root, text="Print Results", command=print_results)
print_button.pack(pady=5)

plot_button = Button(root, text="Plot Results", command=plot_results)
plot_button.pack(pady=5)

# Run
root.mainloop()
