import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog  


# Load the uploaded CSV file
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
file_path = filedialog.askopenfilename(title="Get CSV File")
df = pd.read_csv(file_path)

# Show the first few rows and basic info
df_info = df.info()
df_head = df.head()

df_info, df_head


# Step 1: Select relevant features
features = ['pow vol diff', 'aw vol diff', 'pow sphereicity diff', 'aw sphereicity diff']
X = df[features]

# Step 2: Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Step 3: Determine optimal number of clusters using silhouette score
silhouette_scores = []
k_range = range(2, 11)  # Try K from 2 to 10

for k in k_range:
    kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
    labels = kmeans.fit_predict(X_scaled)
    score = silhouette_score(X_scaled, labels)
    silhouette_scores.append(score)

# Find the best K
best_k = k_range[np.argmax(silhouette_scores)]

# Step 4: Run KMeans with the best number of clusters
kmeans_final = KMeans(n_clusters=best_k, n_init=10, random_state=42)
df['cluster'] = kmeans_final.fit_predict(X_scaled)

# Step 5: Dimensionality reduction for visualization
pca = PCA(n_components=)
X_pca = pca.fit_transform(X_scaled)
df['pca1'] = X_pca[:, 0]
df['pca2'] = X_pca[:, 1]



# Display the clustered data
print("\nClustered Atomic Differences:")
print(df)

# Output best number of clusters
print(f"\nBest number of clusters: {best_k}")
