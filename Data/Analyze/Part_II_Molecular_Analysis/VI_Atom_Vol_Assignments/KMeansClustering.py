import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
import tkinter as tk
from tkinter import filedialog

# Load Data
root = tk.Tk()
root.withdraw()
root.wm_attributes()
file_path = filedialog.askopenfilename()
df = pd.read_csv(file_path)

# Separate hydrogen atoms for plotting but not clustering
hydrogens = df[df['name'].str.startswith('H')]
df = df[~df['name'].str.startswith('H')]

# Select relevant features for clustering
features = ['pow sphereicity diff', 'pow vol diff']
X = df[features]

# Standardize numerical features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Encode categorical atom names to help clustering
encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
X_names_encoded = encoder.fit_transform(df[['name']])

# Combine numerical and categorical features
X_augmented = np.hstack((X_scaled, X_names_encoded))

# Run K-Means clustering with 12 clusters
best_k = 12
kmeans = KMeans(n_clusters=best_k, n_init=10, random_state=42)
df['cluster'] = kmeans.fit_predict(X_augmented)

# Save the full dataset with cluster labels and coordinates
df.to_csv("atomic_comparisons_with_clusters.csv", index=False)

# Create summary CSV with counts of each atom/residue per cluster, including coordinates
cluster_summary = (
    df.groupby(['cluster', 'name', 'residue'])
    .size()
    .reset_index(name='count')
)

# Add coordinates to the cluster summary if they exist in the dataset
if 'x' in df.columns and 'y' in df.columns and 'z' in df.columns:
    cluster_summary = cluster_summary.merge(
        df[['cluster', 'name', 'residue', 'x', 'y', 'z']].drop_duplicates(),
        on=['cluster', 'name', 'residue'],
        how='left'
    )

# Save the full dataset with cluster labels and coordinates
output_folder = os.path.dirname(file_path)
df.to_csv(output_folder + "/atomic_comparisons_with_clusters.csv", index=False)

# Plot the clusters
plt.figure(figsize=(10, 8))
ax = plt.gca()

color_map = {'C': 'grey', 'N': 'blue', 'P': 'orange', 'O': 'red', 'S': 'yellow'}
df['element_color'] = df['name'].str[0].map(color_map).fillna('black')

sns.scatterplot(
    data=hydrogens,
    x='pow sphereicity diff',
    y='pow vol diff',
    color='pink',
    s=20,
    alpha=0.5,
    ax=ax
)

sns.scatterplot(
    data=df,
    x='pow sphereicity diff',
    y='pow vol diff',
    hue=df['element_color'],  # Color by first letter of atom name (element)
    palette=df['element_color'].unique(),
    s=60,
    alpha=0.8,
    ax=ax
)

# Draw bounding boxes around clusters
for cluster_id, group in df.groupby('cluster'):
    x_center = group['pow sphereicity diff'].mean()
    y_center = group['pow vol diff'].mean()
    x_range = (group['pow sphereicity diff'].max() - group['pow sphereicity diff'].min()) * 2/3
    y_range = (group['pow vol diff'].max() - group['pow vol diff'].min()) * 2/3
    xmin, ymin = x_center - x_range / 2, y_center - y_range / 2
    xmax, ymax = x_center + x_range / 2, y_center + y_range / 2
    rect = patches.Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin,
        linewidth=1.5, edgecolor='black', facecolor='none', linestyle='--'
    )
    ax.add_patch(rect)
    ax.text(xmin, ymax, str(cluster_id), fontsize=9, color='black')

plt.title(f'p53tet Clustering', fontsize=36)
plt.xlabel('Power Sphericity Difference', fontsize=36)
plt.ylabel('Power Volume Difference', fontsize=36)
handles, labels = ax.get_legend_handles_labels()
legend_labels = {v: k for k, v in color_map.items() if v in labels}
new_labels = [legend_labels[label] if label in legend_labels else label for label in labels]
ax.legend(handles, new_labels, bbox_to_anchor=(1.05, 1), loc='upper left', title='Element', fontsize=27, title_fontsize=30)
plt.grid(False)
plt.tight_layout()
plt.show()
