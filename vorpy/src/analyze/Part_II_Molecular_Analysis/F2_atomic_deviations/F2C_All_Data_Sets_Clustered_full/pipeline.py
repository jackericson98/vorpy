#
# import os
# import sys
# import tkinter as tk
# from tkinter import filedialog
# from typing import Dict, List, Optional
#
# import matplotlib.pyplot as plt
# import matplotlib.patheffects as pe
# from matplotlib.patches import Ellipse
# import numpy as np
# import pandas as pd
# # from sklearn.decomposition import PCA
# # from sklearn.preprocessing import StandardScaler
#
# from .config import (
#     SMALL_MOL_ATOM_ALIASES,
#     PROTEIN_ATOM_ALIASES,
#     RNA_ATOM_ALIASES,
#     DNA_ATOM_ALIASES,
#     DIRECT_LABEL_GROUPS,
#     FORCE_INCLUDE_GROUPS,
# )
#
# # Get the path to the root vorpy folder
# vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..', '..'))
# if vorpy_root not in sys.path:
#     sys.path.append(vorpy_root)
#
# from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
# from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import compare_manual_vs_ml
# from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import make_sol_facing_binary
# from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import print_cluster_summary
# from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import run_clustering
# from vorpy.src.analyze.Part_II_Molecular_Analysis.F2_atomic_deviations.F2C_Helper_Auto_Cluster import summarize_clusters
#
#
# def run_pca_analysis(df: pd.DataFrame, numeric_cols: List[str], label: str = ""):
#     print(f"\n=== PCA ANALYSIS: {label} ===")
#
#     if not numeric_cols:
#         print("No numeric columns supplied for PCA.")
#         return None, None
#
#     work_df = df.dropna(subset=numeric_cols).copy()
#
#     if len(work_df) == 0:
#         print("No data available for PCA.")
#         return None, None
#
#     X = work_df[numeric_cols].values
#     # X_scaled = StandardScaler().fit_transform(X)
#     #
#     # pca = PCA()
#     pca.fit(X_scaled)
#
#     explained = pca.explained_variance_ratio_
#
#     print("\nExplained variance ratio:")
#     for i, var in enumerate(explained, start=1):
#         print(f"PC{i}: {var:.4f}")
#
#     print("\nCumulative variance:")
#     print(np.cumsum(explained))
#
#     loadings = pd.DataFrame(
#         pca.components_.T,
#         columns=[f"PC{i+1}" for i in range(len(numeric_cols))],
#         index=numeric_cols,
#     )
#
#     print("\nPCA loadings:")
#     print(loadings.to_string())
#
#     return pca, loadings
#
#
# def summarize_sol_facing(df: pd.DataFrame):
#     if 'SolFacingPct' not in df.columns:
#         print("\nSolFacingPct column missing.")
#         return
#
#     print("\n=== SOL-FACING SUMMARY ===")
#     print(df['SolFacingPct'].describe())
#     print("Non-null count:", int(df['SolFacingPct'].notna().sum()))
#
#     if 'SolFacingBinary' in df.columns:
#         print("\nSolFacingBinary counts:")
#         print(df['SolFacingBinary'].value_counts(dropna=False).to_string())
#
#
# def resolve_group_folders(data_root: str, relative_folders: List[str]) -> List[str]:
#     folders = []
#
#     for rel_path in relative_folders:
#         full_path = os.path.join(data_root, rel_path)
#         if os.path.isdir(full_path):
#             folders.append(full_path)
#         else:
#             print(f"WARNING: folder not found -> {full_path}")
#
#     return folders
#
#
# def canonicalize_atom_name(
#     atom_name: str,
#     molecule_class: str,
#     residue_name: str = ''
# ) -> str:
#     atom = str(atom_name).strip().upper()
#     mol = str(molecule_class).strip().lower()
#
#     if mol == 'protein':
#         return PROTEIN_ATOM_ALIASES.get(atom, atom)
#
#     if mol == 'dna':
#         return DNA_ATOM_ALIASES.get(atom, atom)
#
#     if mol == 'rna':
#         return RNA_ATOM_ALIASES.get(atom, atom)
#
#     if mol == 'small_molecule':
#         return SMALL_MOL_ATOM_ALIASES.get(atom, atom)
#
#     return atom
#
#
# def infer_molecule_class(atoms_df: pd.DataFrame) -> str:
#     residue_names = set(
#         str(x).strip().upper()
#         for x in atoms_df['Residue'].dropna().unique()
#     )
#
#     dna_residues = {'DA', 'DC', 'DG', 'DT', 'DI'}
#     rna_residues = {'A', 'C', 'G', 'U', 'RA', 'RC', 'RG', 'RU'}
#     protein_markers = {'ALA', 'GLY', 'VAL', 'LEU', 'SER', 'THR', 'ASP', 'GLU', 'LYS', 'ARG'}
#
#     if residue_names & protein_markers:
#         return 'protein'
#
#     if residue_names & dna_residues:
#         return 'dna'
#
#     if residue_names & rna_residues:
#         return 'rna'
#
#     return 'small_molecule'
#
#
# def _read_scheme_logs(folder: str):
#     """
#     Reads aw/pow/prm logs from either:
#       folder/{aw_logs.csv,pow_logs.csv,prm_logs.csv}
#     or
#       folder/{aw/aw_logs.csv, pow/pow_logs.csv, prm/prm_logs.csv}
#     """
#     try:
#         aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
#         pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
#         prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)
#
#     except FileNotFoundError:
#         aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
#         pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
#         prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)
#
#     return aw_logs, pow_logs, prm_logs
#
#
# def select_folders_multi(title: str = "Select a folder (Cancel to finish)") -> List[str]:
#     root = tk.Tk()
#     root.withdraw()
#
#     folders = []
#
#     while True:
#         folder = filedialog.askdirectory(title=title)
#         if not folder:
#             break
#         folders.append(folder)
#
#     root.destroy()
#     return folders
#
#
# def collect_atom_volume_points(
#     folders: List[str],
#     atom_name_field: str = 'Name',
#     molecule_class: str = 'protein',
#     volume_range: Optional[tuple] = None
# ) -> pd.DataFrame:
#     records: List[Dict[str, object]] = []
#
#     for folder in folders:
#         aw_logs, pow_logs, _ = _read_scheme_logs(folder)
#
#         aw_atoms = aw_logs['atoms']
#         pow_atoms = pow_logs['atoms']
#
#         pow_lookup = {
#             int(row['Index']): row
#             for _, row in pow_atoms.iterrows()
#         }
#
#         for _, atom in aw_atoms.iterrows():
#             idx = int(atom['Index'])
#             chain = atom['Chain']
#             res_seq = atom['Residue Sequence']
#
#             if idx not in pow_lookup:
#                 continue
#
#             pow_atom = pow_lookup[idx]
#
#             aw_v = float(atom['Volume'])
#             pow_v = float(pow_atom['Volume'])
#
#             if volume_range is not None:
#                 vmin, vmax = volume_range
#                 if aw_v < vmin or aw_v > vmax or pow_v < vmin or pow_v > vmax:
#                     continue
#
#             atom_name = str(atom.get(atom_name_field, '')).strip().upper()
#
#             if not atom_name:
#                 continue
#
#             residue_name = str(atom.get('Residue', atom.get('Residue Name', atom.get('ResName', '')))).strip().upper()
#             canonical_name = canonicalize_atom_name(
#                 atom_name=atom_name,
#                 molecule_class=molecule_class,
#                 residue_name=residue_name
#             )
#
#             records.append({
#                 'Folder': folder,
#                 'Index': idx,
#                 'Chain': chain,
#                 'ResidueName': residue_name,
#                 'Residue Sequence': res_seq,
#                 'AtomName': atom_name,
#                 'CanonicalName': canonical_name,
#                 'AW': aw_v,
#                 'Pow': pow_v,
#                 'x': float(atom['x']) if 'x' in atom.index else np.nan,
#                 'y': float(atom['y']) if 'y' in atom.index else np.nan,
#                 'SolFacingPct': float(atom['SolFacingPct']) if 'SolFacingPct' in atom.index else np.nan,
#             })
#
#     return pd.DataFrame(records)
#
#
# def analyze_group_outliers(df, min_group_size=3, z_thresh=2.5, top_k=20):
#     """
#     Identify:
#     - Small / ungrouped clusters
#     - Within-group outliers
#     - Global outliers
#
#     Returns:
#         dict of DataFrames
#     """
#
#     results = {}
#
#     # -----------------------
#     # 1. Group stats
#     # -----------------------
#     group_stats = df.groupby('CanonicalName').agg(
#         n=('AW', 'count'),
#         aw_mean=('AW', 'mean'),
#         pow_mean=('Pow', 'mean'),
#         aw_std=('AW', 'std'),
#         pow_std=('Pow', 'std')
#     ).reset_index()
#
#     results['group_stats'] = group_stats
#
#     # -----------------------
#     # 2. Small / ungrouped groups
#     # -----------------------
#     small_groups = group_stats[group_stats['n'] <= min_group_size].copy()
#     results['small_groups'] = small_groups.sort_values('n')
#
#     # -----------------------
#     # 3. Within-group z-score outliers
#     # -----------------------
#     df = df.copy()
#
#     df = df.merge(group_stats, on='CanonicalName', how='left')
#
#     df['z_aw'] = (df['AW'] - df['aw_mean']) / (df['aw_std'].replace(0, np.nan))
#     df['z_pow'] = (df['Pow'] - df['pow_mean']) / (df['pow_std'].replace(0, np.nan))
#
#     df['z_total'] = np.sqrt(df['z_aw']**2 + df['z_pow']**2)
#
#     within_outliers = df[df['z_total'] > z_thresh].copy()
#     results['within_group_outliers'] = within_outliers.sort_values('z_total', ascending=False)
#
#     global_mean = df[['AW', 'Pow']].mean().values
#
#     df['global_dist'] = np.sqrt(
#         (df['AW'] - global_mean[0])**2 +
#         (df['Pow'] - global_mean[1])**2
#     )
#
#     global_outliers = df.sort_values('global_dist', ascending=False).head(top_k)
#     results['global_outliers'] = global_outliers
#
#     results['annotated_df'] = df
#
#     return results
#
#
# def find_ungrouped_atoms(
#     atom_df: pd.DataFrame,
#     aw_range: Optional[tuple] = None,
#     pow_range: Optional[tuple] = None,
#     top_k: int = 50
# ):
#     """
#     Find atoms whose canonical name is unchanged, meaning they were not grouped
#     by the alias dictionary.
#
#     Optional AW / Pow window lets you inspect a suspicious patch.
#     """
#     df = atom_df.copy()
#
#     # ungrouped = canonicalization did nothing
#     df = df[df['CanonicalName'] == df['AtomName']].copy()
#
#     if aw_range is not None:
#         df = df[(df['AW'] >= aw_range[0]) & (df['AW'] <= aw_range[1])]
#
#     if pow_range is not None:
#         df = df[(df['Pow'] >= pow_range[0]) & (df['Pow'] <= pow_range[1])]
#
#     if len(df) == 0:
#         print("\nNo ungrouped atoms found in this selection.")
#         return df, pd.DataFrame(), pd.DataFrame()
#
#     summary_by_atom = (
#         df.groupby('AtomName')
#         .agg(
#             Count=('AtomName', 'size'),
#             Mean_AW=('AW', 'mean'),
#             Mean_Pow=('Pow', 'mean'),
#             Residues=('ResidueName', lambda x: ', '.join(sorted(set(map(str, x))))),
#         )
#         .reset_index()
#         .sort_values(['Count', 'Mean_AW', 'Mean_Pow'], ascending=[False, True, True])
#     )
#
#     summary_by_atom_residue = (
#         df.groupby(['AtomName', 'ResidueName'])
#         .agg(
#             Count=('AtomName', 'size'),
#             Mean_AW=('AW', 'mean'),
#             Mean_Pow=('Pow', 'mean'),
#         )
#         .reset_index()
#         .sort_values(['Count', 'Mean_AW', 'Mean_Pow'], ascending=[False, True, True])
#     )
#
#     print("\n=== UNGROUPED ATOMS: BY ATOM NAME ===")
#     print(summary_by_atom.head(top_k).to_string(index=False))
#
#     print("\n=== UNGROUPED ATOMS: BY ATOM NAME + RESIDUE ===")
#     print(summary_by_atom_residue.head(top_k).to_string(index=False))
#
#     return df, summary_by_atom, summary_by_atom_residue
#
#
# def extract_patch_atoms(
#     atom_df: pd.DataFrame,
#     aw_range: tuple,
#     pow_range: tuple,
#     min_count: int = 1
# ):
#     """
#     Extract atoms in a specific AW/Pow region and print
#     PyMOL-ready index selections.
#     """
#
#     df = atom_df.copy()
#
#     # filter by region
#     df = df[
#         (df['AW'] >= aw_range[0]) & (df['AW'] <= aw_range[1]) &
#         (df['Pow'] >= pow_range[0]) & (df['Pow'] <= pow_range[1])
#     ].copy()
#
#     if len(df) == 0:
#         print("No atoms found in patch.")
#         return df
#
#     # summary
#     summary = (
#         df.groupby(['AtomName', 'ResidueName'])
#         .size()
#         .reset_index(name='Count')
#         .sort_values('Count', ascending=False)
#     )
#
#     summary = summary[summary['Count'] >= min_count]
#
#     print("\n=== PATCH SUMMARY ===")
#     print(summary.to_string(index=False))
#
#     # print raw atoms
#     print("\n=== PATCH ATOMS ===")
#     cols = ['Index', 'AtomName', 'ResidueName', 'Chain', 'Residue Sequence', 'AW', 'Pow']
#     cols = [c for c in cols if c in df.columns]
#
#     print(df[cols].sort_values(['ResidueName', 'Residue Sequence']).to_string(index=False))
#
#     # build PyMOL selections
#     print("\n=== PYMOL SELECTIONS ===")
#
#     # 1. Precise atom-level selection using chain / resi / name
#     required_cols = {'AtomName', 'ResidueName', 'Index'}
#     has_chain = 'Chain' in df.columns
#     has_resi = 'Residue Sequence' in df.columns
#
#     if has_chain and has_resi and 'AtomName' in df.columns:
#         selection_terms = []
#
#         for _, row in df.sort_values(['Chain', 'Residue Sequence', 'AtomName']).iterrows():
#             chain = str(row['Chain']).strip()
#             resi = str(row['Residue Sequence']).strip()
#             atom_name = str(row['AtomName']).strip()
#
#             # PyMOL atom selection: chain X and resi 12 and name C4
#             term = f"(chain {chain} and resi {resi} and name {atom_name})"
#             selection_terms.append(term)
#
#         selection_str = " or ".join(selection_terms)
#
#         print("\nPyMOL precise selection:")
#         print(f"select patch_atoms, ({selection_str})")
#
#     # 2. Optional debug print using your internal Index shifted to PyMOL-style +1
#     elif 'Index' in df.columns:
#         pymol_indices = sorted((df['Index'] + 1).unique())
#         index_str = " or ".join([f"index {i}" for i in pymol_indices])
#
#         print("\nPyMOL index-only selection (+1 corrected):")
#         print(f"select patch_atoms_idx, ({index_str})")
#
#     return df
#
#
# def compute_name_volume_stats(df: pd.DataFrame, group_col: str = 'CanonicalName') -> pd.DataFrame:
#     if len(df) == 0:
#         return pd.DataFrame(
#             columns=[
#                 'GroupName', 'Count', 'Mean_AW', 'Mean_Pow',
#                 'SD_AW', 'SD_Pow', 'Var_AW', 'Var_Pow', 'Spread', 'Members'
#             ]
#         )
#
#     grouped_rows = []
#
#     for group_name, subdf in df.groupby(group_col):
#         member_names = sorted(subdf['AtomName'].unique().tolist())
#
#         var_aw = subdf['AW'].var(ddof=1) if len(subdf) > 1 else 0.0
#         var_pow = subdf['Pow'].var(ddof=1) if len(subdf) > 1 else 0.0
#
#         grouped_rows.append({
#             'GroupName': group_name,
#             'Count': len(subdf),
#             'Mean_AW': subdf['AW'].mean(),
#             'Mean_Pow': subdf['Pow'].mean(),
#             'SD_AW': subdf['AW'].std(ddof=1) if len(subdf) > 1 else 0.0,
#             'SD_Pow': subdf['Pow'].std(ddof=1) if len(subdf) > 1 else 0.0,
#             'Var_AW': var_aw,
#             'Var_Pow': var_pow,
#             'Spread': np.sqrt(var_aw + var_pow),
#             'Members': ', '.join(member_names),
#         })
#
#     stats_df = pd.DataFrame(grouped_rows)
#
#     stats_df['SD_AW'] = stats_df['SD_AW'].fillna(0.0)
#     stats_df['SD_Pow'] = stats_df['SD_Pow'].fillna(0.0)
#     stats_df['Spread'] = stats_df['Spread'].fillna(0.0)
#
#     stats_df = stats_df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)
#
#     return stats_df
#
#
# def filter_plot_groups(
#     stats_df: pd.DataFrame,
#     plot_min_count: int = 50,
#     max_spread: Optional[float] = None,
#     force_include: Optional[Set[str]] = None
# ) -> pd.DataFrame:
#     if force_include is None:
#         force_include = set()
#
#     mask = stats_df['Count'] >= plot_min_count
#
#     if max_spread is not None:
#         mask &= stats_df['Spread'] <= max_spread
#
#     mask |= stats_df['GroupName'].isin(force_include)
#
#     plot_stats_df = stats_df[mask].copy()
#     plot_stats_df = plot_stats_df.sort_values(['Mean_AW', 'Mean_Pow']).reset_index(drop=True)
#     plot_stats_df['PlotNumber'] = np.arange(1, len(plot_stats_df) + 1)
#
#     return plot_stats_df
#
#
# def print_name_volume_stats(stats_df: pd.DataFrame):
#     print("\nGrouped atom-name centers and standard deviations in AW vs Pow volume space:\n")
#
#     for _, row in stats_df.iterrows():
#         print(
#             f"{row['GroupName']:>14s} | "
#             f"AW={row['Mean_AW']:>6.3f} ± {row['SD_AW']:.3f} | "
#             f"Pow={row['Mean_Pow']:>6.3f} ± {row['SD_Pow']:.3f} | "
#             f"Spread={row['Spread']:.3f} | "
#             f"n={int(row['Count']):>5d} | "
#             f"Members: {row['Members']}"
#         )
#
#
# def print_group_plot_table(plot_stats_df: pd.DataFrame):
#     print("\nPlotted groups:\n")
#
#     for _, row in plot_stats_df.iterrows():
#         print(
#             f"{int(row['PlotNumber']):>2d} | "
#             f"{row['GroupName']:<18s} | "
#             f"AW={row['Mean_AW']:>6.3f}, "
#             f"Pow={row['Mean_Pow']:>6.3f} | "
#             f"Spread={row['Spread']:.3f} | "
#             f"n={int(row['Count']):>5d} | "
#             f"Members: {row['Members']}"
#         )
#
#
# def save_group_plot_table(plot_stats_df: pd.DataFrame, out_path: str):
#     cols = ['PlotNumber', 'GroupName', 'Mean_AW', 'Mean_Pow', 'Spread', 'Count', 'Members']
#     plot_stats_df[cols].to_csv(out_path, index=False)
#     print(f"\nSaved plotted group table -> {out_path}")
#
#
# def get_plot_color(name: str) -> str:
#     if name.startswith('H'):
#         return '#1f77b4'
#
#     if name.startswith('C'):
#         return '#ff7f0e'
#
#     if name.startswith('N'):
#         return '#2ca02c'
#
#     if name.startswith('O'):
#         return '#d62728'
#
#     return 'gray'
#
#
# def add_covariance_ellipse(
#     ax,
#     x,
#     y,
#     color,
#     n_std: float = 1.5,
#     face_alpha: float = 0.10,
#     edge_alpha: float = 0.75,
#     linewidth: float = 2.0,
#     zorder: float = 2
# ):
#     x = np.asarray(x, dtype=float)
#     y = np.asarray(y, dtype=float)
#
#     if len(x) < 3:
#         return
#
#     cov = np.cov(x, y)
#
#     if np.any(~np.isfinite(cov)):
#         return
#
#     vals, vecs = np.linalg.eigh(cov)
#
#     order = vals.argsort()[::-1]
#     vals = vals[order]
#     vecs = vecs[:, order]
#
#     if np.any(vals < 0):
#         return
#
#     angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
#     width = 2.0 * n_std * np.sqrt(vals[0])
#     height = 2.0 * n_std * np.sqrt(vals[1])
#
#     fill = Ellipse(
#         xy=(np.mean(x), np.mean(y)),
#         width=width,
#         height=height,
#         angle=angle,
#         facecolor=color,
#         edgecolor=color,
#         linewidth=linewidth,
#         alpha=face_alpha,
#         zorder=zorder
#     )
#     ax.add_patch(fill)
#
#     edge = Ellipse(
#         xy=(np.mean(x), np.mean(y)),
#         width=width,
#         height=height,
#         angle=angle,
#         facecolor='none',
#         edgecolor=color,
#         linewidth=linewidth,
#         alpha=edge_alpha,
#         zorder=zorder + 0.1
#     )
#
#     edge.set_path_effects([
#         pe.Stroke(linewidth=linewidth + 2.5, foreground='white'),
#         pe.Normal()
#     ])
#
#     ax.add_patch(edge)
#
#
# def add_outside_group_list(fig, plot_stats_df: pd.DataFrame):
#     x_num = 0.80
#     x_name = 0.835
#     y_start = 0.90
#     y_step = 0.024
#
#     for i, (_, row) in enumerate(plot_stats_df.iterrows()):
#         y = y_start - i * y_step
#         name = row['GroupName']
#         color = plt.cm.tab20(int(name) % 20)
#
#         if y < 0.06:
#             break
#
#         fig.text(
#             x_num,
#             y,
#             f"{int(row['PlotNumber']):>2d}",
#             color=color,
#             fontsize=10,
#             fontweight='bold',
#             ha='left',
#             va='center',
#             path_effects=[
#                 pe.Stroke(linewidth=2.5, foreground='white'),
#                 pe.Normal()
#             ]
#         )
#
#         fig.text(
#             x_name,
#             y,
#             row['GroupName'],
#             color='black',
#             fontsize=10,
#             ha='left',
#             va='center'
#         )
#
#
# def plot_ml_clusters_aw_pow(
#     df: pd.DataFrame,
#     title: str,
#     save: Optional[str] = None,
#     show: bool = True
# ):
#     fig, ax = plt.subplots(figsize=(10, 8))
#
#     cluster_ids = sorted(df['MLCluster'].dropna().unique().tolist())
#
#     for cluster_id in cluster_ids:
#         subdf = df[df['MLCluster'] == cluster_id]
#
#         label = f"Cluster {cluster_id}"
#         if cluster_id == -1:
#             label = "Noise"
#
#         ax.scatter(
#             subdf['AW'],
#             subdf['Pow'],
#             s=16,
#             alpha=0.35,
#             label=label
#         )
#
#     xmin = min(df['AW'].min(), df['Pow'].min())
#     xmax = max(df['AW'].max(), df['Pow'].max())
#     pad = 0.05 * (xmax - xmin if xmax > xmin else 1.0)
#
#     ax.plot(
#         [xmin - pad, xmax + pad],
#         [xmin - pad, xmax + pad],
#         linestyle='--',
#         linewidth=2.5,
#         color='black'
#     )
#
#     ax.set_xlim(xmin - pad, xmax + pad)
#     ax.set_ylim(xmin - pad, xmax + pad)
#
#     ax.set_xlabel('AW Volume', fontsize=22)
#     ax.set_ylabel('Pow Volume', fontsize=22)
#     ax.set_title(title, fontsize=20)
#     ax.tick_params(axis='both', which='major', labelsize=16, width=2.0, length=8)
#
#     for spine in ax.spines.values():
#         spine.set_linewidth(2)
#
#     ax.set_aspect('equal', adjustable='box')
#     # ax.legend(fontsize=10, loc='best')
#
#     if save is not None:
#         plt.savefig(save, dpi=300, bbox_inches='tight')
#
#     if show:
#         plt.show()
#
#     plt.close(fig)
#
#
# def make_settings_name(
#     mode: str,
#     ml_method: str,
#     n_clusters: Optional[int],
#     eps: Optional[float],
#     min_samples: Optional[int],
#     min_cluster_size: Optional[int],
#     numeric_cols: List[str],
#     categorical_cols: List[str],
#     use_sol_binary: bool,
#     point_alpha: float,
#     ellipse_n_std: float
# ) -> str:
#     num_tag = '-'.join(numeric_cols) if numeric_cols else 'none'
#     cat_tag = '-'.join(categorical_cols) if categorical_cols else 'none'
#
#     parts = [
#         mode,
#         ml_method,
#         f"k{n_clusters}" if n_clusters is not None else None,
#         f"eps{eps}" if eps is not None else None,
#         f"minsamp{min_samples}" if min_samples is not None else None,
#         f"minclust{min_cluster_size}" if min_cluster_size is not None else None,
#         f"num_{num_tag}",
#         f"cat_{cat_tag}",
#         f"solbin_{int(use_sol_binary)}",
#         f"alpha_{point_alpha}",
#         f"ell_{ellipse_n_std}",
#     ]
#
#     safe = "_".join([str(p) for p in parts if p is not None])
#     safe = safe.replace("'", "").replace(" ", "")
#     return safe
#
#
# def plot_name_volume_groups(
#     atom_df: pd.DataFrame,
#     plot_stats_df: pd.DataFrame,
#     title: str,
#     group_col: str = 'CanonicalName',
#     volume_range: Optional[tuple] = None,
#     show_points: bool = True,
#     show_numbers: bool = True,
#     annotate_direct_groups: bool = True,
#     point_size: float = 18,
#     point_alpha: float = 0.10,
#     number_fontsize: int = 14,
#     direct_label_fontsize: int = 11,
#     show_ellipses: bool = True,
#     ellipse_n_std: float = 1.5,
#     ellipse_min_count: int = 50,
#     ellipse_max_spread: Optional[float] = None,
#     save_png: Optional[str] = None,
#     save_svg: Optional[str] = None,
#     show: bool = True
# ):
#     fig, ax = plt.subplots(figsize=(12, 9))
#
#     print_group_plot_table(plot_stats_df)
#
#     for _, row in plot_stats_df.iterrows():
#         name = row['GroupName']
#         if group_col == 'MLCluster':
#             color = plt.cm.tab20(int(name) % 20)
#         else:
#             color = get_plot_color(str(name))
#
#         group_df = atom_df[atom_df[group_col].astype(str) == str(name)]
#         print(f"\nGroupName from plot_stats_df: {name}")
#         print(f"Matching rows in atom_df['CanonicalName'] == {name}: {len(group_df)}")
#         if len(group_df) == 0:
#             print("WARNING: empty group_df")
#             print("Sample CanonicalName values:")
#             print(atom_df['CanonicalName'].dropna().astype(str).unique()[:20])
#
#         if show_points:
#             ax.scatter(
#                 group_df['AW'],
#                 group_df['Pow'],
#                 s=point_size,
#                 alpha=point_alpha,
#                 color=color,
#                 zorder=1
#             )
#
#         can_draw_ellipse = (
#             show_ellipses and
#             len(group_df) >= ellipse_min_count and
#             (ellipse_max_spread is None or row['Spread'] <= ellipse_max_spread)
#         )
#
#         if can_draw_ellipse:
#             add_covariance_ellipse(
#                 ax=ax,
#                 x=group_df['AW'].to_numpy(),
#                 y=group_df['Pow'].to_numpy(),
#                 color=color,
#                 n_std=ellipse_n_std,
#                 face_alpha=0.20,
#                 edge_alpha=1.0,
#                 linewidth=1.5,
#                 zorder=3
#             )
#
#         if show_numbers:
#             ax.text(
#                 row['Mean_AW'],
#                 row['Mean_Pow'],
#                 str(int(row['PlotNumber'])),
#                 color=color,
#                 fontsize=number_fontsize,
#                 fontweight='bold',
#                 ha='center',
#                 va='center',
#                 zorder=6,
#                 path_effects=[
#                     pe.Stroke(linewidth=3.0, foreground='white'),
#                     pe.Normal()
#                 ]
#             )
#
#         if annotate_direct_groups and row['GroupName'] in DIRECT_LABEL_GROUPS:
#             ax.text(
#                 row['Mean_AW'] + 0.18,
#                 row['Mean_Pow'] + 0.18,
#                 row['GroupName'],
#                 fontsize=direct_label_fontsize,
#                 color='black',
#                 ha='left',
#                 va='bottom',
#                 zorder=6
#             )
#
#     if volume_range is not None:
#         vmin, vmax = volume_range
#
#         ax.set_xlim(vmin, vmax)
#         ax.set_ylim(vmin, vmax)
#
#         ax.plot(
#             [vmin, vmax],
#             [vmin, vmax],
#             linestyle='--',
#             linewidth=3.5,
#             color='black',
#             alpha=0.9,
#             zorder=0
#         )
#
#     else:
#         xmin = min(atom_df['AW'].min(), atom_df['Pow'].min())
#         xmax = max(atom_df['AW'].max(), atom_df['Pow'].max())
#         pad = 0.05 * (xmax - xmin if xmax > xmin else 1.0)
#
#         ax.set_xlim(xmin - pad, xmax + pad)
#         ax.set_ylim(xmin - pad, xmax + pad)
#
#         ax.plot(
#             [xmin - pad, xmax + pad],
#             [xmin - pad, xmax + pad],
#             linestyle='--',
#             linewidth=3.5,
#             color='black',
#             alpha=0.9,
#             zorder=0
#         )
#
#     ax.set_xlabel('AW Volume', fontsize=24)
#     ax.set_ylabel('Pow Volume', fontsize=24)
#     ax.set_title(title, fontsize=22)
#
#     ax.tick_params(axis='both', which='major', labelsize=20, width=2.5, length=10)
#
#     for spine in ax.spines.values():
#         spine.set_linewidth(2)
#
#     ax.set_aspect('equal', adjustable='box')
#
#     fig.subplots_adjust(left=0.12, right=0.77, bottom=0.12, top=0.90)
#
#     add_outside_group_list(fig, plot_stats_df)
#
#     if save_png is not None:
#         plt.savefig(save_png, dpi=300, bbox_inches='tight')
#
#     if save_svg is not None:
#         plt.savefig(save_svg, bbox_inches='tight')
#
#     if show:
#         plt.show()
#
#     plt.close(fig)
