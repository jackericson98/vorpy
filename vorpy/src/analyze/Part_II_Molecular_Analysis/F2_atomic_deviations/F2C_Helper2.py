import os
import sys
import math
import tkinter as tk
from tkinter import filedialog

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SOL_RESIDUES = {
    'SOL', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'SPC', 'SPCE', 'H2O'
}

NUCLEIC_RESIDUES = {
    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT', 'DU',
    'ADE', 'CYT', 'GUA', 'THY', 'URA'
}


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)


def choose_csv() -> str | None:
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(
        title='Select patch_atom_aw_pow_pair_compare.csv',
        filetypes=[('CSV files', '*.csv'), ('All files', '*.*')]
    )
    root.destroy()
    return path if path else None


def classify_neighbor_type(residue: str) -> str:
    residue = str(residue).strip().upper()
    if residue in SOL_RESIDUES:
        return 'SOL'
    if residue in NUCLEIC_RESIDUES:
        return 'NUCLEIC'
    return 'OTHER'



def shannon_entropy(values: pd.Series) -> float:
    vals = values.dropna().astype(str)
    if len(vals) == 0:
        return np.nan

    counts = vals.value_counts(normalize=True)
    return float(-(counts * np.log2(counts)).sum())



def safe_std(values: pd.Series) -> float:
    vals = pd.to_numeric(values, errors='coerce').dropna()
    if len(vals) < 2:
        return np.nan
    return float(vals.std(ddof=1))



def safe_mean(values: pd.Series) -> float:
    vals = pd.to_numeric(values, errors='coerce').dropna()
    if len(vals) == 0:
        return np.nan
    return float(vals.mean())



def safe_sum(values: pd.Series) -> float:
    vals = pd.to_numeric(values, errors='coerce').dropna()
    if len(vals) == 0:
        return 0.0
    return float(vals.sum())



def unique_count(values: pd.Series) -> int:
    return int(values.dropna().nunique())



def build_patch_level_features(pair_df: pd.DataFrame) -> pd.DataFrame:
    df = pair_df.copy()
    df['Neighbor Class'] = df['Neighbor Residue'].apply(classify_neighbor_type)

    numeric_cols = [
        'AW Shared Surface Area',
        'AW Shared Contact Area',
        'AW Shared Mean Curvature Mean',
        'AW Neighbor Surface Area',
        'AW Neighbor Volume',
        'AW Neighbor Number of Neighbors',
        'Pow Shared Surface Area',
        'Pow Shared Contact Area',
        'Pow Shared Mean Curvature Mean',
        'Pow Neighbor Surface Area',
        'Pow Neighbor Volume',
        'Pow Neighbor Number of Neighbors',
        'Delta Volume',
        'Delta Volume Abs',
        'Delta Shared Surface Area',
        'Delta Shared Contact Area',
        'Delta Neighbor Surface Area',
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    rows = []

    for patch_idx, sub in df.groupby('Patch Index', sort=True):
        row = {
            'Patch Index': int(patch_idx),
            'Patch Atom Name': sub['Patch Atom Name'].iloc[0] if 'Patch Atom Name' in sub.columns else None,
            'Patch Residue': sub['Patch Residue'].iloc[0] if 'Patch Residue' in sub.columns else None,
            'Patch Residue Sequence': sub['Patch Residue Sequence'].iloc[0] if 'Patch Residue Sequence' in sub.columns else np.nan,
            'Patch Chain': sub['Patch Chain'].iloc[0] if 'Patch Chain' in sub.columns else None,
            'Pair Count': int(len(sub)),
            'Unique Neighbor Count': unique_count(sub['Neighbor Index']) if 'Neighbor Index' in sub.columns else int(len(sub)),
            'Unique Neighbor Residue Count': unique_count(sub['Neighbor Residue']) if 'Neighbor Residue' in sub.columns else np.nan,
            'Unique Neighbor Element Count': unique_count(sub['Neighbor Element']) if 'Neighbor Element' in sub.columns else np.nan,
            'Neighbor Residue Entropy': shannon_entropy(sub['Neighbor Residue']) if 'Neighbor Residue' in sub.columns else np.nan,
            'Neighbor Element Entropy': shannon_entropy(sub['Neighbor Element']) if 'Neighbor Element' in sub.columns else np.nan,
            'Neighbor Class Entropy': shannon_entropy(sub['Neighbor Class']),
            'Delta Volume': safe_mean(sub['Delta Volume']) if 'Delta Volume' in sub.columns else np.nan,
            'Delta Volume Abs': safe_mean(sub['Delta Volume Abs']) if 'Delta Volume Abs' in sub.columns else np.nan,
            'Outlier By Atom Name': sub['Outlier By Atom Name'].iloc[0] if 'Outlier By Atom Name' in sub.columns else None,
        }

        # SOL / nucleic / other neighbor counts and fractions
        for cls in ['SOL', 'NUCLEIC', 'OTHER']:
            count = int((sub['Neighbor Class'] == cls).sum())
            row[f'{cls} Neighbor Count'] = count
            row[f'{cls} Neighbor Fraction'] = float(count / len(sub)) if len(sub) else np.nan

        # Shared surface and contact split by neighbor class
        for prefix in ['AW', 'Pow']:
            for cls in ['SOL', 'NUCLEIC', 'OTHER']:
                mask = sub['Neighbor Class'] == cls
                row[f'{prefix} {cls} Shared Surface Area Sum'] = safe_sum(sub.loc[mask, f'{prefix} Shared Surface Area']) if f'{prefix} Shared Surface Area' in sub.columns else np.nan
                row[f'{prefix} {cls} Shared Contact Area Sum'] = safe_sum(sub.loc[mask, f'{prefix} Shared Contact Area']) if f'{prefix} Shared Contact Area' in sub.columns else np.nan

        # Curvature features
        row['AW Curvature Mean'] = safe_mean(sub['AW Shared Mean Curvature Mean']) if 'AW Shared Mean Curvature Mean' in sub.columns else np.nan
        row['AW Curvature Abs Mean'] = safe_mean(sub['AW Shared Mean Curvature Mean'].abs()) if 'AW Shared Mean Curvature Mean' in sub.columns else np.nan
        row['AW Curvature Abs Max'] = float(sub['AW Shared Mean Curvature Mean'].abs().max()) if 'AW Shared Mean Curvature Mean' in sub.columns and sub['AW Shared Mean Curvature Mean'].notna().any() else np.nan
        row['AW Curvature Std'] = safe_std(sub['AW Shared Mean Curvature Mean']) if 'AW Shared Mean Curvature Mean' in sub.columns else np.nan

        row['Pow Curvature Mean'] = safe_mean(sub['Pow Shared Mean Curvature Mean']) if 'Pow Shared Mean Curvature Mean' in sub.columns else np.nan
        row['Pow Curvature Abs Mean'] = safe_mean(sub['Pow Shared Mean Curvature Mean'].abs()) if 'Pow Shared Mean Curvature Mean' in sub.columns else np.nan
        row['Pow Curvature Abs Max'] = float(sub['Pow Shared Mean Curvature Mean'].abs().max()) if 'Pow Shared Mean Curvature Mean' in sub.columns and sub['Pow Shared Mean Curvature Mean'].notna().any() else np.nan
        row['Pow Curvature Std'] = safe_std(sub['Pow Shared Mean Curvature Mean']) if 'Pow Shared Mean Curvature Mean' in sub.columns else np.nan

        # Heterogeneity features
        for prefix in ['AW', 'Pow']:
            row[f'{prefix} Neighbor Surface Area Mean'] = safe_mean(sub[f'{prefix} Neighbor Surface Area']) if f'{prefix} Neighbor Surface Area' in sub.columns else np.nan
            row[f'{prefix} Neighbor Surface Area Std'] = safe_std(sub[f'{prefix} Neighbor Surface Area']) if f'{prefix} Neighbor Surface Area' in sub.columns else np.nan
            row[f'{prefix} Neighbor Volume Mean'] = safe_mean(sub[f'{prefix} Neighbor Volume']) if f'{prefix} Neighbor Volume' in sub.columns else np.nan
            row[f'{prefix} Neighbor Volume Std'] = safe_std(sub[f'{prefix} Neighbor Volume']) if f'{prefix} Neighbor Volume' in sub.columns else np.nan
            row[f'{prefix} Neighbor Degree Mean'] = safe_mean(sub[f'{prefix} Neighbor Number of Neighbors']) if f'{prefix} Neighbor Number of Neighbors' in sub.columns else np.nan
            row[f'{prefix} Neighbor Degree Std'] = safe_std(sub[f'{prefix} Neighbor Number of Neighbors']) if f'{prefix} Neighbor Number of Neighbors' in sub.columns else np.nan
            row[f'{prefix} Shared Surface Area Sum'] = safe_sum(sub[f'{prefix} Shared Surface Area']) if f'{prefix} Shared Surface Area' in sub.columns else np.nan
            row[f'{prefix} Shared Surface Area Mean'] = safe_mean(sub[f'{prefix} Shared Surface Area']) if f'{prefix} Shared Surface Area' in sub.columns else np.nan
            row[f'{prefix} Shared Surface Area Std'] = safe_std(sub[f'{prefix} Shared Surface Area']) if f'{prefix} Shared Surface Area' in sub.columns else np.nan
            row[f'{prefix} Shared Contact Area Sum'] = safe_sum(sub[f'{prefix} Shared Contact Area']) if f'{prefix} Shared Contact Area' in sub.columns else np.nan
            row[f'{prefix} Shared Contact Area Mean'] = safe_mean(sub[f'{prefix} Shared Contact Area']) if f'{prefix} Shared Contact Area' in sub.columns else np.nan
            row[f'{prefix} Shared Contact Area Std'] = safe_std(sub[f'{prefix} Shared Contact Area']) if f'{prefix} Shared Contact Area' in sub.columns else np.nan

        # Delta geometry features
        row['Delta Shared Surface Area Sum'] = safe_sum(sub['Delta Shared Surface Area']) if 'Delta Shared Surface Area' in sub.columns else np.nan
        row['Delta Shared Surface Area Mean'] = safe_mean(sub['Delta Shared Surface Area']) if 'Delta Shared Surface Area' in sub.columns else np.nan
        row['Delta Shared Surface Area Abs Sum'] = safe_sum(sub['Delta Shared Surface Area'].abs()) if 'Delta Shared Surface Area' in sub.columns else np.nan
        row['Delta Shared Surface Area Std'] = safe_std(sub['Delta Shared Surface Area']) if 'Delta Shared Surface Area' in sub.columns else np.nan

        row['Delta Shared Contact Area Sum'] = safe_sum(sub['Delta Shared Contact Area']) if 'Delta Shared Contact Area' in sub.columns else np.nan
        row['Delta Shared Contact Area Mean'] = safe_mean(sub['Delta Shared Contact Area']) if 'Delta Shared Contact Area' in sub.columns else np.nan
        row['Delta Shared Contact Area Abs Sum'] = safe_sum(sub['Delta Shared Contact Area'].abs()) if 'Delta Shared Contact Area' in sub.columns else np.nan
        row['Delta Shared Contact Area Std'] = safe_std(sub['Delta Shared Contact Area']) if 'Delta Shared Contact Area' in sub.columns else np.nan

        row['Delta Neighbor Surface Area Sum'] = safe_sum(sub['Delta Neighbor Surface Area']) if 'Delta Neighbor Surface Area' in sub.columns else np.nan
        row['Delta Neighbor Surface Area Mean'] = safe_mean(sub['Delta Neighbor Surface Area']) if 'Delta Neighbor Surface Area' in sub.columns else np.nan
        row['Delta Neighbor Surface Area Abs Sum'] = safe_sum(sub['Delta Neighbor Surface Area'].abs()) if 'Delta Neighbor Surface Area' in sub.columns else np.nan
        row['Delta Neighbor Surface Area Std'] = safe_std(sub['Delta Neighbor Surface Area']) if 'Delta Neighbor Surface Area' in sub.columns else np.nan

        rows.append(row)

    patch_df = pd.DataFrame(rows)

    # Higher-level composite features
    if 'AW Neighbor Surface Area Std' in patch_df.columns and 'AW Neighbor Volume Std' in patch_df.columns:
        patch_df['AW Heterogeneity Score'] = patch_df[['AW Neighbor Surface Area Std', 'AW Neighbor Volume Std', 'AW Neighbor Degree Std']].mean(axis=1, skipna=True)
        patch_df['Pow Heterogeneity Score'] = patch_df[['Pow Neighbor Surface Area Std', 'Pow Neighbor Volume Std', 'Pow Neighbor Degree Std']].mean(axis=1, skipna=True)

    if 'AW SOL Shared Surface Area Sum' in patch_df.columns and 'AW Shared Surface Area Sum' in patch_df.columns:
        denom = patch_df['AW Shared Surface Area Sum'].replace(0, np.nan)
        patch_df['AW SOL Shared Surface Fraction'] = patch_df['AW SOL Shared Surface Area Sum'] / denom

    if 'Pow SOL Shared Surface Area Sum' in patch_df.columns and 'Pow Shared Surface Area Sum' in patch_df.columns:
        denom = patch_df['Pow Shared Surface Area Sum'].replace(0, np.nan)
        patch_df['Pow SOL Shared Surface Fraction'] = patch_df['Pow SOL Shared Surface Area Sum'] / denom

    patch_df['Curvature_vs_Heterogeneity'] = patch_df[['AW Curvature Abs Mean', 'AW Heterogeneity Score']].mean(axis=1, skipna=True)

    return patch_df



def safe_corr(x: pd.Series, y: pd.Series, method: str = 'pearson') -> float:
    pair = pd.concat([x, y], axis=1).dropna()
    if len(pair) < 3:
        return np.nan

    x_clean = pd.to_numeric(pair.iloc[:, 0], errors='coerce')
    y_clean = pd.to_numeric(pair.iloc[:, 1], errors='coerce')
    valid = pd.concat([x_clean, y_clean], axis=1).dropna()

    if len(valid) < 3:
        return np.nan

    x_clean = valid.iloc[:, 0]
    y_clean = valid.iloc[:, 1]

    if x_clean.nunique() < 2 or y_clean.nunique() < 2:
        return np.nan

    return float(x_clean.corr(y_clean, method=method))



def build_correlation_table(patch_df: pd.DataFrame, target_col: str = 'Delta Volume Abs') -> pd.DataFrame:
    exclude_cols = {
        'Patch Index', 'Patch Atom Name', 'Patch Residue', 'Patch Residue Sequence',
        'Patch Chain', 'Outlier By Atom Name'
    }

    feature_cols = [
        col for col in patch_df.columns
        if col not in exclude_cols and col != target_col
    ]

    rows = []
    for feature in feature_cols:
        pearson = safe_corr(patch_df[target_col], patch_df[feature], method='pearson')
        spearman = safe_corr(patch_df[target_col], patch_df[feature], method='spearman')
        n = int(pd.concat([patch_df[target_col], patch_df[feature]], axis=1).dropna().shape[0])

        rows.append({
            'Feature': feature,
            'N': n,
            'Pearson': pearson,
            'Spearman': spearman,
            'Abs Pearson': abs(pearson) if pd.notna(pearson) else np.nan,
            'Abs Spearman': abs(spearman) if pd.notna(spearman) else np.nan,
        })

    corr_df = pd.DataFrame(rows)
    corr_df = corr_df.sort_values(['Abs Pearson', 'Abs Spearman', 'Feature'], ascending=[False, False, True]).reset_index(drop=True)
    return corr_df



def build_heatmap_matrix(patch_df: pd.DataFrame, selected_features: list[str], target_col: str = 'Delta Volume Abs') -> pd.DataFrame:
    cols = [target_col] + [c for c in selected_features if c in patch_df.columns]
    data = patch_df[cols].apply(pd.to_numeric, errors='coerce')
    return data.corr(method='pearson')



def plot_heatmap(corr_matrix: pd.DataFrame, out_path: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))

    matrix = corr_matrix.values
    im = ax.imshow(matrix, aspect='auto')

    ax.set_xticks(range(len(corr_matrix.columns)))
    ax.set_yticks(range(len(corr_matrix.index)))
    ax.set_xticklabels(corr_matrix.columns, rotation=90, fontsize=10)
    ax.set_yticklabels(corr_matrix.index, fontsize=10)
    ax.set_title(title, fontsize=14)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            if pd.notna(val):
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=8)

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Pearson r', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)



def plot_target_scatter_grid(patch_df: pd.DataFrame, target_col: str, features: list[str], out_path: str) -> None:
    features = [f for f in features if f in patch_df.columns]
    if not features:
        return

    n = len(features)
    ncols = 2
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4.5 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, feature in zip(axes, features):
        sub = patch_df[[feature, target_col]].apply(pd.to_numeric, errors='coerce').dropna()
        ax.scatter(sub[feature], sub[target_col], alpha=0.7, s=35)
        ax.set_xlabel(feature)
        ax.set_ylabel(target_col)
        ax.set_title(feature)

        if len(sub) >= 2 and sub[feature].nunique() >= 2:
            coeffs = np.polyfit(sub[feature], sub[target_col], deg=1)
            xs = np.linspace(sub[feature].min(), sub[feature].max(), 200)
            ys = coeffs[0] * xs + coeffs[1]
            ax.plot(xs, ys, linewidth=2)

    for ax in axes[n:]:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)



def main() -> None:
    csv_path = choose_csv()
    if not csv_path:
        print('No CSV selected.')
        return

    folder = os.path.dirname(csv_path)
    pair_df = pd.read_csv(csv_path)

    patch_df = build_patch_level_features(pair_df)
    corr_df = build_correlation_table(patch_df, target_col='Delta Volume Abs')

    # Pick a focused set for the heatmap and scatter grid
    preferred_features = [
        'AW Curvature Abs Mean',
        'AW Curvature Abs Max',
        'AW Heterogeneity Score',
        'Pow Heterogeneity Score',
        'Neighbor Class Entropy',
        'Neighbor Residue Entropy',
        'SOL Neighbor Fraction',
        'AW SOL Shared Surface Fraction',
        'Pow SOL Shared Surface Fraction',
        'AW SOL Shared Surface Area Sum',
        'Pow SOL Shared Surface Area Sum',
        'AW Shared Surface Area Sum',
        'Pow Shared Surface Area Sum',
        'Delta Shared Surface Area Abs Sum',
        'Delta SharedContactDummy',
        'Delta Shared Contact Area Abs Sum',
        'Delta Neighbor Surface Area Abs Sum',
        'Unique Neighbor Count',
    ]

    # normalize the preferred list in case one dummy slipped in
    preferred_features = [f for f in preferred_features if f in patch_df.columns]

    heatmap_matrix = build_heatmap_matrix(
        patch_df,
        selected_features=preferred_features,
        target_col='Delta Volume Abs'
    )

    patch_csv = os.path.join(folder, 'patch_level_features.csv')
    corr_csv = os.path.join(folder, 'patch_level_feature_correlations.csv')
    heatmap_png = os.path.join(folder, 'patch_level_feature_heatmap.png')
    scatter_png = os.path.join(folder, 'patch_level_feature_scatter_grid.png')
    summary_txt = os.path.join(folder, 'patch_level_feature_summary.txt')

    patch_df.to_csv(patch_csv, index=False)
    corr_df.to_csv(corr_csv, index=False)
    plot_heatmap(heatmap_matrix, heatmap_png, 'Patch-level Feature Correlation Heatmap')

    top_scatter_features = corr_df['Feature'].dropna().head(6).tolist()
    plot_target_scatter_grid(patch_df, 'Delta Volume Abs', top_scatter_features, scatter_png)

    with open(summary_txt, 'w', encoding='utf-8') as handle:
        handle.write('PATCH-LEVEL FEATURE SUMMARY\n')
        handle.write('=' * 80 + '\n\n')

        handle.write('Top correlations with Delta Volume Abs:\n')
        handle.write(corr_df.head(15).to_string(index=False))
        handle.write('\n\n')

        if 'SOL Neighbor Fraction' in patch_df.columns:
            nonzero_sol = int((patch_df['SOL Neighbor Fraction'].fillna(0) > 0).sum())
            handle.write(f'Patch atoms with non-zero SOL Neighbor Fraction: {nonzero_sol} / {len(patch_df)}\n')
        else:
            handle.write('SOL Neighbor Fraction was not present in the patch-level table.\n')

        unique_neighbor_classes = sorted(pair_df['Neighbor Residue'].dropna().astype(str).str.upper().unique().tolist()) if 'Neighbor Residue' in pair_df.columns else []
        sol_residues_present = [r for r in unique_neighbor_classes if r in SOL_RESIDUES]
        handle.write(f'SOL-like residue names present in pair CSV: {sol_residues_present}\n')

    print('\n=== WRITTEN FILES ===')
    print(patch_csv)
    print(corr_csv)
    print(heatmap_png)
    print(scatter_png)
    print(summary_txt)

    print('\n=== TOP CORRELATIONS ===')
    print(corr_df.head(15).to_string(index=False))

    if 'SOL Neighbor Fraction' in patch_df.columns:
        nonzero_sol = int((patch_df['SOL Neighbor Fraction'].fillna(0) > 0).sum())
        print(f'\nPatch atoms with non-zero SOL Neighbor Fraction: {nonzero_sol} / {len(patch_df)}')

    if 'Neighbor Residue' in pair_df.columns:
        unique_neighbor_classes = sorted(pair_df['Neighbor Residue'].dropna().astype(str).str.upper().unique().tolist())
        sol_residues_present = [r for r in unique_neighbor_classes if r in SOL_RESIDUES]
        print(f'SOL-like residue names present in pair CSV: {sol_residues_present}')


if __name__ == '__main__':
    main()
