import os
import sys
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from typing import List, Optional, Dict


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def _read_scheme_logs(folder: str):
    """
    Reads aw/pow/prm logs from either:
      folder/{aw_logs.csv,pow_logs.csv,prm_logs.csv}
    or
      folder/{aw/aw_logs.csv, pow/pow_logs.csv, prm/prm_logs.csv}
    """
    try:
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)

    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)

    return aw_logs, pow_logs, prm_logs


def select_folders_multi(title: str = "Select a folder (Cancel to finish)"):
    """
    Tkinter doesn’t provide a native multi-directory picker, so this prompts repeatedly.
    Select folders one at a time; press Cancel when you’re done.
    """
    root = tk.Tk()
    root.withdraw()

    folders = []
    while True:
        folder = filedialog.askdirectory(title=title)
        if not folder:
            break
        folders.append(folder)

    root.destroy()
    return folders


def rotate_to_diagonal(xs: np.ndarray, ys: np.ndarray):
    """
    Rotate (x,y) so that:
      u = along y=x
      v = perpendicular to y=x (deviation / bias)
    """
    u = (xs + ys) / np.sqrt(2.0)
    v = (ys - xs) / np.sqrt(2.0)

    return u, v


def two_pass_dbscan_diagonal(
    df: pd.DataFrame,
    eps1: float,
    min_samples1: int,
    eps2: float,
    min_samples2: Optional[int] = None
):
    """
    Pass 1: DBSCAN(eps1) on all points.
    Pass 2: DBSCAN(eps2) only on points labeled -1 from pass 1.
    Merges labels so pass-2 clusters get new ids appended after pass-1 ids.

    Returns:
      final_labels (np.ndarray), labels1 (np.ndarray), labels2 (np.ndarray on outliers only, else -999)
    """
    if min_samples2 is None:
        min_samples2 = min_samples1

    # ---- Build diagonal coords once ----
    xs = df['AW'].to_numpy(dtype=float)
    ys = df['Pow'].to_numpy(dtype=float)
    u, v = rotate_to_diagonal(xs, ys)

    Z = StandardScaler().fit_transform(np.column_stack([u, v]))

    # ---- Pass 1 ----
    labels1 = DBSCAN(eps=eps1, min_samples=min_samples1).fit_predict(Z)
    final_labels = labels1.copy()

    out_mask = labels1 == -1
    if not np.any(out_mask):
        labels2_full = np.full_like(final_labels, fill_value=-999)
        return final_labels, labels1, labels2_full, u, v

    # ---- Pass 2 (only outliers) ----
    Z2 = Z[out_mask]
    labels2 = DBSCAN(eps=eps2, min_samples=min_samples2).fit_predict(Z2)

    # Remap pass-2 cluster ids so they don’t collide with pass-1 ids
    existing = sorted(set(labels1.tolist()))
    max_id = max([lab for lab in existing if lab != -1], default=-1)
    next_id = max_id + 1

    labels2_mapped = labels2.copy()
    for lab in sorted(set(labels2.tolist())):
        if lab == -1:
            continue
        labels2_mapped[labels2 == lab] = next_id
        next_id += 1

    # Merge: replace only where pass2 formed a cluster
    out_idx = np.where(out_mask)[0]
    for k, i in enumerate(out_idx):
        if labels2_mapped[k] != -1:
            final_labels[i] = labels2_mapped[k]  # promote to new cluster id

    labels2_full = np.full_like(final_labels, fill_value=-999)
    labels2_full[out_mask] = labels2  # store original pass2 labels for debugging

    return final_labels, labels1, labels2_full, u, v


def cluster_dbscan_diagonal(df: pd.DataFrame, eps=0.35, min_samples=10):
    """
    DBSCAN in diagonal coordinates (u along y=x, v perpendicular).
    Expects df columns: ['AW', 'Pow'].
    Returns: labels (np.ndarray), u (np.ndarray), v (np.ndarray)
    """
    xs = df['AW'].to_numpy(dtype=float)
    ys = df['Pow'].to_numpy(dtype=float)

    u, v = rotate_to_diagonal(xs, ys)

    Z = np.column_stack([u, v])
    Z = StandardScaler().fit_transform(Z)

    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(Z)

    return labels, u, v


def report_outliers(df: pd.DataFrame, out_path: Optional[str] = None, top_names: int = 15):
    out_df = df[df['cluster'] == -1].copy()

    n_total = len(df)
    n_out = len(out_df)
    frac = 0.0 if n_total == 0 else n_out / n_total

    print(f"\nOutliers: n={n_out} / {n_total} ({frac:.2%})")

    if n_out == 0:
        return

    # What are they?
    print("\nTop outlier AtomName counts:")
    print(out_df['AtomName'].value_counts().head(top_names).to_string())

    # Print a small sample for manual sanity checking
    cols = [c for c in ['Folder', 'Index', 'AtomName', 'AW', 'Pow', 'u', 'v'] if c in out_df.columns]
    print("\nSample outliers (first 30):")
    print(out_df[cols].head(30).to_string(index=False))

    if out_path is not None:
        out_df.to_csv(out_path, index=False)
        print(f"\nSaved outliers CSV -> {out_path}")


def print_cluster_name_audit(
    df: pd.DataFrame,
    name_col: str = 'AtomName',
    top_k: int = 8,
    show_names: Optional[List[str]] = None
):
    """
    Prints:
      1) Top-k name fractions per cluster
      2) Optional: distribution of specific names across clusters (e.g., ['CA'])
    Expects df columns: ['cluster', name_col]
    """
    df_in = df[df['cluster'] != -1].copy()
    if len(df_in) == 0:
        print("\nNo non-outlier clusters to audit (everything is -1).")
        return

    comp = (
        df_in
        .groupby(['cluster', name_col])
        .size()
        .reset_index(name='n')
    )

    comp['frac'] = comp['n'] / comp.groupby('cluster')['n'].transform('sum')

    top = (
        comp
        .sort_values(['cluster', 'frac'], ascending=[True, False])
        .groupby('cluster')
        .head(top_k)
    )

    print("\nTop atom-name fractions per cluster:")
    print(top.to_string(index=False))

    if show_names is not None:
        for nm in show_names:
            mask = df[name_col] == nm
            if not mask.any():
                print(f"\nName '{nm}' not present in data.")
                continue

            print(f"\n'{nm}' atoms per cluster:")
            print(df.loc[mask, 'cluster'].value_counts(dropna=False).sort_index())


def plot_clusters_with_name_overlay(
    df: pd.DataFrame,
    vmin: float,
    vmax: float,
    title: str,
    overlay_names: Optional[List[str]] = None,
    name_col: str = 'AtomName',
    show: bool = True,
    save: Optional[str] = None
):
    """
    Plots (AW, Pow) colored by cluster.
    Optionally overlays specific AtomName groups with hollow circles (e.g., ['CA']).
    Uses manual layout (subplots_adjust) to avoid constrained_layout/tight_layout warnings.
    """
    fig, ax = plt.subplots(figsize=(10, 8))  # <-- NO constrained_layout

    ax.plot([vmin, vmax], [vmin, vmax], color='black', linestyle='--', linewidth=3, alpha=0.7)

    xs = df['AW'].to_numpy(dtype=float)
    ys = df['Pow'].to_numpy(dtype=float)
    clabs = df['cluster'].to_numpy()

    unique_labels = sorted(set(clabs.tolist()))

    for lab in unique_labels:
        mask = clabs == lab

        if lab == -1:
            ax.scatter(xs[mask], ys[mask], s=30, alpha=0.35, label='Outliers', marker='x')
        else:
            ax.scatter(xs[mask], ys[mask], s=40, alpha=0.6, label=f'Cluster {lab}')

    if overlay_names is not None:
        for nm in overlay_names:
            mask_nm = df[name_col] == nm
            if not mask_nm.any():
                continue

            ax.scatter(
                df.loc[mask_nm, 'AW'],
                df.loc[mask_nm, 'Pow'],
                s=140,
                facecolors='none',
                linewidths=2.5,
                label=f"{nm} overlay"
            )

    ax.set_xlabel('AW Volume', fontsize=25)
    ax.set_ylabel('Pow Volume', fontsize=25)
    ax.set_title(title, fontsize=22)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_xticks([5, 10, 15, 20])
    ax.set_yticks([5, 10, 15, 20])
    ax.tick_params(axis='both', which='major', labelsize=25, width=3, length=12)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # --- Make legend manageable: show outliers + overlays + top-N clusters by size ---
    cluster_sizes = (
        df[df['cluster'] != -1]
        .groupby('cluster')
        .size()
        .sort_values(ascending=False)
    )

    report_outliers(df, out_path=None)  # or provide a path

    top_n = 8
    top_clusters = set(cluster_sizes.head(top_n).index.tolist())

    handles, texts = ax.get_legend_handles_labels()
    new_handles = []
    new_texts = []

    for h, t in zip(handles, texts):
        if t == 'Outliers':
            new_handles.append(h)
            new_texts.append(t)
            continue

        if t.endswith('overlay'):
            new_handles.append(h)
            new_texts.append(t)
            continue

        if t.startswith('Cluster '):
            try:
                lab = int(t.replace('Cluster ', '').strip())
            except ValueError:
                continue
            if lab in top_clusters:
                new_handles.append(h)
                new_texts.append(t)

    ax.legend(
        new_handles,
        new_texts,
        fontsize=12,
        frameon=True,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0
    )

    # Manual spacing (now valid because no layout engine is active)
    fig.subplots_adjust(left=0.14, right=0.75, bottom=0.14, top=0.88)

    if save is not None:
        plt.savefig(save, dpi=300)

    if show:
        plt.show()

    plt.close(fig)


def plot_vols(
    volume_range=(3, 22),
    cluster: bool = True,
    cluster_eps: float = 0.35,
    cluster_min_samples: int = 10,
    audit_names: bool = True,
    overlay_names: Optional[List[str]] = None,
    atom_name_field: str = 'Name',
    cluster_multiplier: float = 2.0
):
    """
    Pooled AW vs Pow volume plot across multiple selected folders.
    If cluster=True: DBSCAN in diagonal coords + name audit + optional overlays.

    overlay_names example: ['CA', 'CB', 'O', 'N']
    atom_name_field: column in your logs' atoms df to use as the "name convention" (usually 'Name')
    """
    folders = select_folders_multi()
    if len(folders) == 0:
        return

    vmin, vmax = volume_range

    records: List[Dict[str, object]] = []

    for folder in folders:
        aw_logs, pow_logs, prm_logs = _read_scheme_logs(folder)

        for _, atom in aw_logs['atoms'].iterrows():
            pow_match = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == atom['Index']]
            prm_match = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == atom['Index']]

            if len(pow_match) == 0 or len(prm_match) == 0:
                continue

            pow_atom = pow_match.to_dict(orient='records')[0]
            prm_atom = prm_match.to_dict(orient='records')[0]

            aw_v = float(atom['Volume'])
            pow_v = float(pow_atom['Volume'])
            prm_v = float(prm_atom['Volume'])

            if (
                aw_v < vmin or aw_v > vmax or
                pow_v < vmin or pow_v > vmax or
                prm_v < vmin or prm_v > vmax
            ):
                continue

            raw_name = atom.get(atom_name_field, '')
            atom_name = str(raw_name).strip().upper()

            records.append({
                'Folder': folder,
                'Index': int(atom['Index']),
                'AtomName': atom_name,
                'AW': aw_v,
                'Pow': pow_v
            })

    if len(records) == 0:
        print("No atoms passed filters / matching across selected folders.")
        return

    df = pd.DataFrame(records)

    title_base = f"Pooled Volume Comparison (n_folders={len(folders)}, n_atoms={len(df)})"

    if cluster:
        final_labels, labels1, labels2_full, u, v = two_pass_dbscan_diagonal(
            df=df,
            eps1=cluster_eps,  # your strict eps
            min_samples1=cluster_min_samples,
            eps2=cluster_eps * cluster_multiplier,  # catch pass (tune this)
            min_samples2=max(6, cluster_min_samples // 2)
        )

        df['cluster'] = final_labels
        df['u'] = u
        df['v'] = v

        print("\nTwo-pass DBSCAN summary:")
        n_total = len(df)
        n_out1 = int(np.sum(labels1 == -1))
        n_out_final = int(np.sum(final_labels == -1))
        print(f"  pass1 outliers: {n_out1}/{n_total} ({n_out1 / n_total:.2%})")
        print(f"  final outliers: {n_out_final}/{n_total} ({n_out_final / n_total:.2%})")

        # Use FINAL labels (after both passes)
        labels = df['cluster'].to_numpy()

        cluster_ids = set(labels.tolist())

        n_clusters = len(cluster_ids - {-1})
        n_outliers = int(np.sum(labels == -1))

        title = (
            f"{title_base} | "
            f"DBSCAN clusters={n_clusters}, outliers={n_outliers}"
        )

        if audit_names:
            print_cluster_name_audit(
                df=df,
                name_col='AtomName',
                top_k=8,
                show_names=overlay_names
            )

        plot_clusters_with_name_overlay(
            df=df,
            vmin=vmin,
            vmax=vmax,
            title=title,
            overlay_names=overlay_names,
            name_col='AtomName',
            show=True
        )

        return

    print("cluster=False path not implemented in this simplified audit-focused version.")


if __name__ == "__main__":
    # Example: audit whether CA atoms concentrate in one cluster
    plot_vols(
        volume_range=(3, 22),
        cluster=True,
        cluster_eps=0.04,
        cluster_multiplier=2,
        cluster_min_samples=50,
        audit_names=True,
        overlay_names=None,
        atom_name_field='Name'
    )
