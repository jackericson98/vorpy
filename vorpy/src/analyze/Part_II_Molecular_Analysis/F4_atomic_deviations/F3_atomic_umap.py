"""Joint AW/Power atom-descriptor UMAP and exploratory HDBSCAN clustering.

All usable scalar measurements and inertia eigenvalues enter one standardized
feature space. Atom identity, coordinates, cell-status flags and neighbor IDs
remain metadata. Clustering uses a five-dimensional UMAP; an independent 2-D
UMAP of the same features supplies the figure coordinates.
"""
from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
from vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale.F3_common import (
    discover_systems, read_pair,
)


IDENTITY = ['Index', 'Chain', 'Residue', 'Residue Sequence', 'Atom']
EXCLUDED = {
    'Index', 'Name', 'Atom', 'Atom Name', 'Chain', 'Residue', 'Residue Sequence',
    'X', 'Y', 'Z', 'Complete Cell?', 'Inner Ball?', 'Closest Neighbor', 'neighbors',
    'Solvent Layer', 'Center of Mass', 'Moment of Inertia Tensor', 'Bounding Box',
}


def prepare_atoms(atoms, scheme):
    """Preserve compound atom identity and select physical numeric measurements."""
    name = next((column for column in ('Name', 'Atom', 'Atom Name') if column in atoms), None)
    if 'Index' not in atoms or name is None:
        raise ValueError('Atom logs require Index and atom-name columns.')
    identity = pd.DataFrame(index=atoms.index)
    identity['Index'] = pd.to_numeric(atoms['Index'], errors='raise')
    if not np.isfinite(identity['Index']).all() or (identity['Index'] % 1 != 0).any():
        raise ValueError('Atom indices must be finite integers.')
    for column, default in [('Chain', ''), ('Residue', ''), ('Residue Sequence', -1)]:
        identity[column] = atoms[column].fillna(default) if column in atoms else default
    identity['Residue Sequence'] = pd.to_numeric(identity['Residue Sequence'], errors='raise')
    for column in ('Chain', 'Residue'):
        identity[column] = identity[column].astype(str).str.strip().str.upper()
    identity['Atom'] = atoms[name].astype(str).str.strip().str.upper()
    key = pd.MultiIndex.from_frame(identity[IDENTITY])
    if key.has_duplicates:
        raise ValueError(f'Duplicate {scheme} atom identities prevent one-to-one matching.')
    complete = pd.Series(True, index=atoms.index)
    if 'Complete Cell?' in atoms:
        complete = atoms['Complete Cell?'].astype(str).str.strip().str.lower().isin(['true', '1', '1.0'])
    output = pd.DataFrame(index=key)
    excluded = {}
    for column in atoms:
        if column in EXCLUDED or pd.api.types.is_bool_dtype(atoms[column]):
            excluded[column] = 'identity, status, coordinates, or structured metadata'
            continue
        values = pd.to_numeric(atoms[column], errors='coerce')
        if values.notna().any() or atoms[column].isna().all():
            output[f'{scheme} {column}'] = values.to_numpy()
        else:
            excluded[column] = 'not a scalar numeric measurement'
    if 'Moment of Inertia Tensor' in atoms:
        eigenvalues = []
        for value in atoms['Moment of Inertia Tensor']:
            try:
                tensor = np.asarray(value, dtype=float)
                if (tensor.shape != (3, 3) or not np.isfinite(tensor).all()
                        or not np.array_equal(tensor, tensor.T)):
                    raise ValueError('Invalid inertia tensor')
                eigenvalues.append(np.linalg.eigvalsh(tensor))
            except (TypeError, ValueError, np.linalg.LinAlgError):
                eigenvalues.append([np.nan] * 3)
        for axis in range(3):
            output[f'{scheme} Inertia Eigenvalue {axis + 1}'] = np.asarray(eigenvalues).reshape(-1, 3)[:, axis]
        excluded['Moment of Inertia Tensor'] = 'represented by rotation-invariant eigenvalues'
    return output.loc[complete.to_numpy()], {
        'source_atoms': len(atoms), 'incomplete_atoms_excluded': int((~complete).sum()),
        'completeness_available': 'Complete Cell?' in atoms, 'excluded_columns': excluded,
    }


def collect_descriptors(data_root, exclude_keys):
    frames, reports = [], []
    for paths in discover_systems(str(data_root), exclude_keys=exclude_keys):
        aw_logs, power_logs = read_pair(paths, need_surfs=False)
        aw, aw_report = prepare_atoms(aw_logs['atoms'], 'AW')
        power, power_report = prepare_atoms(power_logs['atoms'], 'Power')
        joined = aw.join(power, how='inner', validate='one_to_one').reset_index()
        joined.insert(0, 'System', paths.name)
        frames.append(joined)
        reports.append({'system': paths.name, 'matched_atoms': len(joined),
                        'aw_unmatched': len(aw) - len(joined), 'power_unmatched': len(power) - len(joined),
                        'AW': aw_report, 'Power': power_report})
        print(f'{paths.name}: {len(joined):,} matched complete atoms', flush=True)
    if not frames or not any(len(frame) for frame in frames):
        raise ValueError('No matched complete AW/Power atoms were found.')
    return pd.concat(frames, ignore_index=True), reports


def prepare_features(atoms, missing='median'):
    columns = [c for c in atoms if c.startswith(('AW ', 'Power '))]
    raw = atoms[columns].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)
    missing_counts = raw.isna().sum().to_dict()
    all_missing = raw.columns[raw.isna().all()].tolist()
    raw = raw.drop(columns=all_missing)
    keep = raw.notna().all(axis=1) if missing == 'complete' else raw.notna().any(axis=1)
    raw = raw.loc[keep]
    if raw.empty:
        raise ValueError('No usable atom measurements remain after missing-value filtering.')
    values = pd.DataFrame(SimpleImputer(strategy='median').fit_transform(raw),
                          index=raw.index, columns=raw.columns)
    constant = values.columns[values.nunique() <= 1].tolist()
    values = values.drop(columns=constant)
    duplicates, seen = {}, []
    for column in values:
        previous = next((other for other in seen if values[column].equals(values[other])), None)
        if previous is None:
            seen.append(column)
        else:
            duplicates[column] = previous
    values = values.drop(columns=list(duplicates))
    if values.shape[1] < 2 or len(values) < 4:
        raise ValueError('UMAP requires at least four usable atoms and two varying features.')
    scaler = StandardScaler()
    scaled = scaler.fit_transform(values)
    report = {'features': list(values.columns), 'missing_policy': missing,
              'missing_values': {k: int(v) for k, v in missing_counts.items()},
              'all_missing_features': all_missing, 'constant_features': constant,
              'duplicate_features': duplicates, 'excluded_rows': int((~keep).sum()),
              'scaler_mean': scaler.mean_.tolist(), 'scaler_scale': scaler.scale_.tolist()}
    return atoms.loc[keep].reset_index(drop=True), scaled, report


def fit_clusters(scaled, *, neighbors=30, min_cluster_size=15, seed=42):
    if not 2 <= min_cluster_size <= len(scaled):
        raise ValueError('min_cluster_size must be between 2 and the usable atom count.')
    if neighbors < 2:
        raise ValueError('neighbors must be at least 2.')
    try:
        from umap import UMAP
    except ImportError as error:
        raise RuntimeError('This analysis requires umap-learn; install the optional analysis dependencies.') from error
    common = dict(n_neighbors=min(neighbors, len(scaled) - 1), metric='euclidean',
                  random_state=seed, n_jobs=1, init='random')
    print('Fitting 5-D UMAP for clustering...', flush=True)
    embedding = UMAP(n_components=5, min_dist=0.0, **common).fit_transform(scaled)
    if not np.isfinite(embedding).all():
        raise ValueError('UMAP produced non-finite clustering coordinates; check disconnected observations.')
    clustering = HDBSCAN(min_cluster_size=min_cluster_size).fit(embedding)
    print('Fitting 2-D UMAP for visualization...', flush=True)
    visual = UMAP(n_components=2, min_dist=0.1, **common).fit_transform(scaled)
    if not np.isfinite(visual).all():
        raise ValueError('UMAP produced non-finite visualization coordinates; check disconnected observations.')
    return embedding, visual, clustering.labels_, clustering.probabilities_


def export_results(atoms, scaled, feature_report, collection, output, *, neighbors=30,
                   min_cluster_size=15, seed=42):
    embedding, visual, labels, probabilities = fit_clusters(
        scaled, neighbors=neighbors, min_cluster_size=min_cluster_size, seed=seed)
    results = atoms.copy()
    results[['UMAP1', 'UMAP2']] = visual
    for index in range(5):
        results[f'Cluster UMAP{index + 1}'] = embedding[:, index]
    results['Cluster'] = labels
    results['Cluster Probability'] = probabilities
    output.mkdir(parents=True, exist_ok=True)
    results.to_csv(output / 'F3_atomic_umap_assignments.csv', index=False)
    pd.DataFrame(scaled, columns=feature_report['features']).to_csv(
        output / 'F3_atomic_umap_standardized.csv', index=False)
    results.groupby('Cluster')[feature_report['features']].agg(['count', 'mean', 'std']).to_csv(
        output / 'F3_atomic_umap_cluster_summary.csv')
    fig, ax = plt.subplots(figsize=(8, 6))
    points = ax.scatter(*visual.T, c=labels, cmap='tab20', s=3, alpha=0.6, rasterized=True)
    fig.colorbar(points, ax=ax, label='HDBSCAN cluster (-1 = noise)')
    ax.set(xlabel='UMAP 1', ylabel='UMAP 2', title='Atomic descriptor clusters (fit in 5-D UMAP)')
    fig.tight_layout()
    fig.savefig(output / 'F3_atomic_umap_clusters.png', dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
    for row, scheme in enumerate(('AW', 'Power')):
        for column, metric in enumerate(('Volume', 'Surface Area', 'Number of Neighbors')):
            ax = axes[row, column]
            name = f'{scheme} {metric}'
            if name in results:
                points = ax.scatter(*visual.T, c=results[name], cmap='viridis', s=2,
                                    alpha=0.6, rasterized=True)
                fig.colorbar(points, ax=ax)
            ax.set(title=f'{scheme}: {"Contacts" if metric == "Number of Neighbors" else metric}',
                   xlabel='UMAP 1', ylabel='UMAP 2')
    fig.tight_layout()
    fig.savefig(output / 'F3_atomic_umap_metrics.png', dpi=300)
    plt.close(fig)
    report = dict(feature_report, collection=collection, atoms=len(results),
                  neighbors_requested=neighbors, neighbors_used=min(neighbors, len(results) - 1),
                  min_cluster_size=min_cluster_size, seed=seed, clustering_dimensions=5,
                  noise_atoms=int((labels == -1).sum()),
                  packages={name: version(name) for name in ('umap-learn', 'scikit-learn', 'numpy')})
    (output / 'F3_atomic_umap_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return results


def analyze(data_root, output, *, exclude_keys=('A', 'B', 'C'), neighbors=30,
            min_cluster_size=15, seed=42, missing='median'):
    atoms, collection = collect_descriptors(data_root, exclude_keys)
    metadata, scaled, report = prepare_features(atoms, missing)
    return export_results(metadata, scaled, report, collection, Path(output), neighbors=neighbors,
                          min_cluster_size=min_cluster_size, seed=seed)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--exclude-keys', nargs='*', default=['A', 'B', 'C'])
    parser.add_argument('--neighbors', type=int, default=30)
    parser.add_argument('--min-cluster-size', type=int, default=15)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--missing', choices=['median', 'complete'], default='median')
    args = parser.parse_args(argv)
    if args.neighbors < 2 or args.min_cluster_size < 2:
        parser.error('--neighbors and --min-cluster-size must be at least 2')
    analyze(args.data_root, args.output, exclude_keys=args.exclude_keys, neighbors=args.neighbors,
            min_cluster_size=args.min_cluster_size, seed=args.seed, missing=args.missing)
    print(f'UMAP results written to {args.output.resolve()}')


if __name__ == '__main__':
    main()
