
import os
import sys
import math
import tkinter as tk
from tkinter import filedialog
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


SOLVENT_RESNAMES = {'SOL', 'HOH', 'WAT', 'TIP3', 'TIP3P', 'SPC', 'SPCE', 'OPC', 'H2O'}


def choose_folder() -> Optional[str]:
    root = tk.Tk()
    root.withdraw()
    folder = filedialog.askdirectory(
        title='Select folder containing patch_atoms.csv, aw_logs.csv, pow_logs.csv, and a PDB file'
    )
    root.destroy()
    return folder if folder else None


def find_first_existing(folder: str, candidates: List[str]) -> Optional[str]:
    for rel in candidates:
        path = os.path.join(folder, rel)
        if os.path.exists(path):
            return path
    return None


def require_path(path: Optional[str], description: str) -> str:
    if path is None:
        raise FileNotFoundError(f'Could not find {description}.')
    return path


def load_patch_atoms_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)

    if 'Index' not in df.columns:
        if len(df.columns) == 1:
            df.columns = ['Index']
        else:
            raise ValueError("Patch atom CSV must contain an 'Index' column or be a single-column file.")

    df = df.copy()
    df['Index'] = pd.to_numeric(df['Index'], errors='coerce')
    df = df.dropna(subset=['Index']).copy()
    df['Index'] = df['Index'].astype(int)
    df = df.drop_duplicates(subset=['Index']).sort_values('Index').reset_index(drop=True)
    return df


def parse_pdb_atoms(pdb_path: str) -> pd.DataFrame:
    rows = []

    with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as handle:
        pdb_row_index = 0

        for line in handle:
            record = line[:6].strip()
            if record not in {'ATOM', 'HETATM'}:
                continue

            serial_text = line[6:11].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21:22].strip()
            res_seq_text = line[22:26].strip()
            x_text = line[30:38].strip()
            y_text = line[38:46].strip()
            z_text = line[46:54].strip()
            element = line[76:78].strip()

            try:
                serial = int(serial_text) if serial_text else np.nan
            except ValueError:
                serial = np.nan

            try:
                res_seq = int(res_seq_text) if res_seq_text else np.nan
            except ValueError:
                res_seq = np.nan

            try:
                x = float(x_text)
                y = float(y_text)
                z = float(z_text)
            except ValueError:
                x = y = z = np.nan

            if not element:
                element = ''.join([c for c in atom_name if c.isalpha()])[:2].strip().title()

            rows.append({
                'PDB Row Index': pdb_row_index,
                'PDB Serial': serial,
                'PDB Atom Name': atom_name,
                'PDB Residue': res_name,
                'PDB Chain': chain,
                'PDB Residue Sequence': res_seq,
                'PDB X': x,
                'PDB Y': y,
                'PDB Z': z,
                'PDB Element': element,
                'Is SOL': str(res_name).strip().upper() in SOLVENT_RESNAMES,
                'Record Type': record,
            })
            pdb_row_index += 1

    if not rows:
        raise ValueError(f'No ATOM/HETATM records found in {pdb_path}')

    return pd.DataFrame(rows)


def normalize_atom_df(atom_df: pd.DataFrame) -> pd.DataFrame:
    df = atom_df.copy()

    for col in ['Index', 'Residue Sequence']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    for col in ['X', 'Y', 'Z']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    for col in ['Name', 'Residue', 'Chain']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    return df


def _norm_str(value) -> str:
    if pd.isna(value):
        return ''
    return str(value).strip().upper()


def _coord_key(x: float, y: float, z: float, decimals: int = 3) -> Tuple[float, float, float]:
    return (round(float(x), decimals), round(float(y), decimals), round(float(z), decimals))


def build_pdb_index_map(log_atoms: pd.DataFrame, pdb_df: pd.DataFrame, coord_decimals: int = 3) -> Tuple[Dict[int, Optional[int]], pd.DataFrame]:
    """
    Build a robust map from log Index -> PDB Row Index.

    Matching order:
      1. exact rounded coordinate key + exact atom/residue/chain/resseq identity
      2. exact rounded coordinate key only
      3. nearest coordinate within tolerance with identity preference
    """
    log_atoms = normalize_atom_df(log_atoms)
    pdb_df = pdb_df.copy()

    coord_buckets: Dict[Tuple[float, float, float], List[int]] = {}
    for _, row in pdb_df.iterrows():
        if pd.isna(row['PDB X']) or pd.isna(row['PDB Y']) or pd.isna(row['PDB Z']):
            continue
        key = _coord_key(row['PDB X'], row['PDB Y'], row['PDB Z'], decimals=coord_decimals)
        coord_buckets.setdefault(key, []).append(int(row['PDB Row Index']))

    pdb_small = pdb_df.set_index('PDB Row Index', drop=False)
    mapping: Dict[int, Optional[int]] = {}
    diagnostics = []

    pdb_coords = pdb_df[['PDB X', 'PDB Y', 'PDB Z']].to_numpy(dtype=float)

    for _, atom in log_atoms.iterrows():
        log_idx = int(atom['Index'])
        log_name = _norm_str(atom.get('Name'))
        log_res = _norm_str(atom.get('Residue'))
        log_chain = _norm_str(atom.get('Chain'))
        log_resseq = atom.get('Residue Sequence', np.nan)

        match_idx = None
        match_method = 'unmatched'
        dist = np.nan

        if not any(pd.isna(atom[c]) for c in ['X', 'Y', 'Z']):
            key = _coord_key(atom['X'], atom['Y'], atom['Z'], decimals=coord_decimals)
            candidates = coord_buckets.get(key, [])

            if candidates:
                best_identity = []
                for cand_idx in candidates:
                    cand = pdb_small.loc[cand_idx]
                    same_name = _norm_str(cand['PDB Atom Name']) == log_name
                    same_res = _norm_str(cand['PDB Residue']) == log_res
                    same_chain = _norm_str(cand['PDB Chain']) == log_chain
                    same_resseq = (
                        (pd.isna(log_resseq) and pd.isna(cand['PDB Residue Sequence']))
                        or (not pd.isna(log_resseq) and not pd.isna(cand['PDB Residue Sequence']) and int(log_resseq) == int(cand['PDB Residue Sequence']))
                    )
                    score = int(same_name) + int(same_res) + int(same_chain) + int(same_resseq)
                    best_identity.append((score, cand_idx))

                best_identity.sort(reverse=True)
                best_score, best_idx = best_identity[0]
                match_idx = int(best_idx)
                match_method = 'coord+identity' if best_score >= 2 else 'coord-only'
                cand = pdb_small.loc[match_idx]
                dist = math.sqrt(
                    (float(atom['X']) - float(cand['PDB X'])) ** 2 +
                    (float(atom['Y']) - float(cand['PDB Y'])) ** 2 +
                    (float(atom['Z']) - float(cand['PDB Z'])) ** 2
                )

            else:
                diffs = pdb_coords - np.array([[float(atom['X']), float(atom['Y']), float(atom['Z'])]])
                dists = np.sqrt(np.sum(diffs ** 2, axis=1))
                if len(dists):
                    nearest_pos = int(np.argmin(dists))
                    nearest_dist = float(dists[nearest_pos])

                    if nearest_dist <= 0.05:
                        candidate = pdb_df.iloc[nearest_pos]
                        match_idx = int(candidate['PDB Row Index'])
                        match_method = 'nearest<=0.05'
                        dist = nearest_dist

        mapping[log_idx] = match_idx

        if match_idx is None:
            diagnostics.append({
                'Index': log_idx,
                'Match Method': match_method,
                'Distance': dist,
                'Log Name': atom.get('Name'),
                'Log Residue': atom.get('Residue'),
                'Log Chain': atom.get('Chain'),
                'Log Residue Sequence': atom.get('Residue Sequence'),
                'PDB Row Index': np.nan,
                'PDB Atom Name': np.nan,
                'PDB Residue': np.nan,
                'PDB Chain': np.nan,
                'PDB Residue Sequence': np.nan,
                'Is SOL': np.nan,
            })
        else:
            cand = pdb_small.loc[match_idx]
            diagnostics.append({
                'Index': log_idx,
                'Match Method': match_method,
                'Distance': dist,
                'Log Name': atom.get('Name'),
                'Log Residue': atom.get('Residue'),
                'Log Chain': atom.get('Chain'),
                'Log Residue Sequence': atom.get('Residue Sequence'),
                'PDB Row Index': match_idx,
                'PDB Atom Name': cand['PDB Atom Name'],
                'PDB Residue': cand['PDB Residue'],
                'PDB Chain': cand['PDB Chain'],
                'PDB Residue Sequence': cand['PDB Residue Sequence'],
                'Is SOL': cand['Is SOL'],
            })

    return mapping, pd.DataFrame(diagnostics)


def append_pdb_labels(atom_df: pd.DataFrame, pdb_df: pd.DataFrame, index_map: Dict[int, Optional[int]]) -> pd.DataFrame:
    atom_df = atom_df.copy()
    pdb_small = pdb_df.set_index('PDB Row Index', drop=False)

    def fetch(log_idx: int, key: str):
        pdb_idx = index_map.get(int(log_idx), None)
        if pdb_idx is None or pdb_idx not in pdb_small.index:
            return np.nan
        return pdb_small.loc[pdb_idx, key]

    atom_df['PDB Row Index'] = atom_df['Index'].apply(lambda idx: index_map.get(int(idx), None))
    atom_df['PDB Atom Name'] = atom_df['Index'].apply(lambda idx: fetch(idx, 'PDB Atom Name'))
    atom_df['PDB Residue'] = atom_df['Index'].apply(lambda idx: fetch(idx, 'PDB Residue'))
    atom_df['PDB Chain'] = atom_df['Index'].apply(lambda idx: fetch(idx, 'PDB Chain'))
    atom_df['PDB Residue Sequence'] = atom_df['Index'].apply(lambda idx: fetch(idx, 'PDB Residue Sequence'))
    atom_df['PDB Serial'] = atom_df['Index'].apply(lambda idx: fetch(idx, 'PDB Serial'))
    atom_df['Is SOL'] = atom_df['Index'].apply(
        lambda idx: bool(fetch(idx, 'Is SOL')) if pd.notna(fetch(idx, 'Is SOL')) else False
    )

    return atom_df


def build_surface_lookup(surfs_df: pd.DataFrame) -> Dict[Tuple[int, int], List[Dict]]:
    lookup: Dict[Tuple[int, int], List[Dict]] = {}
    if surfs_df is None or len(surfs_df) == 0:
        return lookup

    for _, row in surfs_df.iterrows():
        balls = row.get('Balls', None)
        if not isinstance(balls, list) or len(balls) != 2:
            continue
        key = tuple(sorted((int(balls[0]), int(balls[1]))))
        lookup.setdefault(key, []).append(row.to_dict())

    return lookup


def get_shared_surfaces(surface_lookup: Dict[Tuple[int, int], List[Dict]], idx1: int, idx2: int) -> List[Dict]:
    return surface_lookup.get(tuple(sorted((int(idx1), int(idx2)))), [])


def summarize_shared_surfaces(shared_surfs: List[Dict]) -> Dict[str, float]:
    if not shared_surfs:
        return {
            'Shared Surface Count': 0,
            'Shared Surface Area Total': 0.0,
            'Shared Surface Area Mean': np.nan,
            'Shared Surface Area Max': np.nan,
            'Shared Contact Area Total': 0.0,
            'Shared Overlap Total': 0.0,
            'Shared Mean Curvature Mean': np.nan,
            'Shared Mean Curvature MaxAbs': np.nan,
        }

    area_vals = [float(s.get('Surface Area', np.nan)) for s in shared_surfs]
    contact_vals = [float(s.get('Contact Area', np.nan)) for s in shared_surfs]
    overlap_vals = [float(s.get('Overlap', np.nan)) for s in shared_surfs]

    curv_vals = []
    for s in shared_surfs:
        if 'Mean Curvature' in s and pd.notna(s['Mean Curvature']):
            curv_vals.append(float(s['Mean Curvature']))
        elif 'Curvature' in s and pd.notna(s['Curvature']):
            curv_vals.append(float(s['Curvature']))

    return {
        'Shared Surface Count': len(shared_surfs),
        'Shared Surface Area Total': float(np.nansum(area_vals)),
        'Shared Surface Area Mean': float(np.nanmean(area_vals)) if len(area_vals) else np.nan,
        'Shared Surface Area Max': float(np.nanmax(area_vals)) if len(area_vals) else np.nan,
        'Shared Contact Area Total': float(np.nansum(contact_vals)),
        'Shared Overlap Total': float(np.nansum(overlap_vals)),
        'Shared Mean Curvature Mean': float(np.nanmean(curv_vals)) if len(curv_vals) else np.nan,
        'Shared Mean Curvature MaxAbs': float(np.nanmax(np.abs(curv_vals))) if len(curv_vals) else np.nan,
    }


def collect_neighbor_pair_rows(
    patch_indices: List[int],
    scheme_name: str,
    atom_df: pd.DataFrame,
    surface_lookup: Dict[Tuple[int, int], List[Dict]],
    pdb_df: pd.DataFrame,
) -> pd.DataFrame:
    atom_indexed = atom_df.set_index('Index', drop=False)
    pdb_row_lookup = build_pdb_row_lookup(pdb_df)
    pdb_serial_lookup = build_pdb_serial_lookup(pdb_df)

    rows = []

    for patch_idx in patch_indices:
        if patch_idx not in atom_indexed.index:
            continue

        patch_atom = atom_indexed.loc[patch_idx]
        neighbors = patch_atom.get('Neighbors', [])

        if not isinstance(neighbors, list):
            continue

        for nbr_idx in neighbors:
            nbr_info = resolve_neighbor_from_log_or_pdb(
                nbr_idx=nbr_idx,
                atom_indexed=atom_indexed,
                pdb_row_lookup=pdb_row_lookup,
                pdb_serial_lookup=pdb_serial_lookup,
            )

            surf_stats = summarize_shared_surfaces(
                get_shared_surfaces(surface_lookup, patch_idx, nbr_idx)
            )

            row = {
                'Scheme': scheme_name,
                'Patch Index': int(patch_idx),
                'Patch Atom Name': patch_atom.get('Name'),
                'Patch Residue': patch_atom.get('Residue'),
                'Patch Chain': patch_atom.get('Chain'),
                'Patch Residue Sequence': patch_atom.get('Residue Sequence'),
                'Patch Volume': patch_atom.get('Volume'),
                'Patch Surface Area': patch_atom.get('Surface Area'),
                'Patch Number of Neighbors': patch_atom.get('Number of Neighbors'),
                'Patch Is SOL': bool(patch_atom.get('Is SOL', False)),
                'Patch PDB Residue': patch_atom.get('PDB Residue'),
                'Patch PDB Atom Name': patch_atom.get('PDB Atom Name'),
            }

            row.update(nbr_info)
            row['Neighbor Type'] = 'SOL' if bool(row['Neighbor Is SOL']) else 'NON-SOL'
            row['Patch-Neighbor Type'] = (
                f"{'SOL' if bool(row['Patch Is SOL']) else 'NON-SOL'}--"
                f"{'SOL' if bool(row['Neighbor Is SOL']) else 'NON-SOL'}"
            )
            row.update(surf_stats)

            rows.append(row)

    return pd.DataFrame(rows)


def build_pdb_row_lookup(pdb_df: pd.DataFrame) -> Dict[int, Dict]:
    """
    Map 0-based PDB atom row index -> PDB atom metadata.
    """
    lookup = {}
    for _, row in pdb_df.iterrows():
        lookup[int(row['PDB Row Index'])] = row.to_dict()
    return lookup


def build_pdb_serial_lookup(pdb_df: pd.DataFrame) -> Dict[int, Dict]:
    """
    Map PDB serial number -> PDB atom metadata.
    """
    lookup = {}
    for _, row in pdb_df.iterrows():
        serial = row.get('PDB Serial', None)
        if pd.notna(serial):
            lookup[int(serial)] = row.to_dict()
    return lookup


def resolve_neighbor_from_log_or_pdb(
    nbr_idx: int,
    atom_indexed: pd.DataFrame,
    pdb_row_lookup: Dict[int, Dict],
    pdb_serial_lookup: Dict[int, Dict],
) -> Dict:
    """
    Resolve neighbor identity.

    Priority:
    1. log atom table by Index
    2. PDB row index fallback (0-based)
    3. PDB serial fallback
    """
    nbr_idx = int(nbr_idx)

    if nbr_idx in atom_indexed.index:
        row = atom_indexed.loc[nbr_idx]
        return {
            'Neighbor Found': True,
            'Neighbor Source': 'LOG',
            'Neighbor Index': nbr_idx,
            'Neighbor Atom Name': row.get('Name'),
            'Neighbor Residue': row.get('Residue'),
            'Neighbor Chain': row.get('Chain'),
            'Neighbor Residue Sequence': row.get('Residue Sequence'),
            'Neighbor Volume': row.get('Volume'),
            'Neighbor Surface Area': row.get('Surface Area'),
            'Neighbor Number of Neighbors': row.get('Number of Neighbors'),
            'Neighbor Is SOL': bool(row.get('Is SOL', False)),
            'Neighbor PDB Atom Name': row.get('PDB Atom Name'),
            'Neighbor PDB Residue': row.get('PDB Residue'),
            'Neighbor PDB Chain': row.get('PDB Chain'),
            'Neighbor PDB Residue Sequence': row.get('PDB Residue Sequence'),
            'Neighbor PDB Serial': row.get('PDB Serial'),
        }

    if nbr_idx in pdb_row_lookup:
        row = pdb_row_lookup[nbr_idx]
        return {
            'Neighbor Found': True,
            'Neighbor Source': 'PDB_ROW',
            'Neighbor Index': nbr_idx,
            'Neighbor Atom Name': row.get('PDB Atom Name'),
            'Neighbor Residue': row.get('PDB Residue'),
            'Neighbor Chain': row.get('PDB Chain'),
            'Neighbor Residue Sequence': row.get('PDB Residue Sequence'),
            'Neighbor Volume': np.nan,
            'Neighbor Surface Area': np.nan,
            'Neighbor Number of Neighbors': np.nan,
            'Neighbor Is SOL': bool(row.get('Is SOL', False)),
            'Neighbor PDB Atom Name': row.get('PDB Atom Name'),
            'Neighbor PDB Residue': row.get('PDB Residue'),
            'Neighbor PDB Chain': row.get('PDB Chain'),
            'Neighbor PDB Residue Sequence': row.get('PDB Residue Sequence'),
            'Neighbor PDB Serial': row.get('PDB Serial'),
        }

    if nbr_idx in pdb_serial_lookup:
        row = pdb_serial_lookup[nbr_idx]
        return {
            'Neighbor Found': True,
            'Neighbor Source': 'PDB_SERIAL',
            'Neighbor Index': nbr_idx,
            'Neighbor Atom Name': row.get('PDB Atom Name'),
            'Neighbor Residue': row.get('PDB Residue'),
            'Neighbor Chain': row.get('PDB Chain'),
            'Neighbor Residue Sequence': row.get('PDB Residue Sequence'),
            'Neighbor Volume': np.nan,
            'Neighbor Surface Area': np.nan,
            'Neighbor Number of Neighbors': np.nan,
            'Neighbor Is SOL': bool(row.get('Is SOL', False)),
            'Neighbor PDB Atom Name': row.get('PDB Atom Name'),
            'Neighbor PDB Residue': row.get('PDB Residue'),
            'Neighbor PDB Chain': row.get('PDB Chain'),
            'Neighbor PDB Residue Sequence': row.get('PDB Residue Sequence'),
            'Neighbor PDB Serial': row.get('PDB Serial'),
        }

    return {
        'Neighbor Found': False,
        'Neighbor Source': 'MISSING',
        'Neighbor Index': nbr_idx,
        'Neighbor Atom Name': None,
        'Neighbor Residue': None,
        'Neighbor Chain': None,
        'Neighbor Residue Sequence': None,
        'Neighbor Volume': np.nan,
        'Neighbor Surface Area': np.nan,
        'Neighbor Number of Neighbors': np.nan,
        'Neighbor Is SOL': False,
        'Neighbor PDB Atom Name': None,
        'Neighbor PDB Residue': None,
        'Neighbor PDB Chain': None,
        'Neighbor PDB Residue Sequence': None,
        'Neighbor PDB Serial': None,
    }


def summarize_patch_atoms(
    patch_indices: List[int],
    aw_atoms: pd.DataFrame,
    pow_atoms: pd.DataFrame,
    aw_pairs: pd.DataFrame,
    pow_pairs: pd.DataFrame,
) -> pd.DataFrame:
    aw_indexed = aw_atoms.set_index('Index', drop=False)
    pow_indexed = pow_atoms.set_index('Index', drop=False)

    rows = []
    for idx in patch_indices:
        if idx not in aw_indexed.index or idx not in pow_indexed.index:
            continue

        aw_atom = aw_indexed.loc[idx]
        pow_atom = pow_indexed.loc[idx]

        aw_sub = aw_pairs[aw_pairs['Patch Index'] == idx].copy() if len(aw_pairs) else pd.DataFrame()
        pow_sub = pow_pairs[pow_pairs['Patch Index'] == idx].copy() if len(pow_pairs) else pd.DataFrame()

        row = {
            'Index': int(idx),
            'Atom Name': aw_atom.get('Name'),
            'Residue': aw_atom.get('Residue'),
            'Chain': aw_atom.get('Chain'),
            'Residue Sequence': aw_atom.get('Residue Sequence'),
            'PDB Residue': aw_atom.get('PDB Residue'),
            'Is SOL': bool(aw_atom.get('Is SOL', False)),
            'AW Volume': aw_atom.get('Volume'),
            'Pow Volume': pow_atom.get('Volume'),
            'Delta Volume': float(pow_atom.get('Volume', np.nan) - aw_atom.get('Volume', np.nan)),
            'Delta Volume Abs': float(abs(pow_atom.get('Volume', np.nan) - aw_atom.get('Volume', np.nan))),
            'AW Surface Area': aw_atom.get('Surface Area'),
            'Pow Surface Area': pow_atom.get('Surface Area'),
            'Delta Surface Area': float(pow_atom.get('Surface Area', np.nan) - aw_atom.get('Surface Area', np.nan)),
            'AW Neighbor Count': aw_atom.get('Number of Neighbors'),
            'Pow Neighbor Count': pow_atom.get('Number of Neighbors'),
            'Delta Neighbor Count': float(pow_atom.get('Number of Neighbors', np.nan) - aw_atom.get('Number of Neighbors', np.nan)),
            'AW Neighbor Surface Area Sum': float(aw_sub['Neighbor Surface Area'].sum()) if len(aw_sub) else 0.0,
            'Pow Neighbor Surface Area Sum': float(pow_sub['Neighbor Surface Area'].sum()) if len(pow_sub) else 0.0,
            'AW Shared Surface Area Sum': float(aw_sub['Shared Surface Area Total'].sum()) if len(aw_sub) else 0.0,
            'Pow Shared Surface Area Sum': float(pow_sub['Shared Surface Area Total'].sum()) if len(pow_sub) else 0.0,
            'AW Shared Contact Area Sum': float(aw_sub['Shared Contact Area Total'].sum()) if len(aw_sub) else 0.0,
            'Pow Shared Contact Area Sum': float(pow_sub['Shared Contact Area Total'].sum()) if len(pow_sub) else 0.0,
            'AW SOL Neighbor Count': int((aw_sub['Neighbor Is SOL'] == True).sum()) if len(aw_sub) else 0,
            'Pow SOL Neighbor Count': int((pow_sub['Neighbor Is SOL'] == True).sum()) if len(pow_sub) else 0,
            'AW Non-SOL Neighbor Count': int((aw_sub['Neighbor Is SOL'] == False).sum()) if len(aw_sub) else 0,
            'Pow Non-SOL Neighbor Count': int((pow_sub['Neighbor Is SOL'] == False).sum()) if len(pow_sub) else 0,
            'AW SOL Shared Surface Area Sum': float(aw_sub.loc[aw_sub['Neighbor Is SOL'] == True, 'Shared Surface Area Total'].sum()) if len(aw_sub) else 0.0,
            'Pow SOL Shared Surface Area Sum': float(pow_sub.loc[pow_sub['Neighbor Is SOL'] == True, 'Shared Surface Area Total'].sum()) if len(pow_sub) else 0.0,
            'AW Non-SOL Shared Surface Area Sum': float(aw_sub.loc[aw_sub['Neighbor Is SOL'] == False, 'Shared Surface Area Total'].sum()) if len(aw_sub) else 0.0,
            'Pow Non-SOL Shared Surface Area Sum': float(pow_sub.loc[pow_sub['Neighbor Is SOL'] == False, 'Shared Surface Area Total'].sum()) if len(pow_sub) else 0.0,
        }
        rows.append(row)

    return pd.DataFrame(rows)


def compare_aw_pow_pairs(aw_pairs: pd.DataFrame, pow_pairs: pd.DataFrame) -> pd.DataFrame:
    key_cols = ['Patch Index', 'Neighbor Index']

    aw_cmp = aw_pairs.copy()
    pow_cmp = pow_pairs.copy()

    aw_keep = key_cols + [
        'Patch Atom Name', 'Patch Residue', 'Neighbor Atom Name', 'Neighbor Residue',
        'Neighbor PDB Residue', 'Neighbor Is SOL', 'Neighbor Type',
        'Shared Surface Count', 'Shared Surface Area Total', 'Shared Contact Area Total',
        'Shared Overlap Total', 'Shared Mean Curvature Mean', 'Shared Mean Curvature MaxAbs'
    ]
    pow_keep = key_cols + [
        'Shared Surface Count', 'Shared Surface Area Total', 'Shared Contact Area Total',
        'Shared Overlap Total', 'Shared Mean Curvature Mean', 'Shared Mean Curvature MaxAbs'
    ]

    aw_cmp = aw_cmp[aw_keep].rename(columns={
        'Shared Surface Count': 'AW Shared Surface Count',
        'Shared Surface Area Total': 'AW Shared Surface Area Total',
        'Shared Contact Area Total': 'AW Shared Contact Area Total',
        'Shared Overlap Total': 'AW Shared Overlap Total',
        'Shared Mean Curvature Mean': 'AW Shared Mean Curvature Mean',
        'Shared Mean Curvature MaxAbs': 'AW Shared Mean Curvature MaxAbs',
    })

    pow_cmp = pow_cmp[pow_keep].rename(columns={
        'Shared Surface Count': 'Pow Shared Surface Count',
        'Shared Surface Area Total': 'Pow Shared Surface Area Total',
        'Shared Contact Area Total': 'Pow Shared Contact Area Total',
        'Shared Overlap Total': 'Pow Shared Overlap Total',
        'Shared Mean Curvature Mean': 'Pow Shared Mean Curvature Mean',
        'Shared Mean Curvature MaxAbs': 'Pow Shared Mean Curvature MaxAbs',
    })

    merged = aw_cmp.merge(pow_cmp, on=key_cols, how='outer')

    numeric_pairs = [
        ('AW Shared Surface Count', 'Pow Shared Surface Count', 'Delta Shared Surface Count'),
        ('AW Shared Surface Area Total', 'Pow Shared Surface Area Total', 'Delta Shared Surface Area Total'),
        ('AW Shared Contact Area Total', 'Pow Shared Contact Area Total', 'Delta Shared Contact Area Total'),
        ('AW Shared Overlap Total', 'Pow Shared Overlap Total', 'Delta Shared Overlap Total'),
        ('AW Shared Mean Curvature Mean', 'Pow Shared Mean Curvature Mean', 'Delta Shared Mean Curvature Mean'),
        ('AW Shared Mean Curvature MaxAbs', 'Pow Shared Mean Curvature MaxAbs', 'Delta Shared Mean Curvature MaxAbs'),
    ]

    for left, right, out in numeric_pairs:
        merged[left] = pd.to_numeric(merged[left], errors='coerce').fillna(0.0)
        merged[right] = pd.to_numeric(merged[right], errors='coerce').fillna(0.0)
        merged[out] = merged[right] - merged[left]

    return merged


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


def debug_patch_neighbors(
    patch_indices: List[int],
    aw_atoms: pd.DataFrame,
    pow_atoms: pd.DataFrame,
    max_patch_atoms: int = 5,
    max_neighbors_per_patch: int = 20,
) -> None:
    aw_indexed = aw_atoms.set_index('Index', drop=False)
    pow_indexed = pow_atoms.set_index('Index', drop=False)

    print("\n" + "=" * 80)
    print("DEBUG NEIGHBOR CHECK")
    print("=" * 80)

    shown = 0

    for patch_idx in patch_indices:
        if patch_idx not in aw_indexed.index:
            print(f"\nPatch index {patch_idx} not found in AW atoms.")
            continue

        aw_patch = aw_indexed.loc[patch_idx]

        print("\n" + "-" * 80)
        print(f"PATCH INDEX: {patch_idx}")
        print(
            f"AW PATCH: Name={aw_patch.get('Name')}, "
            f"Residue={aw_patch.get('Residue')}, "
            f"Chain={aw_patch.get('Chain')}, "
            f"ResSeq={aw_patch.get('Residue Sequence')}, "
            f"IsSOL={aw_patch.get('Is SOL')}, "
            f"PDB_Residue={aw_patch.get('PDB Residue')}, "
            f"PDB_Atom={aw_patch.get('PDB Atom Name')}, "
            f"PDB_Serial={aw_patch.get('PDB Serial')}"
        )

        aw_neighbors = aw_patch.get('Neighbors', [])
        print(f"RAW AW NEIGHBORS FIELD: {aw_neighbors}")
        print(f"AW NUMBER OF NEIGHBORS FIELD: {aw_patch.get('Number of Neighbors')}")

        if not isinstance(aw_neighbors, list):
            print("AW neighbors is not a list. Skipping.")
            continue

        if len(aw_neighbors) == 0:
            print("No AW neighbors recorded.")
        else:
            print("\nAW NEIGHBOR DETAILS:")
            for nbr_idx in aw_neighbors[:max_neighbors_per_patch]:
                if nbr_idx in aw_indexed.index:
                    nbr = aw_indexed.loc[nbr_idx]
                    print(
                        f"  Neighbor Index={nbr_idx} | "
                        f"LogName={nbr.get('Name')} | "
                        f"LogResidue={nbr.get('Residue')} | "
                        f"Chain={nbr.get('Chain')} | "
                        f"ResSeq={nbr.get('Residue Sequence')} | "
                        f"PDBResidue={nbr.get('PDB Residue')} | "
                        f"PDBAtom={nbr.get('PDB Atom Name')} | "
                        f"PDBSerial={nbr.get('PDB Serial')} | "
                        f"IsSOL={nbr.get('Is SOL')}"
                    )
                else:
                    print(f"  Neighbor Index={nbr_idx} | NOT FOUND IN AW ATOM TABLE")

        if patch_idx in pow_indexed.index:
            pow_patch = pow_indexed.loc[patch_idx]
            pow_neighbors = pow_patch.get('Neighbors', [])

            print(f"\nRAW POW NEIGHBORS FIELD: {pow_neighbors}")
            print(f"POW NUMBER OF NEIGHBORS FIELD: {pow_patch.get('Number of Neighbors')}")

            if isinstance(pow_neighbors, list) and len(pow_neighbors) > 0:
                print("\nPOW NEIGHBOR DETAILS:")
                for nbr_idx in pow_neighbors[:max_neighbors_per_patch]:
                    if nbr_idx in pow_indexed.index:
                        nbr = pow_indexed.loc[nbr_idx]
                        print(
                            f"  Neighbor Index={nbr_idx} | "
                            f"LogName={nbr.get('Name')} | "
                            f"LogResidue={nbr.get('Residue')} | "
                            f"Chain={nbr.get('Chain')} | "
                            f"ResSeq={nbr.get('Residue Sequence')} | "
                            f"PDBResidue={nbr.get('PDB Residue')} | "
                            f"PDBAtom={nbr.get('PDB Atom Name')} | "
                            f"PDBSerial={nbr.get('PDB Serial')} | "
                            f"IsSOL={nbr.get('Is SOL')}"
                        )
                    else:
                        print(f"  Neighbor Index={nbr_idx} | NOT FOUND IN POW ATOM TABLE")

        shown += 1
        if shown >= max_patch_atoms:
            break


def print_missing_neighbor_resolution(pairs_df: pd.DataFrame, max_rows: int = 50) -> None:
    sub = pairs_df[pairs_df['Neighbor Source'] != 'LOG'].copy()

    print("\n=== MISSING-NEIGHBOR / PDB-FALLBACK CHECK ===")
    if len(sub) == 0:
        print("No neighbors required fallback resolution.")
        return

    cols = [
        'Patch Index',
        'Neighbor Index',
        'Neighbor Source',
        'Neighbor Atom Name',
        'Neighbor Residue',
        'Neighbor Chain',
        'Neighbor Residue Sequence',
        'Neighbor PDB Atom Name',
        'Neighbor PDB Residue',
        'Neighbor PDB Serial',
        'Neighbor Is SOL',
    ]
    cols = [c for c in cols if c in sub.columns]

    print(sub[cols].head(max_rows).to_string(index=False))


def compute_correlations(df: pd.DataFrame, target_col: str, feature_cols: List[str]) -> pd.DataFrame:
    rows = []
    for feature in feature_cols:
        if feature not in df.columns:
            continue
        rows.append({
            'Feature': feature,
            'N': int(pd.concat([df[target_col], df[feature]], axis=1).dropna().shape[0]),
            'Pearson': safe_corr(df[target_col], df[feature], method='pearson'),
            'Spearman': safe_corr(df[target_col], df[feature], method='spearman'),
        })

    out = pd.DataFrame(rows)
    if len(out):
        out['Abs Pearson'] = out['Pearson'].abs()
        out = out.sort_values(['Abs Pearson', 'Feature'], ascending=[False, True]).reset_index(drop=True)

    return out


def robust_outlier_flags(series: pd.Series, z_cutoff: float = 3.5) -> pd.Series:
    values = pd.to_numeric(series, errors='coerce')
    med = values.median()
    mad = (values - med).abs().median()

    if pd.isna(mad) or mad == 0:
        return pd.Series([False] * len(series), index=series.index)

    robust_z = 0.6745 * (values - med) / mad
    return robust_z.abs() >= z_cutoff


def summarize_outliers_by_group(summary_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    df = summary_df.copy()
    df['Is Outlier'] = robust_outlier_flags(df['Delta Volume Abs'])

    grouped = (
        df.groupby(group_col)
        .agg(
            Count=('Index', 'size'),
            Outlier_Count=('Is Outlier', 'sum'),
            Mean_Delta_Volume_Abs=('Delta Volume Abs', 'mean'),
            Median_Delta_Volume_Abs=('Delta Volume Abs', 'median'),
        )
        .reset_index()
    )
    grouped['Outlier_Fraction'] = grouped['Outlier_Count'] / grouped['Count']
    return grouped.sort_values(['Outlier_Fraction', 'Count'], ascending=[False, False]).reset_index(drop=True)


def main() -> None:
    folder = choose_folder()
    if not folder:
        print('No folder selected.')
        return

    patch_path = require_path(find_first_existing(folder, ['patch_atoms.csv']), 'patch_atoms.csv')
    aw_path = require_path(find_first_existing(folder, ['aw_logs.csv', os.path.join('aw', 'aw_logs.csv')]), 'aw_logs.csv')
    pow_path = require_path(find_first_existing(folder, ['pow_logs.csv', os.path.join('pow', 'pow_logs.csv')]), 'pow_logs.csv')

    pdb_candidates = ['structure.pdb', 'model.pdb', 'input.pdb'] + [f for f in os.listdir(folder) if f.lower().endswith('.pdb')]
    pdb_path = require_path(find_first_existing(folder, pdb_candidates), 'a PDB file')

    print('\n=== INPUTS ===')
    print(f'Folder: {folder}')
    print(f'Patch atoms: {patch_path}')
    print(f'AW logs: {aw_path}')
    print(f'Pow logs: {pow_path}')
    print(f'PDB: {pdb_path}')

    patch_df = load_patch_atoms_csv(patch_path)
    patch_indices = patch_df['Index'].tolist()

    aw_logs = read_logs2(aw_path, all_=False, balls=True, surfs=True)
    pow_logs = read_logs2(pow_path, all_=False, balls=True, surfs=True)

    aw_atoms = normalize_atom_df(aw_logs['atoms'].copy())
    pow_atoms = normalize_atom_df(pow_logs['atoms'].copy())
    aw_surfs = aw_logs['surfs'].copy()
    pow_surfs = pow_logs['surfs'].copy()

    pdb_df = parse_pdb_atoms(pdb_path)

    # Build mapping from AW atoms to PDB, then reuse by log Index for Pow.
    pdb_index_map, mapping_diag = build_pdb_index_map(aw_atoms, pdb_df, coord_decimals=3)
    aw_atoms = append_pdb_labels(aw_atoms, pdb_df, pdb_index_map)
    pow_atoms = append_pdb_labels(pow_atoms, pdb_df, pdb_index_map)

    debug_patch_neighbors(
        patch_indices=patch_indices,
        aw_atoms=aw_atoms,
        pow_atoms=pow_atoms,
        max_patch_atoms=5,
        max_neighbors_per_patch=20,
    )

    mapping_diag.to_csv(os.path.join(folder, 'pdb_mapping_diagnostics.csv'), index=False)

    mapping_summary = (
        mapping_diag.groupby('Match Method')
        .size()
        .reset_index(name='Count')
        .sort_values('Count', ascending=False)
    )
    mapping_summary.to_csv(os.path.join(folder, 'pdb_mapping_summary.csv'), index=False)

    matched_count = int(mapping_diag['PDB Row Index'].notna().sum())
    total_count = int(len(mapping_diag))
    print('\n=== PDB MAPPING SUMMARY ===')
    print(mapping_summary.to_string(index=False))
    print(f'\nMatched {matched_count} of {total_count} AW atoms to PDB rows ({matched_count / total_count:.2%}).')

    missing_aw = sorted(set(patch_indices) - set(aw_atoms['Index'].tolist()))
    missing_pow = sorted(set(patch_indices) - set(pow_atoms['Index'].tolist()))

    print('\n=== PATCH INDEX CHECK ===')
    print(f'Patch atoms requested: {len(patch_indices)}')
    print(f'Missing in AW: {len(missing_aw)}')
    print(f'Missing in Pow: {len(missing_pow)}')

    if missing_aw:
        print('First missing AW indices:', missing_aw[:10])
    if missing_pow:
        print('First missing Pow indices:', missing_pow[:10])

    valid_patch_indices = [idx for idx in patch_indices if idx in set(aw_atoms['Index']) and idx in set(pow_atoms['Index'])]

    aw_surface_lookup = build_surface_lookup(aw_surfs)
    pow_surface_lookup = build_surface_lookup(pow_surfs)

    aw_pairs = collect_neighbor_pair_rows(
        valid_patch_indices,
        'AW',
        aw_atoms,
        aw_surface_lookup,
        pdb_df,
    )

    print_missing_neighbor_resolution(aw_pairs, max_rows=50)

    pow_pairs = collect_neighbor_pair_rows(
        valid_patch_indices,
        'Pow',
        pow_atoms,
        pow_surface_lookup,
        pdb_df,
    )

    pair_compare = compare_aw_pow_pairs(aw_pairs, pow_pairs)
    patch_summary = summarize_patch_atoms(valid_patch_indices, aw_atoms, pow_atoms, aw_pairs, pow_pairs)

    feature_cols = [
        'AW Surface Area', 'Pow Surface Area', 'Delta Surface Area',
        'AW Neighbor Count', 'Pow Neighbor Count', 'Delta Neighbor Count',
        'AW Neighbor Surface Area Sum', 'Pow Neighbor Surface Area Sum',
        'AW Shared Surface Area Sum', 'Pow Shared Surface Area Sum',
        'AW Shared Contact Area Sum', 'Pow Shared Contact Area Sum',
        'AW SOL Neighbor Count', 'Pow SOL Neighbor Count',
        'AW Non-SOL Neighbor Count', 'Pow Non-SOL Neighbor Count',
        'AW SOL Shared Surface Area Sum', 'Pow SOL Shared Surface Area Sum',
        'AW Non-SOL Shared Surface Area Sum', 'Pow Non-SOL Shared Surface Area Sum',
    ]
    feature_corr = compute_correlations(patch_summary, 'Delta Volume Abs', feature_cols)

    pair_feature_cols = [
        'AW Shared Surface Count', 'Pow Shared Surface Count', 'Delta Shared Surface Count',
        'AW Shared Surface Area Total', 'Pow Shared Surface Area Total', 'Delta Shared Surface Area Total',
        'AW Shared Contact Area Total', 'Pow Shared Contact Area Total', 'Delta Shared Contact Area Total',
        'AW Shared Overlap Total', 'Pow Shared Overlap Total', 'Delta Shared Overlap Total',
        'AW Shared Mean Curvature Mean', 'Pow Shared Mean Curvature Mean', 'Delta Shared Mean Curvature Mean',
        'AW Shared Mean Curvature MaxAbs', 'Pow Shared Mean Curvature MaxAbs', 'Delta Shared Mean Curvature MaxAbs',
    ]

    pair_join = pair_compare.merge(
        patch_summary[['Index', 'Delta Volume', 'Delta Volume Abs', 'Is SOL']],
        left_on='Patch Index',
        right_on='Index',
        how='left'
    )
    pair_feature_corr = compute_correlations(pair_join, 'Delta Volume Abs', pair_feature_cols)

    outlier_by_atom = summarize_outliers_by_group(
        patch_summary.rename(columns={'Atom Name': 'Atom Name Group'}),
        'Atom Name Group'
    )

    if len(aw_pairs):
        aw_neighbor_enrichment = (
            aw_pairs.groupby(['Neighbor Type', 'Neighbor PDB Residue'])
            .agg(
                Count=('Neighbor Index', 'size'),
                Mean_Shared_Surface_Area=('Shared Surface Area Total', 'mean'),
                Mean_Shared_Contact_Area=('Shared Contact Area Total', 'mean'),
            )
            .reset_index()
            .sort_values('Count', ascending=False)
        )
    else:
        aw_neighbor_enrichment = pd.DataFrame(
            columns=['Neighbor Type', 'Neighbor PDB Residue', 'Count', 'Mean_Shared_Surface_Area', 'Mean_Shared_Contact_Area']
        )

    patch_summary.to_csv(os.path.join(folder, 'patch_atom_volume_neighbor_summary.csv'), index=False)
    aw_pairs.to_csv(os.path.join(folder, 'patch_atom_aw_neighbor_pairs.csv'), index=False)
    pow_pairs.to_csv(os.path.join(folder, 'patch_atom_pow_neighbor_pairs.csv'), index=False)
    pair_compare.to_csv(os.path.join(folder, 'patch_atom_aw_pow_pair_compare.csv'), index=False)
    feature_corr.to_csv(os.path.join(folder, 'patch_atom_feature_correlations.csv'), index=False)
    pair_feature_corr.to_csv(os.path.join(folder, 'patch_atom_pair_feature_correlations.csv'), index=False)
    outlier_by_atom.to_csv(os.path.join(folder, 'patch_atom_outlier_summary_by_atom_name.csv'), index=False)
    aw_neighbor_enrichment.to_csv(os.path.join(folder, 'patch_atom_outlier_summary_by_neighbor_type.csv'), index=False)

    print('\n=== OUTPUTS WRITTEN ===')
    for name in [
        'pdb_mapping_diagnostics.csv',
        'pdb_mapping_summary.csv',
        'patch_atom_volume_neighbor_summary.csv',
        'patch_atom_aw_neighbor_pairs.csv',
        'patch_atom_pow_neighbor_pairs.csv',
        'patch_atom_aw_pow_pair_compare.csv',
        'patch_atom_feature_correlations.csv',
        'patch_atom_pair_feature_correlations.csv',
        'patch_atom_outlier_summary_by_atom_name.csv',
        'patch_atom_outlier_summary_by_neighbor_type.csv',
    ]:
        print(os.path.join(folder, name))


if __name__ == '__main__':
    main()
