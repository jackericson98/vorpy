import ast
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd



DEFAULT_SOLVENT_RESIDUES: Set[str] = {
    'SOL',
    'HOH',
    'WAT',
    'TIP3',
    'TIP3P',
    'SPC',
    'SPCE',
    'H2O',
}



def _find_first_matching_column(
    df: pd.DataFrame,
    candidates: List[str],
    df_name: str
) -> str:
    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        f"Could not find any of the expected columns in {df_name}: {candidates}"
    )



def _coerce_index_value(value) -> Optional[int]:
    if pd.isna(value):
        return None

    try:
        return int(value)

    except Exception:
        try:
            return int(float(value))

        except Exception:
            return None



def _normalize_residue_name(value) -> str:
    if pd.isna(value):
        return ''

    return str(value).strip().upper()



def _parse_neighbor_list(value) -> List[int]:
    if pd.isna(value):
        return []

    if isinstance(value, (list, tuple, np.ndarray)):
        out = []

        for item in value:
            idx = _coerce_index_value(item)

            if idx is not None:
                out.append(idx)

        return out

    text = str(value).strip()

    if not text:
        return []

    try:
        parsed = ast.literal_eval(text)

        if isinstance(parsed, (list, tuple, np.ndarray)):
            out = []

            for item in parsed:
                idx = _coerce_index_value(item)

                if idx is not None:
                    out.append(idx)

            return out

    except Exception:
        pass

    text = text.replace('[', '').replace(']', '').replace('(', '').replace(')', '')
    parts = [part.strip() for part in text.split(',') if part.strip()]

    out = []

    for part in parts:
        idx = _coerce_index_value(part)

        if idx is not None:
            out.append(idx)

    return out



def _build_solvent_lookup(
    atom_df: pd.DataFrame,
    residue_col: str,
    solvent_residues: Set[str]
) -> Dict[int, bool]:
    index_col = _find_first_matching_column(
        atom_df,
        ['Index', 'Atom Index', 'AtomIndex'],
        'atom_df'
    )

    solvent_lookup: Dict[int, bool] = {}

    for _, row in atom_df.iterrows():
        atom_idx = _coerce_index_value(row[index_col])

        if atom_idx is None:
            continue

        residue_name = _normalize_residue_name(row[residue_col])
        solvent_lookup[atom_idx] = residue_name in solvent_residues

    return solvent_lookup



def _resolve_surface_columns(surfs_df: pd.DataFrame) -> Tuple[str, str, str]:
    atom1_col = _find_first_matching_column(
        surfs_df,
        [
            'Ball 1',
            'Ball1',
            'Atom 1',
            'Atom1',
            'Index 1',
            'Index1',
            'i',
            'a',
        ],
        'surfs_df'
    )

    atom2_col = _find_first_matching_column(
        surfs_df,
        [
            'Ball 2',
            'Ball2',
            'Atom 2',
            'Atom2',
            'Index 2',
            'Index2',
            'j',
            'b',
        ],
        'surfs_df'
    )

    area_col = _find_first_matching_column(
        surfs_df,
        [
            'Surface Area',
            'SurfaceArea',
            'Area',
            'SA',
            'Surface_Area',
        ],
        'surfs_df'
    )

    return atom1_col, atom2_col, area_col



def _build_surface_index(
    surfs_df: pd.DataFrame,
    atom1_col: str,
    atom2_col: str,
    area_col: str
) -> Dict[int, List[Tuple[int, float]]]:
    surface_index: Dict[int, List[Tuple[int, float]]] = {}

    for _, row in surfs_df.iterrows():
        atom1 = _coerce_index_value(row[atom1_col])
        atom2 = _coerce_index_value(row[atom2_col])

        if atom1 is None or atom2 is None:
            continue

        try:
            area = float(row[area_col])

        except Exception:
            continue

        if not np.isfinite(area):
            continue

        surface_index.setdefault(atom1, []).append((atom2, area))
        surface_index.setdefault(atom2, []).append((atom1, area))

    return surface_index



def _get_neighbor_map_from_atom_df(
    atom_df: pd.DataFrame
) -> Optional[Dict[int, List[int]]]:
    index_candidates = ['Index', 'Atom Index', 'AtomIndex']
    neighbor_candidates = [
        'Neighbors',
        'Neighbor Indices',
        'NeighborIndices',
        'Adjacent Balls',
        'AdjacentBallIndices',
        'Touching Neighbors',
    ]

    index_col = None
    neighbor_col = None

    for col in index_candidates:
        if col in atom_df.columns:
            index_col = col
            break

    for col in neighbor_candidates:
        if col in atom_df.columns:
            neighbor_col = col
            break

    if index_col is None or neighbor_col is None:
        return None

    neighbor_map: Dict[int, List[int]] = {}

    for _, row in atom_df.iterrows():
        atom_idx = _coerce_index_value(row[index_col])

        if atom_idx is None:
            continue

        neighbor_map[atom_idx] = _parse_neighbor_list(row[neighbor_col])

    return neighbor_map



def compute_sol_facing_percent(
    atom_df: pd.DataFrame,
    surfs_df: pd.DataFrame,
    solvent_residues: Optional[Set[str]] = None,
    include_solvent_atoms: bool = False,
    total_sa_col: Optional[str] = None,
    use_atom_neighbor_column: bool = False
) -> pd.DataFrame:
    """
    Compute solvent-facing surface area and solvent-facing percentage for each atom.

    Strategy:
    - Determine which atoms are solvent by residue name.
    - Find all surfaces touching each atom.
    - Sum total touching surface area.
    - Sum only the touching surface area where the partner atom is solvent.
    - Compute percentage = solvent-facing area / total area * 100.

    Parameters
    ----------
    atom_df : pd.DataFrame
        Atom table. Must contain atom index and residue name columns.

    surfs_df : pd.DataFrame
        Surface table. Must contain the two atom indices defining each surface and the surface area.

    solvent_residues : Optional[Set[str]]
        Residue names that count as solvent. Defaults to common water names.

    include_solvent_atoms : bool
        If False, solvent atoms are skipped from the output.
        If True, solvent atoms are included too.

    total_sa_col : Optional[str]
        Optional atom_df column to use as the denominator for total surface area.
        If None, denominator is computed directly from surfs_df by summing all touching surfaces.

    use_atom_neighbor_column : bool
        If True, attempts to read a neighbor list column from atom_df and restrict
        solvent-facing surfaces to those neighbor relationships. Usually not necessary,
        because surfaces already define the touching pairs.

    Returns
    -------
    pd.DataFrame
        Copy of atom_df with added columns:
        - IsSolventAtom
        - TotalTouchingSurfaceArea
        - SolFacingSurfaceArea
        - SolFacingPct
        - SolFacingNeighborCount
        - SolFacingNeighborIndices
    """
    if solvent_residues is None:
        solvent_residues = DEFAULT_SOLVENT_RESIDUES

    solvent_residues = {str(x).strip().upper() for x in solvent_residues}

    out_df = atom_df.copy()

    index_col = _find_first_matching_column(
        out_df,
        ['Index', 'Atom Index', 'AtomIndex'],
        'atom_df'
    )

    residue_col = _find_first_matching_column(
        out_df,
        ['Residue', 'ResidueName', 'Residue Name', 'ResName'],
        'atom_df'
    )

    atom1_col, atom2_col, area_col = _resolve_surface_columns(surfs_df)

    solvent_lookup = _build_solvent_lookup(
        atom_df=out_df,
        residue_col=residue_col,
        solvent_residues=solvent_residues
    )

    surface_index = _build_surface_index(
        surfs_df=surfs_df,
        atom1_col=atom1_col,
        atom2_col=atom2_col,
        area_col=area_col
    )

    neighbor_map = None

    if use_atom_neighbor_column:
        neighbor_map = _get_neighbor_map_from_atom_df(out_df)

    is_solvent_values = []
    total_area_values = []
    sol_area_values = []
    sol_pct_values = []
    sol_neighbor_count_values = []
    sol_neighbor_indices_values = []

    for _, row in out_df.iterrows():
        atom_idx = _coerce_index_value(row[index_col])

        if atom_idx is None:
            is_solvent_values.append(False)
            total_area_values.append(np.nan)
            sol_area_values.append(np.nan)
            sol_pct_values.append(np.nan)
            sol_neighbor_count_values.append(np.nan)
            sol_neighbor_indices_values.append('')
            continue

        is_solvent_atom = solvent_lookup.get(atom_idx, False)

        if (not include_solvent_atoms) and is_solvent_atom:
            is_solvent_values.append(True)
            total_area_values.append(np.nan)
            sol_area_values.append(np.nan)
            sol_pct_values.append(np.nan)
            sol_neighbor_count_values.append(np.nan)
            sol_neighbor_indices_values.append('')
            continue

        touching_surfaces = surface_index.get(atom_idx, [])

        if neighbor_map is not None:
            allowed_neighbors = set(neighbor_map.get(atom_idx, []))
            touching_surfaces = [
                (nbr_idx, area)
                for nbr_idx, area in touching_surfaces
                if nbr_idx in allowed_neighbors
            ]

        total_touching_area = 0.0
        sol_facing_area = 0.0
        sol_neighbor_indices = set()

        for nbr_idx, area in touching_surfaces:
            if not np.isfinite(area):
                continue

            total_touching_area += area

            if solvent_lookup.get(nbr_idx, False):
                sol_facing_area += area
                sol_neighbor_indices.add(nbr_idx)

        if total_sa_col is not None and total_sa_col in out_df.columns:
            try:
                denominator = float(row[total_sa_col])

            except Exception:
                denominator = total_touching_area

        else:
            denominator = total_touching_area

        if denominator > 0:
            sol_pct = 100.0 * sol_facing_area / denominator

        else:
            sol_pct = 0.0

        is_solvent_values.append(is_solvent_atom)
        total_area_values.append(total_touching_area)
        sol_area_values.append(sol_facing_area)
        sol_pct_values.append(sol_pct)
        sol_neighbor_count_values.append(len(sol_neighbor_indices))
        sol_neighbor_indices_values.append(
            ','.join(str(idx) for idx in sorted(sol_neighbor_indices))
        )

    out_df['IsSolventAtom'] = is_solvent_values
    out_df['TotalTouchingSurfaceArea'] = total_area_values
    out_df['SolFacingSurfaceArea'] = sol_area_values
    out_df['SolFacingPct'] = sol_pct_values
    out_df['SolFacingNeighborCount'] = sol_neighbor_count_values
    out_df['SolFacingNeighborIndices'] = sol_neighbor_indices_values

    return out_df