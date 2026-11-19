import csv
import os
import sys
import traceback
import tkinter as tk
from tkinter import filedialog

import numpy as np


# =============================================================================
# Settings
# =============================================================================

SCHEME_ORDER = [
    ('1_Atom', 'Atom'),
    ('2_Encap', 'Encap'),
    ('3_Encap_SR', 'Encap_SR'),
    ('4_AD', 'AD'),
    ('5_AD_SR', 'AD_SR'),
    ('6_AD_MW', 'AD_MW'),
    ('7_AD_MW_SR', 'AD_MW_SR'),
]

VORONOI_TYPES = ['aw', 'pow']

REFERENCE_SCHEME = '1_Atom'
REFERENCE_VORONOI = 'aw'

RAW_DIAGNOSTIC_FILENAME = 'residue_instance_raw_volume_diagnostic_read_logs2.csv'
PERCENT_DIAGNOSTIC_FILENAME = 'residue_instance_percent_difference_diagnostic_read_logs2.csv'
SUMMARY_FILENAME = 'residue_percent_difference_summary_read_logs2.csv'

PRINT_SAMPLE_RESIDUES = True
SAMPLE_RESIDUE_COUNT = 20


# =============================================================================
# Import helpers
# =============================================================================

def add_vorpy_root_to_path():
    current_dir = os.path.abspath(os.path.dirname(__file__))

    search_dir = current_dir
    for _ in range(12):
        if os.path.isdir(os.path.join(search_dir, 'vorpy')):
            if search_dir not in sys.path:
                sys.path.append(search_dir)
            return search_dir

        parent = os.path.dirname(search_dir)
        if parent == search_dir:
            break
        search_dir = parent

    return None


def import_read_logs2():
    add_vorpy_root_to_path()

    import_attempts = [
        'vorpy.src.analyze.tools.compare.read_logs2',
        'vorpy.src.analyze.tools.compare.compare_files',
        'vorpy.src.analyze.tools.compare.read_logs',
    ]

    last_error = None
    for module_name in import_attempts:
        try:
            module = __import__(module_name, fromlist=['read_logs2'])
            if hasattr(module, 'read_logs2'):
                print(f'Using read_logs2 from: {module_name}')
                return module.read_logs2
        except Exception as exc:
            last_error = exc

    raise ImportError(f'Could not import read_logs2. Last error: {last_error}')


read_logs2 = import_read_logs2()


# =============================================================================
# Column handling
# =============================================================================

def normalize_col_name(name):
    return str(name).strip().lower().replace('_', '').replace('-', '').replace(' ', '')


def find_column(df, candidates, label, required=True):
    exact_lookup = {str(col): col for col in df.columns}
    normalized_lookup = {normalize_col_name(col): col for col in df.columns}

    for candidate in candidates:
        if candidate in exact_lookup:
            return exact_lookup[candidate]

        normalized = normalize_col_name(candidate)
        if normalized in normalized_lookup:
            return normalized_lookup[normalized]

    if required:
        print('\nAvailable atoms dataframe columns:')
        for col in df.columns:
            print(f'  - {col}')
        raise KeyError(f'Could not find required {label}. Tried: {candidates}')

    return None


def get_schema_columns(atoms_df):
    residue_name_col = find_column(
        atoms_df,
        [
            'Residue',
            'Residue Name',
            'ResidueName',
            'residue',
            'resname',
            'ResName',
            'res_name',
        ],
        'residue name column'
    )

    residue_id_col = find_column(
        atoms_df,
        [
            'Residue Sequence',
            'ResidueSequence',
            'residue_sequence',
            'Residue Number',
            'ResidueNumber',
            'residue_number',
            'resid',
            'ResID',
            'Residue ID',
            'ResidueID',
            'res_id',
            'resSeq',
            'ResSeq',
            'Sequence Number',
            'SequenceNumber',
        ],
        'residue ID/number column'
    )

    chain_col = find_column(
        atoms_df,
        [
            'Chain',
            'chain',
            'Chain ID',
            'ChainID',
            'chain_id',
        ],
        'chain column',
        required=False
    )

    volume_col = find_column(
        atoms_df,
        [
            'Volume',
            'volume',
            'Vol',
            'vol',
        ],
        'volume column'
    )

    surface_area_col = find_column(
        atoms_df,
        [
            'Surface Area',
            'SurfaceArea',
            'surface_area',
            'SA',
            'sa',
        ],
        'surface area column'
    )

    return {
        'residue_name': residue_name_col,
        'residue_id': residue_id_col,
        'chain': chain_col,
        'volume': volume_col,
        'surface_area': surface_area_col,
    }


# =============================================================================
# Data loading and aggregation
# =============================================================================

def get_atoms_df_from_logs(log_path):
    log_data = read_logs2(log_path)

    if isinstance(log_data, dict):
        if 'atoms' in log_data:
            return log_data['atoms']

        if 'Atoms' in log_data:
            return log_data['Atoms']

    raise ValueError(f'read_logs2 did not return an atoms dataframe for: {log_path}')


def clean_value(value):
    if value is None:
        return ''

    text = str(value).strip()
    if text.lower() in {'nan', 'none'}:
        return ''

    return text


def make_residue_key(row, cols):
    resname = clean_value(row[cols['residue_name']])
    resid = clean_value(row[cols['residue_id']])

    if cols['chain'] is None:
        chain = ''
    else:
        chain = clean_value(row[cols['chain']])

    if chain:
        key = f'{resname}_{chain}_{resid}'
        label = f'{resname} {chain}:{resid}'
    else:
        key = f'{resname}_{resid}'
        label = f'{resname} {resid}'

    return key, label, resname, resid, chain


def aggregate_specific_residues(atoms_df, scheme_name, vor_type):
    cols = get_schema_columns(atoms_df)
    residue_data = {}

    for _, row in atoms_df.iterrows():
        key, label, resname, resid, chain = make_residue_key(row, cols)

        if not resname or not resid:
            continue

        try:
            volume = float(row[cols['volume']])
        except (TypeError, ValueError):
            volume = 0.0

        try:
            surface_area = float(row[cols['surface_area']])
        except (TypeError, ValueError):
            surface_area = 0.0

        if key not in residue_data:
            residue_data[key] = {
                'label': label,
                'residue': resname,
                'residue_sequence': resid,
                'chain': chain,
                'volume': 0.0,
                'surface_area': 0.0,
                'atom_count': 0,
                'scheme': scheme_name,
                'voronoi': vor_type,
            }

        residue_data[key]['volume'] += volume
        residue_data[key]['surface_area'] += surface_area
        residue_data[key]['atom_count'] += 1

    return residue_data


def load_all_residue_data(model_folder):
    all_data = {}

    for scheme_name, _ in SCHEME_ORDER:
        scheme_path = os.path.join(model_folder, scheme_name)

        if not os.path.isdir(scheme_path):
            print(f'MISSING scheme folder: {scheme_path}')
            continue

        all_data[scheme_name] = {}

        for vor_type in VORONOI_TYPES:
            log_path = os.path.join(scheme_path, vor_type, f'{vor_type}_logs.csv')

            if not os.path.exists(log_path):
                print(f'MISSING log file: {log_path}')
                continue

            print(f'Reading: {log_path}')
            atoms_df = get_atoms_df_from_logs(log_path)
            all_data[scheme_name][vor_type] = aggregate_specific_residues(
                atoms_df=atoms_df,
                scheme_name=scheme_name,
                vor_type=vor_type
            )

            print(f'  residues found: {len(all_data[scheme_name][vor_type])}')

    return all_data


# =============================================================================
# Diagnostics
# =============================================================================

def get_all_residue_keys(all_data):
    keys = set()

    for scheme_name in all_data:
        for vor_type in all_data[scheme_name]:
            keys.update(all_data[scheme_name][vor_type].keys())

    return sorted(keys)


def get_reference_residue_info(all_data, key):
    for scheme_name, _ in SCHEME_ORDER:
        for vor_type in VORONOI_TYPES:
            try:
                return all_data[scheme_name][vor_type][key]
            except KeyError:
                continue

    return {
        'label': key,
        'residue': '',
        'residue_sequence': '',
        'chain': '',
        'atom_count': '',
    }


def write_raw_volume_diagnostic(all_data, output_path):
    header = ['Residue_Key', 'Residue_Label', 'Residue', 'Chain', 'Residue_Sequence']

    for scheme_name, scheme_label in SCHEME_ORDER:
        for vor_type in VORONOI_TYPES:
            header.append(f'{scheme_label}_{vor_type}_volume')

    keys = get_all_residue_keys(all_data)

    with open(output_path, 'w', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(header)

        for key in keys:
            info = get_reference_residue_info(all_data, key)
            row = [
                key,
                info.get('label', key),
                info.get('residue', ''),
                info.get('chain', ''),
                info.get('residue_sequence', ''),
            ]

            for scheme_name, _ in SCHEME_ORDER:
                for vor_type in VORONOI_TYPES:
                    value = ''
                    try:
                        value = all_data[scheme_name][vor_type][key]['volume']
                    except KeyError:
                        pass
                    row.append(value)

            writer.writerow(row)

    print(f'Wrote raw volume diagnostic: {output_path}')


def write_percent_difference_diagnostic(all_data, output_path):
    ref_data = all_data.get(REFERENCE_SCHEME, {}).get(REFERENCE_VORONOI, {})

    header = ['Residue_Key', 'Residue_Label', 'Residue', 'Chain', 'Residue_Sequence', 'Reference_Atom_aw_volume']

    for scheme_name, scheme_label in SCHEME_ORDER:
        for vor_type in VORONOI_TYPES:
            header.append(f'{scheme_label}_{vor_type}_abs_percent_diff_volume')

    keys = sorted(ref_data.keys())

    with open(output_path, 'w', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(header)

        for key in keys:
            ref = ref_data[key]['volume']
            info = ref_data[key]

            row = [
                key,
                info.get('label', key),
                info.get('residue', ''),
                info.get('chain', ''),
                info.get('residue_sequence', ''),
                ref,
            ]

            for scheme_name, _ in SCHEME_ORDER:
                for vor_type in VORONOI_TYPES:
                    value = ''
                    try:
                        comp = all_data[scheme_name][vor_type][key]['volume']
                        if ref != 0:
                            value = abs(comp - ref) / ref * 100.0
                    except KeyError:
                        pass
                    row.append(value)

            writer.writerow(row)

    print(f'Wrote percent-difference diagnostic: {output_path}')


def print_sample_matches(all_data):
    if not PRINT_SAMPLE_RESIDUES:
        return

    ref_data = all_data.get(REFERENCE_SCHEME, {}).get(REFERENCE_VORONOI, {})
    keys = sorted(ref_data.keys())[:SAMPLE_RESIDUE_COUNT]

    print('\n=== Sample residue instance keys from reference Atom/aw ===')
    for key in keys:
        info = ref_data[key]
        print(
            f"{key} | label={info['label']} | "
            f"volume={info['volume']:.3f} | SA={info['surface_area']:.3f} | atoms={info['atom_count']}"
        )


# =============================================================================
# Summary calculation
# =============================================================================

def compute_summary_rows(all_data):
    ref_data = all_data.get(REFERENCE_SCHEME, {}).get(REFERENCE_VORONOI, {})

    if not ref_data:
        raise ValueError('Reference data was not found: 1_Atom/aw')

    rows = []

    for scheme_name, scheme_label in SCHEME_ORDER:
        for vor_type in VORONOI_TYPES:
            comp_data = all_data.get(scheme_name, {}).get(vor_type, {})
            diffs = []

            for key, ref_info in ref_data.items():
                if key not in comp_data:
                    continue

                ref_volume = ref_info['volume']
                comp_volume = comp_data[key]['volume']

                if ref_volume == 0:
                    continue

                diffs.append(abs(comp_volume - ref_volume) / ref_volume * 100.0)

            if len(diffs) == 0:
                mean_diff = ''
                stderr_diff = ''
                n = 0
            else:
                diff_array = np.array(diffs, dtype=float)
                mean_diff = float(np.mean(diff_array))
                stderr_diff = float(np.std(diff_array) / np.sqrt(len(diff_array)))
                n = len(diff_array)

            rows.append({
                'scheme': scheme_name,
                'scheme_label': scheme_label,
                'voronoi': vor_type,
                'mean_abs_percent_diff_volume': mean_diff,
                'stderr_abs_percent_diff_volume': stderr_diff,
                'n_matched_residues': n,
            })

    return rows


def write_summary(rows, output_path):
    header = [
        'scheme',
        'scheme_label',
        'voronoi',
        'mean_abs_percent_diff_volume',
        'stderr_abs_percent_diff_volume',
        'n_matched_residues',
    ]

    with open(output_path, 'w', newline='') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f'Wrote summary: {output_path}')


def print_summary(rows):
    print('\n=== Mean absolute percent difference from Atom/aw reference ===')
    for row in rows:
        mean_diff = row['mean_abs_percent_diff_volume']
        stderr = row['stderr_abs_percent_diff_volume']

        if mean_diff == '':
            mean_text = 'NA'
            stderr_text = 'NA'
        else:
            mean_text = f'{mean_diff:.3f}'
            stderr_text = f'{stderr:.3f}'

        print(
            f"{row['scheme_label']:>8s} {row['voronoi']:>3s}: "
            f"mean={mean_text}% | stderr={stderr_text}% | n={row['n_matched_residues']}"
        )


# =============================================================================
# Main
# =============================================================================

def select_model_folder():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory(title='Select model folder, e.g. K_NCP')
    root.destroy()

    if not folder:
        raise ValueError('No folder selected.')

    return folder


def main():
    model_folder = select_model_folder()
    print(f'\nSelected model folder: {model_folder}')

    all_data = load_all_residue_data(model_folder)
    print_sample_matches(all_data)

    raw_output_path = os.path.join(model_folder, RAW_DIAGNOSTIC_FILENAME)
    percent_output_path = os.path.join(model_folder, PERCENT_DIAGNOSTIC_FILENAME)
    summary_output_path = os.path.join(model_folder, SUMMARY_FILENAME)

    write_raw_volume_diagnostic(all_data, raw_output_path)
    write_percent_difference_diagnostic(all_data, percent_output_path)

    rows = compute_summary_rows(all_data)
    write_summary(rows, summary_output_path)
    print_summary(rows)

    print('\nDone.')


if __name__ == '__main__':
    try:
        main()
    except Exception:
        print('\nSCRIPT FAILED. FULL TRACEBACK:\n')
        traceback.print_exc()
        input('\nPress Enter to close...')
