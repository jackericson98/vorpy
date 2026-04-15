import os
import sys
import json
import tkinter as tk
from tkinter import filedialog
import pandas as pd

import pandas as pd
import matplotlib.pyplot as plt


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


RESIDUE_CLASS_LABEL = {
    # Hydrophobic
    'ALA': 'H', 'VAL': 'H', 'LEU': 'H', 'ILE': 'H',
    'MET': 'H', 'PRO': 'H',

    # Aromatic
    'PHE': 'A', 'TYR': 'A', 'TRP': 'A',

    # Polar
    'SER': 'P', 'THR': 'P', 'ASN': 'P',
    'GLN': 'P', 'CYS': 'P',

    # Positive
    'LYS': '+', 'ARG': '+', 'HIS': '+',

    # Negative
    'ASP': '-', 'GLU': '-',

    # Special
    'GLY': 'S'
}


CLASS_MARKER_MAP = {
    'H': 'o',   # circle
    'A': 's',   # square
    'P': '^',   # triangle up
    '+': 'D',   # diamond
    '-': 'v',   # triangle down
    'S': 'X',   # x-filled
    '?': 'o'
}

RESIDUE_COLOR_MAP = {
    # Hydrophobic
    'ALA': '#1f77b4',   # blue
    'VAL': '#ff7f0e',   # orange
    'LEU': '#2ca02c',   # green
    'ILE': '#d62728',   # red
    'MET': '#9467bd',   # purple
    'PRO': '#8c564b',   # brown

    # Aromatic
    'PHE': '#e377c2',   # pink
    'TYR': '#7f7f7f',   # gray
    'TRP': '#bcbd22',   # olive

    # Polar
    'SER': '#17becf',   # cyan
    'THR': '#393b79',   # dark blue
    'ASN': '#637939',   # olive green
    'GLN': '#8c6d31',   # mustard brown
    'CYS': '#843c39',   # dark brick

    # Positive
    'LYS': '#e41a1c',   # strong red
    'ARG': '#377eb8',   # strong blue
    'HIS': '#4daf4a',   # strong green

    # Negative
    'ASP': '#984ea3',   # strong violet
    'GLU': '#ff1493',   # deep pink

    # Special
    'GLY': '#000000'    # black
}


PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

DNA_RNA_RESIDUES = {
    'A', 'C', 'G', 'T', 'U',
    'DA', 'DC', 'DG', 'DT', 'DU',
    'ADE', 'CYT', 'GUA', 'THY', 'URA'
}

DEFAULT_X_RANGE = (3, 22)
DEFAULT_Y_RANGE = (3, 22)

SAVE_PLOTS = True
SHOW_PLOTS = True
PROMPT_FOR_GROUPING = True
SKIP_ALREADY_REVIEWED = True


def pick_folders():

    folders = []

    while True:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory(title='Pick a system folder (Cancel when done)')
        root.destroy()

        if folder == '' or folder is None:
            break

        folders.append(folder)

    return folders


def get_output_directory():

    root = tk.Tk()
    root.withdraw()
    out_dir = filedialog.askdirectory(title='Pick output folder for plots and grouping CSV')
    root.destroy()

    return out_dir


def safe_system_name(folder):

    return os.path.basename(folder.rstrip('/\\'))


def load_logs_from_folder(folder):

    try:
        aw_logs = read_logs2(os.path.join(folder, 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm_logs.csv'), all_=False, balls=True)
    except FileNotFoundError:
        aw_logs = read_logs2(os.path.join(folder, 'aw', 'aw_logs.csv'), all_=False, balls=True)
        pow_logs = read_logs2(os.path.join(folder, 'pow', 'pow_logs.csv'), all_=False, balls=True)
        prm_logs = read_logs2(os.path.join(folder, 'prm', 'prm_logs.csv'), all_=False, balls=True)

    return aw_logs, pow_logs, prm_logs


def parse_group_string(note):
    """
    Parse grouping notes where:
      - commas separate distinct groups
      - slashes join residues within a group

    Example:
      'LEU, PHE/TYR, LYS/ARG/MET/GLU/PRO/GLN, HIS, ASN, ASP'
    ->
      [['LEU'], ['PHE', 'TYR'], ['LYS', 'ARG', 'MET', 'GLU', 'PRO', 'GLN'], ['HIS'], ['ASN'], ['ASP']]
    """
    if pd.isna(note):
        return []

    note = str(note).strip()
    if note == '':
        return []

    raw_groups = [x.strip() for x in note.split(',') if x.strip()]

    parsed_groups = []
    for group in raw_groups:
        residues = [r.strip().upper() for r in group.split('/') if r.strip()]
        if residues:
            parsed_groups.append(residues)

    return parsed_groups


def build_atom_group_mapping(progress_csv):
    df = pd.read_csv(progress_csv)

    mapping = {}

    for _, row in df.iterrows():
        atom_name = str(row['AtomName']).strip()
        note = row['Notes']

        mapping[atom_name] = parse_group_string(note)

    return mapping


def standardize_columns(df):


    df = df.copy()
    df.columns = [str(col).strip() for col in df.columns]

    rename_map = {}

    for col in df.columns:
        low = col.lower()

        if low == 'residue name':
            rename_map[col] = 'ResidueName'
        elif low == 'residuename':
            rename_map[col] = 'ResidueName'
        elif low == 'residue':
            rename_map[col] = 'ResidueName'
        elif low == 'name':
            rename_map[col] = 'AtomName'
        elif low == 'atom name':
            rename_map[col] = 'AtomName'
        elif low == 'index':
            rename_map[col] = 'Index'
        elif low == 'volume':
            rename_map[col] = 'Volume'
        elif low == 'maximum mean curvature':
            rename_map[col] = 'MaximumMeanCurvature'

    df = df.rename(columns=rename_map)

    return df


def ensure_required_columns(df, required, label='DataFrame'):

    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(f'{label} is missing required columns: {missing}')


def build_merged_atom_dataframe(folder):

    aw_logs, pow_logs, prm_logs = load_logs_from_folder(folder)

    aw_atoms = standardize_columns(aw_logs['atoms'])
    pow_atoms = standardize_columns(pow_logs['atoms'])
    prm_atoms = standardize_columns(prm_logs['atoms'])

    ensure_required_columns(
        aw_atoms,
        ['Index', 'AtomName', 'ResidueName', 'Volume'],
        label='AW atoms'
    )

    ensure_required_columns(
        pow_atoms,
        ['Index', 'Volume'],
        label='Pow atoms'
    )

    ensure_required_columns(
        prm_atoms,
        ['Index', 'Volume'],
        label='Prm atoms'
    )

    keep_aw_cols = ['Index', 'AtomName', 'ResidueName', 'Volume']
    if 'MaximumMeanCurvature' in aw_atoms.columns:
        keep_aw_cols.append('MaximumMeanCurvature')

    aw_atoms = aw_atoms[keep_aw_cols].copy()
    pow_atoms = pow_atoms[['Index', 'Volume']].copy()
    prm_atoms = prm_atoms[['Index', 'Volume']].copy()

    pow_atoms = pow_atoms.rename(columns={'Volume': 'PowVolume'})
    prm_atoms = prm_atoms.rename(columns={'Volume': 'PrmVolume'})
    aw_atoms = aw_atoms.rename(columns={'Volume': 'AWVolume'})

    merged = aw_atoms.merge(pow_atoms, on='Index', how='inner')
    merged = merged.merge(prm_atoms, on='Index', how='inner')

    merged['ResidueName'] = merged['ResidueName'].astype(str).str.strip().str.upper()
    merged['AtomName'] = merged['AtomName'].astype(str).str.strip()

    merged['System'] = safe_system_name(folder)

    return merged


def filter_protein_atoms(df):

    df = df.copy()

    protein_mask = df['ResidueName'].isin(PROTEIN_RESIDUES)
    df = df[protein_mask]

    df = df[
        df['AWVolume'].between(DEFAULT_X_RANGE[0], DEFAULT_X_RANGE[1]) &
        df['PowVolume'].between(DEFAULT_Y_RANGE[0], DEFAULT_Y_RANGE[1]) &
        df['PrmVolume'].between(DEFAULT_Y_RANGE[0], DEFAULT_Y_RANGE[1])
    ].copy()

    return df


def load_or_initialize_progress_csv(progress_csv, atom_names, counts_dict, residues_dict):

    if os.path.exists(progress_csv):
        progress_df = pd.read_csv(progress_csv)
    else:
        progress_df = pd.DataFrame({
            'AtomName': pd.Series(dtype='string'),
            'Count': pd.Series(dtype='Int64'),
            'ResiduesPresent': pd.Series(dtype='string'),
            'Reviewed': pd.Series(dtype='string'),
            'GroupLabel': pd.Series(dtype='string'),
            'Notes': pd.Series(dtype='string')
        })

    # ensure required columns exist
    for col in ['AtomName', 'Count', 'ResiduesPresent', 'Reviewed', 'GroupLabel', 'Notes']:
        if col not in progress_df.columns:
            if col == 'Count':
                progress_df[col] = pd.Series(dtype='Int64')
            else:
                progress_df[col] = pd.Series(dtype='string')

    # normalize dtypes
    progress_df['AtomName'] = progress_df['AtomName'].astype('string').str.strip()
    progress_df['ResiduesPresent'] = progress_df['ResiduesPresent'].astype('string')
    progress_df['Reviewed'] = progress_df['Reviewed'].astype('string')
    progress_df['GroupLabel'] = progress_df['GroupLabel'].astype('string')
    progress_df['Notes'] = progress_df['Notes'].astype('string')
    progress_df['Count'] = pd.to_numeric(progress_df['Count'], errors='coerce').astype('Int64')

    # drop blank / bad AtomName rows
    progress_df = progress_df.dropna(subset=['AtomName']).copy()
    progress_df = progress_df[progress_df['AtomName'].str.len() > 0].copy()

    existing_atoms = set(progress_df['AtomName'].tolist())

    rows_to_add = []

    for atom_name in atom_names:
        atom_name = str(atom_name).strip()

        if atom_name not in existing_atoms:
            rows_to_add.append({
                'AtomName': atom_name,
                'Count': counts_dict.get(atom_name, 0),
                'ResiduesPresent': ', '.join(residues_dict.get(atom_name, [])),
                'Reviewed': '',
                'GroupLabel': '',
                'Notes': ''
            })

    if rows_to_add:
        progress_df = pd.concat([progress_df, pd.DataFrame(rows_to_add)], ignore_index=True)

    # one more normalization pass after concat
    progress_df['AtomName'] = progress_df['AtomName'].astype('string').str.strip()
    progress_df['ResiduesPresent'] = progress_df['ResiduesPresent'].astype('string')
    progress_df['Reviewed'] = progress_df['Reviewed'].astype('string')
    progress_df['GroupLabel'] = progress_df['GroupLabel'].astype('string')
    progress_df['Notes'] = progress_df['Notes'].astype('string')
    progress_df['Count'] = pd.to_numeric(progress_df['Count'], errors='coerce').astype('Int64')

    progress_df = progress_df.drop_duplicates(subset=['AtomName'], keep='first')
    progress_df = progress_df.sort_values(by=['AtomName']).reset_index(drop=True)

    progress_df.to_csv(progress_csv, index=False)

    return progress_df


def update_progress_row(progress_df, progress_csv, atom_name, reviewed=None, group_label=None, notes=None):

    atom_name = str(atom_name).strip()

    match = progress_df.index[
        progress_df['AtomName'].astype('string').str.strip() == atom_name
    ]

    if len(match) == 0:
        new_row = {
            'AtomName': atom_name,
            'Count': pd.NA,
            'ResiduesPresent': '',
            'Reviewed': '',
            'GroupLabel': '',
            'Notes': ''
        }
        progress_df = pd.concat([progress_df, pd.DataFrame([new_row])], ignore_index=True)
        match = progress_df.index[
            progress_df['AtomName'].astype('string').str.strip() == atom_name
        ]

    idx = match[0]

    if reviewed is not None:
        progress_df.loc[idx, 'Reviewed'] = str(reviewed)

    if group_label is not None:
        progress_df.loc[idx, 'GroupLabel'] = str(group_label)

    if notes is not None:
        progress_df.loc[idx, 'Notes'] = str(notes)

    progress_df.to_csv(progress_csv, index=False)

    return progress_df


def make_atom_plot(atom_df, atom_name, out_path=None):

    fig, ax = plt.subplots(figsize=(8, 6))

    residues = sorted(
        atom_df['ResidueName'].unique(),
        key=lambda r: (RESIDUE_CLASS_LABEL.get(r, '?'), r)
    )

    for residue in residues:
        sub = atom_df[atom_df['ResidueName'] == residue]

        class_label = RESIDUE_CLASS_LABEL.get(residue, '?')
        marker_style = CLASS_MARKER_MAP.get(class_label, 'o')
        point_color = RESIDUE_COLOR_MAP.get(residue, '#000000')
        legend_label = f"{residue} ({class_label})"

        ax.scatter(
            sub['AWVolume'],
            sub['PowVolume'],
            s=95,
            alpha=0.75,
            marker=marker_style,
            c=point_color,
            linewidths=0.7,
            label=legend_label
        )

    ax.plot([DEFAULT_X_RANGE[0], DEFAULT_X_RANGE[1]],
            [DEFAULT_Y_RANGE[0], DEFAULT_Y_RANGE[1]],
            linestyle='--',
            linewidth=2.5,
            color='black',
            alpha=0.8)

    if len(residues) <= 20:
        ax.legend(
            fontsize=10,
            frameon=False,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0
        )

    ax.set_xlim(*DEFAULT_X_RANGE)
    ax.set_ylim(*DEFAULT_Y_RANGE)

    ax.set_xlabel('AW Volume', fontsize=20)
    ax.set_ylabel('Pow Volume', fontsize=20)

    title = f'{atom_name} | n={len(atom_df)} | residues={len(residues)}'
    ax.set_title(title, fontsize=18)

    ax.tick_params(axis='both', which='major', labelsize=16, width=2.5, length=8)

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    if len(residues) <= 12:
        ax.legend(fontsize=10, frameon=False, loc='best')

    plt.tight_layout()

    if out_path is not None:
        plt.savefig(out_path, dpi=300, bbox_inches='tight')

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)


def prompt_for_group(atom_name, atom_df):
    residues = sorted(
        atom_df['ResidueName'].unique(),
        key=lambda r: RESIDUE_CLASS_LABEL.get(r, '?')
    )

    print('\n' + '=' * 80)
    print(f'Atom Name: {atom_name}')
    print(f'Count: {len(atom_df)}')
    print(f'Residues Present: {", ".join(residues)}')
    print('=' * 80)

    group_label = input('Enter group label for this atom name (blank = skip): ').strip()
    notes = input('Enter notes (blank = none): ').strip()
    reviewed = input('Mark reviewed? [y/n, blank=y]: ').strip().lower()

    if reviewed == '':
        reviewed = 'y'

    reviewed = 'yes' if reviewed.startswith('y') else 'no'

    return group_label, notes, reviewed


def save_grouping_summary_json(progress_df, json_path):

    summary = {}

    for _, row in progress_df.iterrows():
        atom_name = row['AtomName']
        summary[str(atom_name)] = {
            'group_label': '' if pd.isna(row['GroupLabel']) else str(row['GroupLabel']),
            'reviewed': '' if pd.isna(row['Reviewed']) else str(row['Reviewed']),
            'notes': '' if pd.isna(row['Notes']) else str(row['Notes']),
            'count': 0 if pd.isna(row['Count']) else int(row['Count']),
            'residues_present': '' if pd.isna(row['ResiduesPresent']) else str(row['ResiduesPresent'])
        }

    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=4)


def main():

    folders = pick_folders()

    if not folders:
        print('No folders selected.')
        return

    out_dir = get_output_directory()

    if not out_dir:
        print('No output folder selected.')
        return

    all_frames = []

    for folder in folders:
        print(f'Loading: {folder}')
        merged = build_merged_atom_dataframe(folder)
        merged = filter_protein_atoms(merged)
        all_frames.append(merged)

    if not all_frames:
        print('No valid protein atom data found.')
        return

    df = pd.concat(all_frames, ignore_index=True)

    atom_names = sorted(df['AtomName'].dropna().unique())
    counts_dict = df.groupby('AtomName').size().to_dict()
    residues_dict = {
        atom_name: sorted(df.loc[df['AtomName'] == atom_name, 'ResidueName'].unique())
        for atom_name in atom_names
    }

    progress_csv = os.path.join(out_dir, 'protein_atom_grouping_progress.csv')
    progress_json = os.path.join(out_dir, 'protein_atom_grouping_progress.json')
    plot_dir = os.path.join(out_dir, 'atom_name_plots')

    os.makedirs(plot_dir, exist_ok=True)

    progress_df = pd.DataFrame({
        'AtomName': pd.Series(dtype='str'),
        'Count': pd.Series(dtype='int'),
        'ResiduesPresent': pd.Series(dtype='str'),
        'Reviewed': pd.Series(dtype='str'),
        'GroupLabel': pd.Series(dtype='str'),
        'Notes': pd.Series(dtype='str')
    })

    for atom_name in atom_names:
        atom_name = str(atom_name).strip()
        match = progress_df.loc[progress_df['AtomName'].astype('string').str.strip() == atom_name]

        if match.empty:
            # repair missing row on the fly
            new_row = pd.DataFrame([{
                'AtomName': atom_name,
                'Count': counts_dict.get(atom_name, 0),
                'ResiduesPresent': ', '.join(residues_dict.get(atom_name, [])),
                'Reviewed': '',
                'GroupLabel': '',
                'Notes': ''
            }])

            progress_df = pd.concat([progress_df, new_row], ignore_index=True)
            progress_df.to_csv(progress_csv, index=False)

            match = progress_df.loc[progress_df['AtomName'].astype('string').str.strip() == atom_name]

        row = match.iloc[0]

        if SKIP_ALREADY_REVIEWED and str(row['Reviewed']).strip().lower() == 'yes':
            continue

        atom_df = df[df['AtomName'] == atom_name].copy()

        if atom_df.empty:
            continue

        out_path = os.path.join(plot_dir, f'{atom_name}.png')

        make_atom_plot(atom_df, atom_name, out_path=out_path if SAVE_PLOTS else None)

        if PROMPT_FOR_GROUPING:
            group_label, notes, reviewed = prompt_for_group(atom_name, atom_df)

            progress_df = update_progress_row(
                progress_df,
                progress_csv,
                atom_name,
                reviewed=reviewed,
                group_label=group_label,
                notes=notes
            )

            save_grouping_summary_json(progress_df, progress_json)

    print('\nDone.')
    print(f'Progress CSV: {progress_csv}')
    print(f'Progress JSON: {progress_json}')
    print(f'Plots folder: {plot_dir}')


if __name__ == '__main__':
    main()