import os
import sys
import warnings
import tkinter as tk

from tkinter import filedialog

import numpy as np
import pandas as pd


warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module=r"PIL|matplotlib\.backends\._backend_tk",
)

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


EXCLUDE_KEYS = ['A', 'B', 'C']

TOP_N_EXTREMES = 100

SEVERE_PERCENT_THRESHOLD = 25.0
SEVERE_ABS_THRESHOLD = 1.0

SAVE_PREFIX = '3C_atom_volume_diagnostic'


def pick_folder(folder=None):
    if folder is not None:
        return folder

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    selected = filedialog.askdirectory()

    root.destroy()

    return selected


def normalize_residue_label(res):
    r = str(res).strip().upper()

    if r in {"A", "C", "G", "U", "T"}:
        if r == "A":
            return "DA"
        if r == "C":
            return "DC"
        if r == "G":
            return "DG"
        if r == "T":
            return "DT"
        if r == "U":
            return "DU"

    return r


def get_first_present(row, candidate_columns, default=''):
    for col in candidate_columns:
        if col in row.index:
            return row[col]

    return default


def normalize_chain_value(value):
    if pd.isna(value):
        return ''

    return str(value).strip()


def normalize_atom_name(value):
    if pd.isna(value):
        return ''

    return str(value).strip().upper()


def normalize_index(value):
    try:
        return int(value)
    except Exception:
        return value


def normalize_resseq(value):
    try:
        return int(value)
    except Exception:
        return value


def get_atom_identity(row):
    chain = normalize_chain_value(
        get_first_present(
            row,
            ['Chain', 'Chain ID', 'ChainID', 'Subunit'],
            default='',
        )
    )

    residue = normalize_residue_label(row['Residue'])
    resseq = normalize_resseq(row['Residue Sequence'])

    atom_name = normalize_atom_name(
        get_first_present(
            row,
            ['Atom', 'Atom Name', 'Name', 'Label'],
            default='',
        )
    )

    index = normalize_index(row['Index'])

    return (
        chain,
        residue,
        resseq,
        atom_name,
        index,
    )


def get_residue_identity(row):
    chain = normalize_chain_value(
        get_first_present(
            row,
            ['Chain', 'Chain ID', 'ChainID', 'Subunit'],
            default='',
        )
    )

    residue = normalize_residue_label(row['Residue'])
    resseq = normalize_resseq(row['Residue Sequence'])

    return (
        chain,
        residue,
        resseq,
    )


def check_duplicate_keys(atom_df, scheme_name, system_name):
    identity_counts = {}
    index_counts = {}

    for _, row in atom_df.iterrows():
        atom_key = get_atom_identity(row)
        idx = normalize_index(row['Index'])

        identity_counts[atom_key] = identity_counts.get(atom_key, 0) + 1
        index_counts[idx] = index_counts.get(idx, 0) + 1

    dup_identity = {k: v for k, v in identity_counts.items() if v > 1}
    dup_index = {k: v for k, v in index_counts.items() if v > 1}

    print(f'\n[{system_name} | {scheme_name}] duplicate check')
    print(f'  atoms total: {len(atom_df)}')
    print(f'  duplicated full atom keys: {len(dup_identity)}')
    print(f'  duplicated raw indices: {len(dup_index)}')

    if len(dup_identity) > 0:
        print('  first duplicated full keys:')
        for i, (k, v) in enumerate(dup_identity.items()):
            if i >= 10:
                break
            print(f'    {k} -> {v}')

    if len(dup_index) > 0:
        print('  first duplicated raw indices:')
        for i, (k, v) in enumerate(dup_index.items()):
            if i >= 10:
                break
            print(f'    {k} -> {v}')

    return dup_identity, dup_index


def build_lookup(atom_df):
    lookup = {}

    for _, row in atom_df.iterrows():
        key = get_atom_identity(row)

        if key in lookup:
            if isinstance(lookup[key], list):
                lookup[key].append(row)
            else:
                lookup[key] = [lookup[key], row]
        else:
            lookup[key] = row

    return lookup


def fetch_match(row, lookup):
    key = get_atom_identity(row)
    match = lookup.get(key, None)

    if isinstance(match, list):
        return match[0], len(match)

    if match is None:
        return None, 0

    return match, 1


def make_record(system_name, aw_row, other_row, scheme_name):
    aw_vol = float(aw_row['Volume'])
    other_vol = float(other_row['Volume'])

    delta = other_vol - aw_vol

    if abs(aw_vol) > 1e-12:
        signed_percent = (delta / aw_vol) * 100.0
    else:
        signed_percent = np.nan

    chain, residue, resseq = get_residue_identity(aw_row)
    atom_key = get_atom_identity(aw_row)

    return {
        'System': system_name,
        'Scheme': scheme_name,
        'Chain': chain,
        'Residue': residue,
        'Residue Sequence': resseq,
        'Atom Name': atom_key[3],
        'Index': atom_key[4],
        'AW Volume': aw_vol,
        f'{scheme_name} Volume': other_vol,
        f'{scheme_name} - AW': delta,
        f'{scheme_name} vs AW %': signed_percent,
        'Abs Delta': abs(delta),
        'Abs Percent': abs(signed_percent) if pd.notna(signed_percent) else np.nan,
    }


def summarize_residue_bias(df, scheme_name):
    if len(df) == 0:
        return pd.DataFrame()

    group_cols = ['Residue']

    summary = (
        df.groupby(group_cols)
        .agg(
            n=('Residue', 'size'),
            mean_percent=(f'{scheme_name} vs AW %', 'mean'),
            median_percent=(f'{scheme_name} vs AW %', 'median'),
            mean_abs_percent=('Abs Percent', 'mean'),
            mean_delta=(f'{scheme_name} - AW', 'mean'),
        )
        .reset_index()
        .sort_values('mean_percent', ascending=False)
    )

    return summary


def summarize_system_bias(df, scheme_name):
    if len(df) == 0:
        return pd.DataFrame()

    summary = (
        df.groupby(['System'])
        .agg(
            n=('System', 'size'),
            mean_percent=(f'{scheme_name} vs AW %', 'mean'),
            median_percent=(f'{scheme_name} vs AW %', 'median'),
            mean_abs_percent=('Abs Percent', 'mean'),
            mean_delta=(f'{scheme_name} - AW', 'mean'),
        )
        .reset_index()
        .sort_values('mean_percent', ascending=False)
    )

    return summary


def analyze_system(system_path, system_name):
    aw_path = os.path.join(system_path, 'aw', 'aw_logs.csv')
    pow_path = os.path.join(system_path, 'pow', 'pow_logs.csv')
    prm_path = os.path.join(system_path, 'prm', 'prm_logs.csv')

    if not (os.path.exists(aw_path) and os.path.exists(pow_path) and os.path.exists(prm_path)):
        print(f'Skipping {system_name}: missing one or more log files')
        return None

    aw_logs = read_logs2(aw_path, all_=False, balls=True, surfs=False)
    pow_logs = read_logs2(pow_path, all_=False, balls=True, surfs=False)
    prm_logs = read_logs2(prm_path, all_=False, balls=True, surfs=False)

    aw_atoms = aw_logs['atoms']
    pow_atoms = pow_logs['atoms']
    prm_atoms = prm_logs['atoms']

    check_duplicate_keys(aw_atoms, 'AW', system_name)
    check_duplicate_keys(pow_atoms, 'POW', system_name)
    check_duplicate_keys(prm_atoms, 'PRM', system_name)

    pow_lookup = build_lookup(pow_atoms)
    prm_lookup = build_lookup(prm_atoms)

    pow_records = []
    prm_records = []

    unmatched_pow = []
    unmatched_prm = []

    multi_pow = 0
    multi_prm = 0

    for _, aw_row in aw_atoms.iterrows():
        pow_match, pow_count = fetch_match(aw_row, pow_lookup)
        prm_match, prm_count = fetch_match(aw_row, prm_lookup)

        if pow_count > 1:
            multi_pow += 1

        if prm_count > 1:
            multi_prm += 1

        if pow_match is None:
            atom_key = get_atom_identity(aw_row)
            unmatched_pow.append({
                'System': system_name,
                'Chain': atom_key[0],
                'Residue': atom_key[1],
                'Residue Sequence': atom_key[2],
                'Atom Name': atom_key[3],
                'Index': atom_key[4],
                'AW Volume': float(aw_row['Volume']),
            })
        else:
            pow_records.append(make_record(system_name, aw_row, pow_match, 'POW'))

        if prm_match is None:
            atom_key = get_atom_identity(aw_row)
            unmatched_prm.append({
                'System': system_name,
                'Chain': atom_key[0],
                'Residue': atom_key[1],
                'Residue Sequence': atom_key[2],
                'Atom Name': atom_key[3],
                'Index': atom_key[4],
                'AW Volume': float(aw_row['Volume']),
            })
        else:
            prm_records.append(make_record(system_name, aw_row, prm_match, 'PRM'))

    print(f'\n[{system_name}] match summary')
    print(f'  AW atoms: {len(aw_atoms)}')
    print(f'  POW matched: {len(pow_records)} | unmatched: {len(unmatched_pow)} | multi-match keys used: {multi_pow}')
    print(f'  PRM matched: {len(prm_records)} | unmatched: {len(unmatched_prm)} | multi-match keys used: {multi_prm}')

    return {
        'pow_df': pd.DataFrame(pow_records),
        'prm_df': pd.DataFrame(prm_records),
        'unmatched_pow_df': pd.DataFrame(unmatched_pow),
        'unmatched_prm_df': pd.DataFrame(unmatched_prm),
    }


def print_extreme_report(df, scheme_name, top_n=20):
    if len(df) == 0:
        print(f'\nNo matched {scheme_name} records to report.')
        return

    print(f'\nTop {top_n} atoms where {scheme_name} is much larger than AW:')
    cols = [
        'System', 'Chain', 'Residue', 'Residue Sequence', 'Atom Name', 'Index',
        'AW Volume', f'{scheme_name} Volume', f'{scheme_name} - AW', f'{scheme_name} vs AW %'
    ]
    print(
        df.sort_values(f'{scheme_name} vs AW %', ascending=False)
        [cols]
        .head(top_n)
        .to_string(index=False)
    )

    print(f'\nTop {top_n} atoms where AW is larger than {scheme_name}:')
    print(
        df.sort_values(f'{scheme_name} vs AW %', ascending=True)
        [cols]
        .head(top_n)
        .to_string(index=False)
    )

    severe_small = df[
        (df[f'{scheme_name} vs AW %'] <= -SEVERE_PERCENT_THRESHOLD) |
        (df[f'{scheme_name} - AW'] <= -SEVERE_ABS_THRESHOLD)
    ].copy()

    severe_large = df[
        (df[f'{scheme_name} vs AW %'] >= SEVERE_PERCENT_THRESHOLD) |
        (df[f'{scheme_name} - AW'] >= SEVERE_ABS_THRESHOLD)
    ].copy()

    print(f'\nSevere cases where AW > {scheme_name}: {len(severe_small)}')
    if len(severe_small) > 0:
        print(severe_small[cols].head(top_n).to_string(index=False))

    print(f'\nSevere cases where {scheme_name} > AW: {len(severe_large)}')
    if len(severe_large) > 0:
        print(severe_large[cols].head(top_n).to_string(index=False))


def main():
    folder = pick_folder()

    if folder is None or folder == '':
        print('No folder selected.')
        return

    all_pow = []
    all_prm = []
    all_unmatched_pow = []
    all_unmatched_prm = []

    for subfolder in os.listdir(folder):
        system_path = os.path.join(folder, subfolder)

        if not os.path.isdir(system_path):
            continue

        sys_key = subfolder.split('_')[0]

        if sys_key in EXCLUDE_KEYS:
            continue

        print('\n' + '=' * 100)
        print(f'Analyzing {subfolder}')
        print('=' * 100)

        result = analyze_system(system_path, subfolder)

        if result is None:
            continue

        if len(result['pow_df']) > 0:
            all_pow.append(result['pow_df'])

        if len(result['prm_df']) > 0:
            all_prm.append(result['prm_df'])

        if len(result['unmatched_pow_df']) > 0:
            all_unmatched_pow.append(result['unmatched_pow_df'])

        if len(result['unmatched_prm_df']) > 0:
            all_unmatched_prm.append(result['unmatched_prm_df'])

    pow_df = pd.concat(all_pow, ignore_index=True) if len(all_pow) > 0 else pd.DataFrame()
    prm_df = pd.concat(all_prm, ignore_index=True) if len(all_prm) > 0 else pd.DataFrame()
    unmatched_pow_df = pd.concat(all_unmatched_pow, ignore_index=True) if len(all_unmatched_pow) > 0 else pd.DataFrame()
    unmatched_prm_df = pd.concat(all_unmatched_prm, ignore_index=True) if len(all_unmatched_prm) > 0 else pd.DataFrame()

    print('\n' + '#' * 100)
    print('GLOBAL SUMMARY')
    print('#' * 100)

    print(f'\nMatched POW atom comparisons: {len(pow_df)}')
    print(f'Matched PRM atom comparisons: {len(prm_df)}')
    print(f'Unmatched AW->POW atoms: {len(unmatched_pow_df)}')
    print(f'Unmatched AW->PRM atoms: {len(unmatched_prm_df)}')

    if len(pow_df) > 0:
        print(f'\nPOW signed percent summary:')
        print(pow_df['POW vs AW %'].describe())

    if len(prm_df) > 0:
        print(f'\nPRM signed percent summary:')
        print(prm_df['PRM vs AW %'].describe())

    print_extreme_report(pow_df, 'POW', top_n=TOP_N_EXTREMES)
    print_extreme_report(prm_df, 'PRM', top_n=TOP_N_EXTREMES)

    pow_residue_summary = summarize_residue_bias(pow_df, 'POW')
    prm_residue_summary = summarize_residue_bias(prm_df, 'PRM')

    pow_system_summary = summarize_system_bias(pow_df, 'POW')
    prm_system_summary = summarize_system_bias(prm_df, 'PRM')

    if len(pow_residue_summary) > 0:
        print('\nResidues with most positive mean POW vs AW %:')
        print(pow_residue_summary.head(20).to_string(index=False))

        print('\nResidues with most negative mean POW vs AW %:')
        print(pow_residue_summary.sort_values('mean_percent', ascending=True).head(20).to_string(index=False))

    if len(prm_residue_summary) > 0:
        print('\nResidues with most positive mean PRM vs AW %:')
        print(prm_residue_summary.head(20).to_string(index=False))

        print('\nResidues with most negative mean PRM vs AW %:')
        print(prm_residue_summary.sort_values('mean_percent', ascending=True).head(20).to_string(index=False))

    if len(pow_system_summary) > 0:
        print('\nSystems sorted by mean POW vs AW %:')
        print(pow_system_summary.to_string(index=False))

    if len(prm_system_summary) > 0:
        print('\nSystems sorted by mean PRM vs AW %:')
        print(prm_system_summary.to_string(index=False))

    pow_csv = os.path.join(folder, f'{SAVE_PREFIX}_pow_atom_comparisons.csv')
    prm_csv = os.path.join(folder, f'{SAVE_PREFIX}_prm_atom_comparisons.csv')
    unmatched_pow_csv = os.path.join(folder, f'{SAVE_PREFIX}_unmatched_aw_to_pow.csv')
    unmatched_prm_csv = os.path.join(folder, f'{SAVE_PREFIX}_unmatched_aw_to_prm.csv')
    pow_residue_csv = os.path.join(folder, f'{SAVE_PREFIX}_pow_residue_summary.csv')
    prm_residue_csv = os.path.join(folder, f'{SAVE_PREFIX}_prm_residue_summary.csv')
    pow_system_csv = os.path.join(folder, f'{SAVE_PREFIX}_pow_system_summary.csv')
    prm_system_csv = os.path.join(folder, f'{SAVE_PREFIX}_prm_system_summary.csv')

    if len(pow_df) > 0:
        pow_df.to_csv(pow_csv, index=False)

    if len(prm_df) > 0:
        prm_df.to_csv(prm_csv, index=False)

    if len(unmatched_pow_df) > 0:
        unmatched_pow_df.to_csv(unmatched_pow_csv, index=False)

    if len(unmatched_prm_df) > 0:
        unmatched_prm_df.to_csv(unmatched_prm_csv, index=False)

    if len(pow_residue_summary) > 0:
        pow_residue_summary.to_csv(pow_residue_csv, index=False)

    if len(prm_residue_summary) > 0:
        prm_residue_summary.to_csv(prm_residue_csv, index=False)

    if len(pow_system_summary) > 0:
        pow_system_summary.to_csv(pow_system_csv, index=False)

    if len(prm_system_summary) > 0:
        prm_system_summary.to_csv(prm_system_csv, index=False)

    print('\nSaved diagnostic files to selected folder.')


if __name__ == '__main__':
    main()