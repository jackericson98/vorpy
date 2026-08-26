import os
import sys
import numpy as np
import tkinter as tk
from tkinter import filedialog

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', '..', '..', '..', '..'))
sys.path.append(vorpy_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def get_res_sa(logs):
    """Calculate external surface area for each residue."""
    res_data, res_key_dict = {}, {}

    for _, atom in logs['atoms'].iterrows():
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        res_data.setdefault(res_key, []).append(atom['Index'])
        res_key_dict[atom['Index']] = res_key

    res_sa = {key: 0.0 for key in res_data}

    for _, surf in logs['surfs'].iterrows():
        ball_1, ball_2 = surf['Balls']
        res_1 = res_key_dict.get(ball_1)
        res_2 = res_key_dict.get(ball_2)

        # Internal to one residue -> don't count
        if res_1 is not None and res_1 == res_2:
            continue

        # Surface contributes to each residue it borders
        if res_1 is not None:
            res_sa[res_1] += surf['Surface Area']
        if res_2 is not None:
            res_sa[res_2] += surf['Surface Area']

    return res_sa


def get_res_vol(logs):
    """Calculate total volume for each residue."""
    res_vol = {}

    for _, atom in logs['atoms'].iterrows():
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        res_vol.setdefault(res_key, 0.0)
        res_vol[res_key] += atom['Volume']

    return res_vol


def main():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    file = filedialog.askopenfilename(title='Select logs file', filetypes=[('CSV files', '*.csv'), ('All files', '*.*')])
    if not file:
        return

    logs = read_logs2(file, all_=False, balls=True, surfs=True)

    res_vol = get_res_vol(logs)
    res_sa = get_res_sa(logs)

    volumes = np.array(list(res_vol.values()))
    surface_areas = np.array(list(res_sa.values()))

    print(f"\nResidues: {len(volumes)}")
    print(f"Mean Volume:       {np.mean(volumes):.3f} ± {np.std(volumes, ddof=1):.3f} Å³")
    print(f"Mean Surface Area: {np.mean(surface_areas):.3f} ± {np.std(surface_areas, ddof=1):.3f} Å²")


if __name__ == '__main__':
    main()