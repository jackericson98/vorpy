import os
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="PIL")


# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.line import line_plot


def get_frames_by_group(folder=None, chains=None):
    # Get the folder if it is not provided
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory()
    # Create the frame dictionary
    frame_dict = {}
    max_diffs = {}
    # Loop through the folders in the folder
    for subfolder in os.listdir(folder):
        if subfolder[0] != 'f':
            continue
        # Get the name of the frame
        frame = int(subfolder[1:])
        # Get the logs
        aw = os.path.join(folder, subfolder, 'aw/aw_logs.csv')
        pw = os.path.join(folder, subfolder, 'pow/pow_logs.csv')
        pm = os.path.join(folder, subfolder, 'prm/prm_logs.csv')
        # read the logs
        aw_logs = read_logs2(aw, all_=False, balls=True)
        pw_logs = read_logs2(pw, all_=False, balls=True)
        pm_logs = read_logs2(pm, all_=False, balls=True)
        # Add the volumes to the frame dictionary
        frame_dict[frame] = {'aw': 0, 'pow': 0, 'prm': 0}
        max_diffs[frame] = {'pow': (-1, 0), 'prm': (-1, 0)}
        # Get the atoms
        for i, aw_atom in aw_logs['atoms'].iterrows():
            if chains is not None and aw_atom['Chain'] not in chains:
                continue
            # Get the pow information
            pw_atom = pw_logs['atoms'].loc[pw_logs['atoms']['Index'] == aw_atom['Index']].to_dict()
            pm_atom = pm_logs['atoms'].loc[pm_logs['atoms']['Index'] == aw_atom['Index']].to_dict()
            # Get the volume
            aw_vol = aw_atom['Volume']
            pw_vol = pw_atom['Volume'][aw_atom['Index']]
            pm_vol = pm_atom['Volume'][aw_atom['Index']]

            max_diffs[frame]['pow'] = (aw_atom['Index'], max(max_diffs[frame]['pow'][1], abs(aw_vol - pw_vol) / aw_vol))
            max_diffs[frame]['prm'] = (aw_atom['Index'], max(max_diffs[frame]['prm'][1], abs(aw_vol - pm_vol) / aw_vol))

            # Add the frames to the frame dictionary
            frame_dict[frame]['aw'] += aw_vol
            frame_dict[frame]['pow'] += pw_vol
            frame_dict[frame]['prm'] += pm_vol
    print(max_diffs)
    return frame_dict


def plot_frame_diffs(frame_dict):
    # Create the lists
    pow_diffs, prm_diffs, frame_names = [], [], []
    # Loop through the frames
    for frame in frame_dict:
        # Get the volumes
        aw_vol = frame_dict[frame]['aw']
        pw_vol = frame_dict[frame]['pow']
        pm_vol = frame_dict[frame]['prm']
        # Calculate the differences
        pow_diff = 100 * (aw_vol - pw_vol) / aw_vol
        prm_diff = 100 * (aw_vol - pm_vol) / aw_vol
        # Add the differences to the lists
        pow_diffs.append(pow_diff)
        prm_diffs.append(prm_diff)
        frame_names.append(frame)

    # INSERT_YOUR_CODE
    # Sort the diffs and frame_names by frame number (ascending)
    sorted_data = sorted(zip(frame_names, pow_diffs, prm_diffs))
    frame_names, pow_diffs, prm_diffs = zip(*sorted_data)
    frame_names = list(frame_names)
    pow_diffs = list(pow_diffs)
    prm_diffs = list(prm_diffs)
    # Plot the differences
    line_plot([frame_names, frame_names], [pow_diffs, prm_diffs], y_label='% Difference', 
              title='% Difference from AW', Show=True, x_label="Frame", labels=["Pow", "Prm"], 
              colors=["red", "purple"], y_ticks=[-2, -1, 0, 1, 2, 3], y_label_size=20, x_label_size=20, 
              title_size=25, xlim=[0.5, 11.5], x_ticks=[2, 4, 6, 8, 10], axis_line_thickness=2,
              linewidth=3, tick_val_size=20)


def plot_frame_deviations(frame_dict):
    # Create the lists
    aw_devs, pow_devs, prm_devs, frame_names = [], [], [], []
    # Get the first frame
    aw_vol = frame_dict[1]['aw']
    pw_vol = frame_dict[1]['pow']
    pm_vol = frame_dict[1]['prm']
    # Loop through the frames
    for frame in frame_dict:
        # Get the volumes
        aw_vol_ = frame_dict[frame]['aw']
        pw_vol_ = frame_dict[frame]['pow']
        pm_vol_ = frame_dict[frame]['prm']
        # Calculate the deviations
        aw_dev = 100 * (aw_vol_ - aw_vol) / aw_vol
        pow_dev = 100 * (pw_vol_ - pw_vol) / pw_vol
        prm_dev = 100 * (pm_vol_ - pm_vol) / pm_vol
        # Add the deviations to the lists
        aw_devs.append(aw_dev)
        pow_devs.append(pow_dev)
        prm_devs.append(prm_dev)
        frame_names.append(frame)
    # Sort the diffs and frame_names by frame number (ascending)
    sorted_data = sorted(zip(frame_names, pow_devs, prm_devs))
    frame_names, pow_devs, prm_devs = zip(*sorted_data)
    frame_names = list(frame_names)
    pow_devs = list(pow_devs)
    prm_devs = list(prm_devs)
    # Plot the deviations
    line_plot([frame_names, frame_names, frame_names], [aw_devs, pow_devs, prm_devs], y_label='% Deviation', 
              title='Frame Volume Deviation', Show=True, x_label="Frame", labels=["AW", "Pow", "Prm"], 
              colors=["blue", "red", "purple"], y_ticks=[-3, -2, -1, 0], y_label_size=20, x_label_size=20, 
              title_size=25, xlim=[0.5, 11.5], x_ticks=[2, 4, 6, 8, 10], axis_line_thickness=2, 
              linewidth=3, tick_val_size=20, alpha=0.8)


if __name__ == "__main__":
    frame_dict = get_frames_by_group(chains=['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'])
    plot_frame_diffs(frame_dict)
    plot_frame_deviations(frame_dict)
