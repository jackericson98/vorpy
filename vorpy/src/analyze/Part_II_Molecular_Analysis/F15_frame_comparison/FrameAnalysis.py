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


def get_frames_by_group(folder=None):
    # Get the folder if it is not provided
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory()
    # Create the frame dictionary
    frame_dict = {}

    # Loop through the folders in the folder
    for subfolder in os.listdir(folder):
        if subfolder[0] != 'f':
            continue
        # Get the name of the frame
        frame = int(subfolder[1:])
        # read the logs
        aw_logs = read_logs2(os.path.join(folder, subfolder, 'aw/aw_logs.csv'), all_=False, balls=False)
        pw_logs = read_logs2(os.path.join(folder, subfolder, 'pow/pow_logs.csv'), all_=False, balls=False)
        pm_logs = read_logs2(os.path.join(folder, subfolder, 'prm/prm_logs.csv'), all_=False, balls=False)
        # Get the volume and SA
        aw_vol, aw_sa = aw_logs['group data']['Volume'], aw_logs['group data']['Surface Area']
        pw_vol, pw_sa = pw_logs['group data']['Volume'], pw_logs['group data']['Surface Area']
        pm_vol, pm_sa = pm_logs['group data']['Volume'], pm_logs['group data']['Surface Area']

        # Add the volumes to the frame dictionary
        frame_dict[frame] = {'vol': {'aw': aw_vol, 'pow': pw_vol, 'prm': pm_vol},
                             'sa': {'aw': aw_sa, 'pow': pw_sa, 'prm': pm_sa}}

    return frame_dict


def plot_frame_diffs(frame_dict):
    # Create the lists
    pow_diffs, prm_diffs, frame_names = [], [], []
    # Loop through the frames
    for frame in frame_dict:
        # Get the volumes
        aw_vol = frame_dict[frame]['vol']['aw']
        pw_vol = frame_dict[frame]['vol']['pow']
        pm_vol = frame_dict[frame]['vol']['prm']
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
    aw_vol = frame_dict[1]['vol']['aw']
    pw_vol = frame_dict[1]['vol']['pow']
    pm_vol = frame_dict[1]['vol']['prm']
    # Loop through the frames
    for frame in frame_dict:
        # Get the volumes
        aw_vol_ = frame_dict[frame]['vol']['aw']
        pw_vol_ = frame_dict[frame]['vol']['pow']
        pm_vol_ = frame_dict[frame]['vol']['prm']
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
    frame_dict = get_frames_by_group()
    plot_frame_diffs(frame_dict)
    plot_frame_deviations(frame_dict)
