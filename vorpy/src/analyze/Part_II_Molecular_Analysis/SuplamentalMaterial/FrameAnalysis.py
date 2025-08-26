import os
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.tools.plot_templates.line import line_plot


def get_frames(folder=None):
    # Get the folder if it is not provided
    if folder is None:
        root = tk.Tk()
        root.withdraw()
        folder = filedialog.askdirectory()
    # Create the frame dictionary
    frame_dict = {}
    # Loop through the folders in the folder
    for subfolder in os.listdir(folder):
        # Get the name of the frame
        frame = int(subfolder[1:])
        # Get the logs
        aw = os.path.join(folder, subfolder, 'aw/aw_logs.csv')
        pw = os.path.join(folder, subfolder, 'pow/pow_logs.csv')
        pm = os.path.join(folder, subfolder, 'prm/prm_logs.csv')
        # read the logs
        aw_logs = read_logs2(aw, all_=False)
        pw_logs = read_logs2(pw, all_=False)
        pm_logs = read_logs2(pm, all_=False)
        # Get the volumes
        aw_vol = aw_logs['group data']['Volume']
        pw_vol = pw_logs['group data']['Volume']
        pm_vol = pm_logs['group data']['Volume']
        # Add the frames to the frame dictionary
        frame_dict[frame] = {'aw': aw_vol, 'pow': pw_vol, 'prm': pm_vol}
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
    
    # Plot the differences
    line_plot([frame_names, frame_names], [pow_diffs, prm_diffs], y_label='% Difference', 
              title='% Difference from AW', Show=True, x_label="Frame", labels=["Pow", "Prm"], 
              colors=["red", "purple"], y_ticks=[-3, -2, -1, 0, 1], y_label_size=20, x_label_size=20, 
              title_size=25, ylim=[-3.5, 1.5], xlim=[0.5, 11.5], x_ticks=[2, 4, 6, 8, 10], axis_line_thickness=2,
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
    # Plot the deviations
    line_plot([frame_names, frame_names, frame_names], [aw_devs, pow_devs, prm_devs], y_label='% Deviation', 
              title='Frame Volume Deviation', Show=True, x_label="Frame", labels=["AW", "Pow", "Prm"], 
              colors=["blue", "red", "purple"], y_ticks=[-3, -2, -1, 0], y_label_size=20, x_label_size=20, 
              title_size=25, ylim=[-4.5, 0.5], xlim=[0.5, 11.5], x_ticks=[2, 4, 6, 8, 10], axis_line_thickness=2, 
              linewidth=3, tick_val_size=20, alpha=0.8)

if __name__ == "__main__":
    frame_dict = get_frames()
    plot_frame_diffs(frame_dict)
    plot_frame_deviations(frame_dict)

