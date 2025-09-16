import os
import sys

import tkinter as tk
from tkinter import filedialog

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../.."))
sys.path.append(project_root)

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def verify_frames(frame_folder):
    # Get the files from the frame folder
    folders = os.listdir(frame_folder)
    frame_verification = {}
    # Get the files from the frame folder
    for folder in folders:
        if folder[-3:] == 'pdb':
            continue
        frame_verification[folder] = {}
        # Check for the aw_logs, pow_logs, and prm_logs
        if os.path.exists(frame_folder + '/' + folder + '/aw/aw_logs.csv'):
            aw_logs = read_logs2(frame_folder + '/' + folder + '/aw/aw_logs.csv', all_=False, balls=True)
            frame_verification[folder]['aw_logs'] = len(aw_logs['atoms'])
        else:
            frame_verification[folder]['aw_logs'] = False
            continue

        if os.path.exists(frame_folder + '/' + folder + '/pow/pow_logs.csv'):
            pow_logs = read_logs2(frame_folder + '/' + folder + '/pow/pow_logs.csv', all_=False, balls=True)
            frame_verification[folder]['pow_logs'] = len(pow_logs['atoms'])
        else:
            frame_verification[folder]['pow_logs'] = False
            continue
        if os.path.exists(frame_folder + '/' + folder + '/prm/prm_logs.csv'):
            prm_logs = read_logs2(frame_folder + '/' + folder + '/prm/prm_logs.csv', all_=False, balls=True)
            frame_verification[folder]['prm_logs'] = len(prm_logs['atoms'])
        else:
            frame_verification[folder]['prm_logs'] = False
            continue
    # Print the breakdown
    for frame in frame_verification:
        print(f"Frame {frame}: AW - {frame_verification[frame]['aw_logs']}, POW - {frame_verification[frame]['pow_logs']}, PRM - {frame_verification[frame]['prm_logs']}")

    return frame_verification


if __name__ == "__main__":
    verify_frames(frame_folder=filedialog.askdirectory())
