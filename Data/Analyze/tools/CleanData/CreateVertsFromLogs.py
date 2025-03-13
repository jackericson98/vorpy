import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import shutil
import tkinter as tk
from tkinter import filedialog
from System.sys_funcs.output.net import write_verts
from Data.Analyze.tools.compare.read_logs2 import read_logs2


def write_log_to_vert(in_file, out_file):
    # Get the logs
    logs = read_logs2(in_file, verts=True)
    # Write the aw_verts
    with open(out_file, 'w') as file:
        # Create a header for the vertices file
        file.write("Vertices - {} vertices, {} atoms, max vert = {}, Net type = {}\n"
                   .format(len(logs['verts']['Index']), 1000, max(logs['verts']['rad']),
                           'aw'))
        # Write the vertices
        for i, vert in logs['verts'].iterrows():
            # Write the vertex
            file.write(" ".join([str(_) for _ in vert['Balls']]) + " " + " ".join([str(_) for _ in vert['loc']]) +
                       " " + str(vert['rad']) + " \n")
        # Write the end line for the file
        file.write("END")


def create_verts_from_logs(folder=None):

    if folder is None:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()

    for subfolder in os.listdir(folder):
        # Create the aw folder
        if not os.path.exists(folder + '/' + subfolder + '/aw'):
            os.mkdir(folder + '/' + subfolder + '/aw')
        # Create the pow folder
        if not os.path.exists(folder + '/' + subfolder + '/pow'):
            os.mkdir(folder + '/' + subfolder + '/pow')

        # For each subfolder we need to first get the logs and then read them
        # aw_path1 = folder + '/' + subfolder + '/aw/aw_logs.csv'
        # aw_path2 = folder + '/' + subfolder + '/aw_logs.csv'
        # aw_path1 = folder + '/' + subfolder + '/aw/aw_logs.csv'
        # aw_path2 = folder + '/' + subfolder + '/aw_logs.csv'
        # if os.path.exists(aw_path1):
        try:
            shutil.move(folder + '/' + subfolder + '/aw_logs.csv', folder + '/' + subfolder + '/aw/aw_logs.csv')
        except FileNotFoundError:
            pass
        try:
            shutil.move(folder + '/' + subfolder + '/pow_logs.csv', folder + '/' + subfolder + '/pow/pow_logs.csv')
        except FileNotFoundError:
            pass
        # Get the verts
        write_log_to_vert(folder + '/' + subfolder + '/aw/aw_logs.csv', folder + '/' + subfolder + '/aw/aw_verts.txt')
        write_log_to_vert(folder + '/' + subfolder + '/pow/pow_logs.csv', folder + '/' + subfolder + '/pow/pow_verts.txt')


if __name__ == '__main__':
    create_verts_from_logs()
