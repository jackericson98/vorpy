import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from System.sys_funcs.calcs.calcs import calc_dist
import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.read_logs2 import read_logs2


def get_information(vpy_fl=None, vta_fl=None, vvv_fl=None, pdb_fl=None):
    # Get the files
    if vpy_fl is None:
        vpy_fl = filedialog.askopenfilename(title="Get Vorpy Logs")
    if vta_fl is None:
        vta_fl = filedialog.askopenfilename(title="Get Voronota Volumes")
    if pdb_fl is None:
        pdb_fl = filedialog.askopenfilename(title="Get pdb file")

    # Go through the vorpy vertices
    with open(vpy_fl, 'r') as vpy_file:
        # Read the logs
        vpy_logs = read_logs2(vpy_file, all_=False, balls=True)
        # Set up the vorpy volumes dictionary
        vpy_vols = {}
        # Go through the logs and make a matching style to voronota
        for i, ball in vpy_logs['atoms'].iterrows():
            # Get the vorpy volume
            vpy_vols[ball['Index']] = ball['Volume']

    # Go through the voronota vertices
    with open(vta_fl, 'r') as vta_file:
        # Create the dictionary
        vta_vols = {}
        # Create a line counter
        line_counter = -1
        # Loop through the lines
        for vta_line in vta_file.readlines():
            # Increment the line counter
            line_counter += 1
            # Split the line
            line_list = [_ for _ in vta_line.split(" ") if _ != ""]
            try:
                # Get the ball index
                ndx = int(line_list[0])
            except ValueError:
                continue
            if ndx >= 1000:
                break
            vta_vols[ndx] = float(line_list[1])

    # Go through the pdb file
    with open(pdb_fl, 'r') as pdb_file:
        # Start the counter
        counter = 0
        # Open the pdb file and get the first line
        for line in pdb_file.readlines():
            # Only get the info from the first line
            if counter > 0:
                break
            # Get the first line information
            line_info = line.split(',')
            # Add the information to the pdb info by hand:
            pdb_info = {'box': float(line_info[0].split(' ')[-1]), 'avg': float(line_info[1].split(' ')[-1]),
                        'cv_': float(line_info[2].split(' ')[-1]), 'num': float(line_info[3].split(' ')[-1]),
                        'den': float(line_info[4].split(' ')[-1]), 'olp': float(line_info[5].split(' ')[-1][:-1]),
                        'dst': line_info[6].split(' ')[-1], 'pbc': line_info[7].split(' ')[-1] == 'True',
                        'sar': line_info[8].split(' ')[-1] == 'True', 'file': pdb_fl}
            # Increment the counter
            counter += 1

    # return the three dictionaries
    return vpy_vols, vta_vols, pdb_info


def compare_vertices(dic1, dic2, pdb, dic1_name='Vorpy', dic2_name=''):
    # Get the overlap indices
    olap_ndxs = set([_ for _ in dic1.keys() if _ in dic2.keys() and _ != 'info'])
    # Create the comparison dictionary
    comparison = {'pdb': pdb,
                  'info': {'dic1 extra verts': len(dic1) - len(olap_ndxs),
                           'dic2 extra verts': len(dic2) - len(olap_ndxs),
                           'dic1 name': dic1_name, 'dic2 name': dic2_name},
                  'diffs': {}}
    # Get the diffs
    for _ in olap_ndxs:
        # Get the locations and radii
        loc1, rad1, dub1 = dic1[_]['loc'], dic1[_]['rad'], dic1[_]['dub']
        loc2, rad2, dub2 = dic2[_]['loc'], dic2[_]['rad'], dic2[_]['dub']
        # Get the location distance and radii difference
        comparison['diffs'][_] = {'loc': calc_dist(loc1, loc2), 'rad': rad1 - rad2, 'dub': dub1 == dub2}
    # Return the dictionary
    return comparison


def plot_comparison(comparison_dictionary, folder=None):
    """
    Plotting the difference in radius and location by vertex and place a venn diagram of the
    differences.
    """

    loc_dist_y, rad_diff_y = [], []
    for ndxs in comparison_dictionary['diffs']:
        if comparison_dictionary['diffs'][ndxs]['loc'] > 1 or comparison_dictionary['diffs'][ndxs]['rad'] > 1:
            continue
        loc_dist_y.append(comparison_dictionary['diffs'][ndxs]['loc'])
        rad_diff_y.append(comparison_dictionary['diffs'][ndxs]['rad'])

    # Create the x values
    xs = np.arange(0, len(loc_dist_y))

    # Create figure and scatter plot
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.scatter(xs, rad_diff_y, label="Tangential Sphere\nRadius Difference", color='blue', alpha=0.3, s=2)
    ax.scatter(xs, loc_dist_y, label="Euclidean Distance", color='red', alpha=0.3, s=2)

    # Setting labels with specified font sizes
    ax.set_xlabel("Vertex", fontsize=20)
    ax.set_ylabel("Difference", fontsize=20)
    ax.set_ylim([-0.1, 0.2])
    ax.set_xticks([])

    # Setting tick parameters directly
    ax.tick_params(axis='y', labelsize=20)

    ax.set_title(f"{comparison_dictionary['info']['dic1 name']} vs {comparison_dictionary['info']['dic2 name']} Vertex Comparison", fontsize=25)
    ax.legend(fontsize=15, loc='upper right', bbox_to_anchor=(1, 0.9), markerscale=5)

    # Create Venn diagram in top right
    venn_ax = fig.add_axes([0.2, 0.55, 0.3, 0.3])  # [left, bottom, width, height]
    venn = venn2(subsets=(np.log(comparison_dictionary['info']['dic1 extra verts']),
                          np.log(comparison_dictionary['info']['dic2 extra verts']),
                          np.log(len(loc_dist_y))),
                 set_labels=(comparison_dictionary['info']['dic1 name'],
                             comparison_dictionary['info']['dic2 name']),
                 ax=venn_ax)
    venn_ax.set_title('Overlapping/Missing Vertices', fontsize=14)

    # Add numbers to the Venn diagram sections
    venn.get_label_by_id('10').set_text(str(comparison_dictionary['info']['dic1 extra verts']))
    venn.get_label_by_id('01').set_text(str(comparison_dictionary['info']['dic2 extra verts']))
    venn.get_label_by_id('11').set_text(str(len(loc_dist_y)))

    # Additional annotations
    plt.text(0.05, y=-2.5, s=f"Retaining Cube Side Length = {comparison_dictionary['pdb']['box']}", fontsize=15)

    # Check if we want to save or show the plot
    if folder is None:
        plt.show()
    else:
        subfolder = os.path.dirname(comparison_dictionary['pdb']['file'])
        print(subfolder)
        plt.savefig(os.path.join(folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.svg"), format='svg')

        plt.savefig(os.path.join(folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.png"), format='png', dpi=600)
        plt.close(fig)


def run_it_all():
    # Get the folder that holds the subfolders
    folder = filedialog.askdirectory()
    # Loop through the directory
    for subfolder in os.listdir(folder):
        # Get the aw verts
        vpy_verts = os.path.join(folder, subfolder, 'Vorpy', 'aw_verts.txt')
        if not os.path.exists(vpy_verts):
            print(os.path.join(folder, subfolder, 'Vorpy', 'aw_verts.txt'))
            vpy_verts = None
        # Get the Voronota verts
        vta_verts = os.path.join(folder, subfolder, 'Voronota', 'vertices.txt')
        if not os.path.exists(vta_verts):
            print(os.path.join(folder, subfolder, 'Voronota', 'vertices.txt'))
            vta_verts = None
        # Get the V verts
        vvv_verts = os.path.join(folder, subfolder, 'V', 'V_Vertices.txt')
        if not os.path.exists(vvv_verts):
            print(os.path.join(folder, subfolder, 'V', 'V_Vertices.txt'))
            vvv_verts = None
        # Get the pdb file
        pdb_file = [os.path.join(folder, subfolder, file) for file in os.listdir(os.path.join(folder, subfolder)) if file[-4:] == '.pdb'][0]
        if vpy_verts is None or vta_verts is None or vvv_verts is None:
            continue
        # Compare the files now
        vpy, vta, pdb = get_information(vpy_fl=vpy_verts, vta_fl=vta_verts, vvv_fl=vvv_verts, pdb_fl=pdb_file)

        # Compare the files
        comparison_dictionary = compare_vertices(vpy, vta, pdb, 'Vorpy', 'Voronota')
        # Plot the comparison
        plot_comparison(comparison_dictionary, folder=os.path.join(folder, subfolder))




if __name__ == '__main__':

    # First open the folder with all the data
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    # vorpy, voronota, V, pdb_dict = get_information()
    # comp_dic = compare_vertices(vorpy, voronota, pdb_dict, 'Vorpy', 'Voronota')
    # plot_comparison(comp_dic)
    run_it_all()

