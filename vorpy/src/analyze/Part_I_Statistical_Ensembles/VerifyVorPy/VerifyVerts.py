import os
import csv
import tkinter as tk
from tkinter import filedialog
import pandas as pd



import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt


def calc_dist(l0, l1):
    """Calculate the Euclidean distance between two points in n-dimensional space.

    Parameters
    ----------
    l0 : array-like
        First point coordinates as an n-dimensional array or list
    l1 : array-like
        Second point coordinates as an n-dimensional array or list with same dimensionality as l0

    Returns
    -------
    float
        The Euclidean distance between the two points

    Examples
    --------
    >>> calc_dist([0, 0, 0], [1, 1, 1])
    1.7320508075688772
    >>> import numpy as np
    >>> calc_dist(np.array([0, 0, 0]), np.array([1, 1, 1]))
    1.7320508075688772
    """

    return np.sqrt(sum(np.square(np.array(l0) - np.array(l1))))


def get_information(vpy_fl=None, vta_fl=None, pdb_fl=None):
    # Get the files
    if vpy_fl is None:
        vpy_fl = filedialog.askopenfilename(title="Get Vorpy Vertices")
    if vta_fl is None:
        vta_fl = filedialog.askopenfilename(title="Get Voronota Vertices")
    if pdb_fl is None:
        pdb_fl = filedialog.askopenfilename(title="Get pdb file")

    # Go through the vorpy vertices
    with open(vpy_fl, 'r') as vpy_file:
        # Create the dictionary
        vpy_vrts = {}
        # Create a line counter
        line_counter = -1
        # Loop through the lines
        for vpy_line in vpy_file.readlines():
            # Increment the line counter
            line_counter += 1
            # Split the line
            line_list = [_ for _ in vpy_line.split(" ") if _ != ""]

            # Read the first line
            if line_counter == 0:
                # Get some vorpy information
                vpy_vrts['info'] = {
                    'Number': int(line_list[2]),
                    'Max Vert': float(line_list[9][:-1])
                }
                # Continue onto the next line
                continue
            if line_list[0] == 'END':
                continue
            # Get the vertex indices
            ndxs = tuple([int(_) for _ in line_list[:4]])
            # Get the vertex location
            loc, rad, dub = tuple([float(_) for _ in line_list[4:7]]), float(line_list[7]), int(line_list[8][0])
            # Check if the vertex is a doublet
            if dub == 1:
                vpy_vrts[ndxs]['loc2'], vpy_vrts['rad2'], vpy_vrts['dub'] = loc, rad, True
                continue
            # Add the vertex to the dictionary
            vpy_vrts[ndxs] = {'loc': loc, 'rad': rad, 'loc2': None, 'rad2': None, 'dub': False}

    # Go through the voronota vertices
    with open(vta_fl, 'r') as vta_file:
        # Create the dictionary
        vta_vrts = {}
        # Create a line counter
        line_counter = -1
        # Create the vta_verts info dictionary
        vta_vrts['info'] = {"Number": 0, "Max Vert": 0, 'Larger Verts Count': 0, 'Larger Verts': {}}
        # Loop through the lines
        for vta_line in vta_file.readlines():
            # Increment the line counter
            line_counter += 1
            # Split the line
            line_list = [_ for _ in vta_line.split(" ") if _ != ""]
            # Get the vertex indices
            try:
                ndxs = tuple([int(_) for _ in line_list[:4]])
            except ValueError:
                continue
            # Check that the vertices have at least one value below 1000
            if not any([_ < 1000 for _ in ndxs]):
                continue
            # Get the vertex location
            loc, rad = tuple([float(_) for _ in line_list[4:7]]), float(line_list[7])
            # Check if the radius is less than the maximum
            if rad > vpy_vrts['info']['Max Vert']:
                if rad > vta_vrts['info']['Max Vert']:
                    vta_vrts['info']['Max Vert'] = rad
                vta_vrts['info']['Larger Verts Count'] += 1
                vta_vrts['info']['Larger Verts'][ndxs] = {'loc': loc, 'rad': rad, 'loc2': None, 'rad2': None, 'dub': False}
            # Check if the vertex is a doublet
            if ndxs in vta_vrts:
                if rad < vta_vrts[ndxs]['rad']:
                    vta_vrts[ndxs] = {'loc': loc, 'rad': rad, 'loc2': vta_vrts[ndxs]['loc'], 'rad2': vta_vrts[ndxs]['rad'], 'dub': True}
                else:
                    vta_vrts[ndxs]['loc2'], vta_vrts[ndxs]['rad2'], vta_vrts[ndxs]['dub'] = loc, rad, True
                continue
            # Add the vertex to the dictionary
            vta_vrts[ndxs] = {'loc': loc, 'rad': rad, 'loc2': None, 'rad2': None, 'dub': False}

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
    return vpy_vrts, vta_vrts, pdb_info


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
        comparison['diffs'][_] = {'loc': calc_dist(loc1, loc2), 'rad': rad1 - rad2,
                                  'rad percent': 100 * (rad1 - rad2) / rad2, 'dub': dub1 == dub2}

        if comparison['diffs'][_]['rad percent'] > 100:
            print(_, comparison['diffs'][_]['rad percent'], rad1, rad2, dub1, dub2)
    # Return the dictionary
    return comparison


def plot_comparison(comparison_dictionary, folder=None, csv_file=None):
    """
    Plotting the difference in radius and location by vertex and place a venn diagram of the
    differences.
    """

    loc_dist_y, rad_diff_y, rad_pcnt_y = [], [], []
    for ndxs in comparison_dictionary['diffs']:
        if comparison_dictionary['diffs'][ndxs]['loc'] > 1 or comparison_dictionary['diffs'][ndxs]['rad'] > 1:
            continue
        loc_dist_y.append(comparison_dictionary['diffs'][ndxs]['loc'])
        rad_diff_y.append(comparison_dictionary['diffs'][ndxs]['rad'])
        rad_pcnt_y.append(comparison_dictionary['diffs'][ndxs]['rad percent'])

    # Create the x values
    xs = np.arange(0, len(loc_dist_y))

    # Create figure and scatter plot
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.scatter(xs, rad_diff_y, label="Tangential Sphere\nRadius Difference", color='blue', alpha=0.3, s=2)
    ax.scatter(xs, loc_dist_y, label="Euclidean Distance", color='red', alpha=0.3, s=2)

    # Setting labels with specified font sizes
    ax.set_xlabel("Vertex", fontsize=20)
    ax.set_ylabel("Difference (e-13)", fontsize=20)
    ax.set_ylim([-0.0000000000001, 0.0000000000002])
    ax.set_xticks([])

    # Setting tick parameters directly
    ax.tick_params(axis='y', labelsize=20)

    ax.set_title(f"{comparison_dictionary['info']['dic1 name']} vs {comparison_dictionary['info']['dic2 name']} Vertex Comparison", fontsize=25)
    ax.legend(fontsize=15, loc='upper right', bbox_to_anchor=(1, 0.9), markerscale=5)

    # Create Venn diagram in top right
    venn_ax = fig.add_axes([0.2, 0.55, 0.3, 0.3])  # [left, bottom, width, height]
    num_shared = len(loc_dist_y)
    dic1_xtra, dic2_xtra = comparison_dictionary['info']['dic1 extra verts'], comparison_dictionary['info']['dic2 extra verts']
    tot_verts = num_shared + dic1_xtra + dic2_xtra
    venn = venn2(subsets=(np.log(dic1_xtra), np.log(dic2_xtra), np.log(num_shared)),
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
        subfolder = os.path.dirname(comparison_dictionary['pdb']['file'])
        sub_info = subfolder.split('_')
        cv, den, olap = sub_info[1], sub_info[3], sub_info[4]
        try:
            num = int(sub_info[-1])
        except ValueError:
            num = 0

        mean_diff_loc, std_loc_diff = np.mean([_ if _ < 0.5 else 0.5 for _ in loc_dist_y]), np.std([_ if _ < 0.5 else 0.5 for _ in loc_dist_y])
        mean_diff_pcnt_rad, std_diff_pcnt_rad = np.mean([abs(_) if abs(_) < 20 else 20 for _ in rad_pcnt_y]), np.std([abs(_) if abs(_) < 20 else 20 for _ in rad_pcnt_y])
        if std_diff_pcnt_rad > 10:
            print(max(rad_pcnt_y))
        line = [comparison_dictionary['info']['dic1 name'], comparison_dictionary['info']['dic2 name'], cv, den, olap,
                num, num_shared, (dic1_xtra, dic2_xtra), 100 * (1 - (num_shared / tot_verts)), mean_diff_loc,
                std_loc_diff, mean_diff_pcnt_rad, std_diff_pcnt_rad]
        if csv_file is not None:
            with open(csv_file, 'a') as my_file:
                csv_writer = csv.writer(my_file)
                csv_writer.writerow(line)
        else:
            print(line)
        new_folder = filedialog.askdirectory(title="Save the plot here")
        print("Saved at ", os.path.join(new_folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.png"))
        plt.savefig(os.path.join(new_folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.png"), format='png', dpi=1200)
        plt.show()
        # plt.show()
    else:
        subfolder = os.path.dirname(comparison_dictionary['pdb']['file'])
        plt.savefig(os.path.join(folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.svg"), format='svg')

        plt.savefig(os.path.join(folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.png"), format='png', dpi=1200)
        print(os.path.join(folder, f"{subfolder}_{comparison_dictionary['info']['dic1 name']}_{comparison_dictionary['info']['dic2 name']}.png"))
        plt.close(fig)


def run_it_all():
    # Get the folder that holds the subfolders
    folder = filedialog.askdirectory()

    choice_csv_folder = filedialog.askdirectory()
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
        # Get the pdb file
        try:
            pdb_file = [os.path.join(folder, subfolder, file) for file in os.listdir(os.path.join(folder, subfolder)) if file[-4:] == '.pdb'][0]
        except IndexError:
            continue
        if vpy_verts is None or vta_verts is None:
            continue
        # Compare the files now
        vpy, vta, pdb = get_information(vpy_fl=vpy_verts, vta_fl=vta_verts, pdb_fl=pdb_file)

        # Create the comparison combos
        

        # Compare the files
        comparison_dictionary = compare_vertices(vpy, vta, pdb, 'Vorpy', 'Voronota')
        # Plot the comparison
        plot_comparison(comparison_dictionary, csv_file=os.path.join(choice_csv_folder, 'all_info.csv'))



if __name__ == '__main__':

    # First open the folder with all the data
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    # vorpy, voronota, V, pdb_dict = get_information()
    # comp_dic = compare_vertices(vorpy, voronota, pdb_dict, 'Vorpy', 'Voronota')
    # plot_comparison(comp_dic)
    run_it_all()

