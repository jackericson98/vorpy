import os
import tkinter as tk
from tkinter import filedialog
from vorpy.src.system.system import System
from System.Group.group import Group
from System.sys_funcs.output.atoms import make_pdb_line
from Data.Analyze.tools.compare.read_logs import read_logs

# Step 1 get the systems
# Step 2 get the logs
# Step three output pdb with curvature at b-factor
# Output pymol script for ease of use


def assign_atom_color(system, values, val='curv', directory=None):
    if directory is None:
        directory = './'
    with open(directory + system.name + '_colored_by_' + val + '.pdb', 'w') as my_pdb:
        for i, a in system.atoms.iterrows():
            # Get the location string
            x, y, z = a['loc']
            # Get the information from the atom in writable format
            try:
                tfact = values[a['num']]
            except KeyError:
                tfact = 0
            chain_name = a['chn'].name
            if a['chn'].name == 'SOL':
                chain_name = 'Z'
                tfact = 0
            # Write the atom information
            my_pdb.write(make_pdb_line(ser_num=i, name=a['name'], res_name=a['res'].name, chain=chain_name,
                                       res_seq=a['res_seq'], x=x, y=y, z=z, tfact=tfact, elem=a['element']))


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory(title='Choose Logs Pdbs Folder')
    # Create the systems
    systems = []
    for root, dir, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb':
                my_sys = System(file=folder + '/' + file)
                my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
                systems.append(my_sys)
    # Sort atoms by number of atoms
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]

    # Set the output folder
    output_folder = filedialog.askdirectory(title='Choose Output Folder') + '/'
    # Create the outputs by system
    my_maxes = []
    for my_sys in systems:
        # Read the logs
        my_log_vals = read_logs(folder + '/' + my_sys.name + '_vor_logs.csv', return_dict=True)
        my_maxes.append(max([_['max curv'] for _ in my_log_vals['atoms']]))
        assign_atom_color(my_sys, {_['num']: _['max curv'] for _ in my_log_vals['atoms']}, directory=output_folder)
    # create the coloring script
    with open(output_folder + 'set_colors.pml', 'w') as pymol_script:
        pymol_script.write('spectrum b, green_yellow_red, minimum=0, maximum={}'.format(max(my_maxes)))
