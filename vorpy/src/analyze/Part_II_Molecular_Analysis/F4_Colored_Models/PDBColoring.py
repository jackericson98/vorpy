import os
import tkinter as tk
from tkinter import filedialog

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)


from vorpy.src.output.pdb import make_pdb_line
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.inputs.pdb import read_pdb_line

def get_atom_dict(logs_file):
    atom_dict = {}
    logs = read_logs2(logs_file)
    # Get the atom dictionary
    for i, _ in logs['atoms'].iterrows():
        atom_dict[(_['Name'], _['Residue'], _['Residue Sequence'])] = (_['Maximum Mean Curvature'], _['Maximum Gaussian Curvature'], _['Average Mean Surface Curvature'], _['Average Gaussian Surface Curvature'])
    return atom_dict


def color_pdb_by_curvature(pdb, values, output_folder=None):
    if output_folder is None:
        output_folder = os.path.dirname(pdb)
    output_pdb = output_folder + '/' + os.path.basename(pdb)[:-4] + '_colored_by_curvature.pdb'
    with open(pdb, 'r') as read_pdb, open(output_pdb, 'w') as write_pdb:
        for line in read_pdb:
            if line[:4] == 'ATOM':
                a = read_pdb_line(line)
                if a['residue_name'] == 'SOL' or a['residue_name'] == 'CL':
                    write_pdb.write(line)
                    continue
                if (a['atom_name'], a['residue_name'], a['residue_sequence_number']) in values:
                    tfact = values[(a['atom_name'], a['residue_name'], a['residue_sequence_number'])][2]
                else:
                    tfact = 0
                
                # Handle empty chain identifier
                chain_id = a['chain_identifier']
                if not chain_id or chain_id.strip() == "":
                    chain_id = "A"  # Default to chain A if no chain is specified
                
                write_pdb.write(make_pdb_line(ser_num=a['atom_serial_number'], name=a['atom_name'], res_name=a['residue_name'], chain=chain_id,
                                       res_seq=a['residue_sequence_number'], x=a['x_coordinate'], y=a['y_coordinate'], z=a['z_coordinate'], tfact=tfact, elem=a['element_symbol']))
    
def write_pymol_script(output_folder, name, low_val=0, high_val=0):
    with open(output_folder + '/' + name + '_set_colors.pml', 'w') as pymol_script:
        pymol_script.write('spectrum b, red_green, minimum=0, maximum={}'.format(high_val))

if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    pdb_file = filedialog.askopenfilename(title='Choose PDB File')
    logs_file = filedialog.askopenfilename(title='Choose Logs File')
    output_folder = filedialog.askdirectory(title='Choose Output Folder')

    atom_dict = get_atom_dict(logs_file)
    max_val = max(atom_dict.values())
    color_pdb_by_curvature(pdb_file, atom_dict, output_folder)
    write_pymol_script(output_folder, name=os.path.basename(pdb_file)[:-4], low_val=0, high_val=max_val)
    print('Done!')

