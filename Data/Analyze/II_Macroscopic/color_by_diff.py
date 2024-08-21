import csv
from os import path
import tkinter as tk
from tkinter import filedialog
from scipy import stats
from matplotlib import pyplot as plt

from Data.Analyze.tools.compare.read_logs import read_logs
from System.sys_funcs.input.pdb import read_pdb_line
from System.sys_funcs.output.atoms import make_pdb_line
import numpy as np


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)

# Take in the PDB
pdb = filedialog.askopenfilename(title='Get the pdb file mf')
# Get the power logs
pow_logs = filedialog.askopenfilename(title='Get the pow logs mf')
# Get the Voronoit Logs
vor_logs = filedialog.askopenfilename(title='Get the vor logs mf')

# def make_logs_dict(logs):
#     # Open the file
#     with open(logs, 'r') as my_logs:
#         # Make the csv_reader
#         lgs = csv.reader(my_logs)
#         # Set up the reading
#         reading = True
#         # Loop through the rows
#         for i, row in enumerate(lgs):
#             # Check if it is an atom line
#             if row[0] == 'Atoms':

pow_ = read_logs(pow_logs)
vor_ = read_logs(vor_logs)

rads, diffs = [], []
# Loop through the pdb file
with (open(pdb, 'r') as pdb_reader, open(pdb[:-4] + '_vor_diff_colored.pdb', 'w') as pdb_writer):
    mini, maxi = np.inf, -np.inf
    # Loop through the lines
    for i, line in enumerate(pdb_reader.readlines()):

        if i == 0:
            pdb_writer.write(line)
            continue
        npl = read_pdb_line(line)
        # Check that the number is actually in the dataframe
        try:
            # Get the pow atom and the vor atom
            pow_atom, vor_atom = pow_['atoms'].loc[pow_['atoms']['num'] == i - 1].to_dict('records')[0], vor_['atoms'].loc[vor_['atoms']['num'] == i - 1].to_dict('records')[0]

            # Calculate the difference in volume
            vol_diff = (pow_atom['volume'] - vor_atom['volume']) / vor_atom['volume']

            # Check for crazy volume difference and trigger a volume difference
            if vol_diff >= 10:
                float('poo')
                print("not triggering a value error")
            diffs.append(vol_diff)
            rads.append(npl['temperature_factor'])
            if vol_diff < mini:
                mini = vol_diff
            if vol_diff > maxi:
                maxi = vol_diff

            new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'],
                                         res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
                                         x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
                                         z=float(npl['z_coordinate']), occ=vol_diff,
                                         tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
            pdb_writer.write(new_pdb_line)

        except IndexError:
            new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'], chain='Z',
                                         res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
                                         x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
                                         z=float(npl['z_coordinate']), occ=0.0,
                                         tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
            pdb_writer.write(new_pdb_line)
        except ValueError:
            new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'], chain='Z',
                                         res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
                                         x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
                                         z=float(npl['z_coordinate']), occ=0.0,
                                         tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
            pdb_writer.write(new_pdb_line)
            print(npl['atom_serial_number'])

# Write the set code
with open(pdb[:-4] + '_set_diff.txt', 'w') as set_color:
    # Write the first line
    set_color.write('spectrum q, green_yellow_red, minimum={}, maximum={}\n'.format(mini, maxi))
    # Select the group to not be colored
    set_color.write('color white, chain Z')

# Plot the radius to difference values
plt.scatter(rads, diffs)
plt.plot([min(rads), max(rads)], [0, 0])
slope, intercept, r_value, p_value, std_err = stats.linregress(rads, diffs)
plt.plot([min(rads), max(rads)], [min(rads) * slope + intercept, max(rads) * slope + intercept])
plt.show()
