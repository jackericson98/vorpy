import csv
import os
from os import path
import tkinter as tk
from tkinter import filedialog
from scipy import stats
from matplotlib import pyplot as plt
from System.system import System
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from Data.Analyze.tools.compare.read_logs import read_logs
from System.sys_funcs.input.pdb import read_pdb_line
from System.sys_funcs.output.atoms import make_pdb_line
import numpy as np
from scipy.optimize import curve_fit


def get_files(folder):
    """Traverse the folder and extract required PDB and log files.
    Logs are identified based on subfolder names containing 'aw' or 'pow'.
    """
    pdb_file, aw_logs, pow_logs = None, None, None

    for root_dir, sub_dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.pdb') and 'atoms' not in file and 'diff' not in file:
                pdb_file = path.join(root_dir, file)

        # Check sub-subfolders for logs
        for sub_dir in sub_dirs:
            if 'aw' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        aw_logs = path.join(root_dir, sub_dir, file)
            elif 'pow' in sub_dir.lower():
                for file in os.listdir(path.join(root_dir, sub_dir)):
                    if file.endswith('.csv'):
                        pow_logs = path.join(root_dir, sub_dir, file)
    return pdb_file, aw_logs, pow_logs


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
sets = {}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
while True:
    folder = filedialog.askdirectory(title='Choose the folder mf')
    pdb, aw_logs, pow_logs = get_files(folder)


    # Take in the PDB
    # pdb = filedialog.askopenfilename(title='Get the pdb file mf')
    my_split_pdb = pdb.split('/')[-1]
    density = my_split_pdb.split('_')[3]
    cv = my_split_pdb.split('_')[1]
    # Get the power logs
    # pow_logs = filedialog.askopenfilename(title='Get the pow logs mf')
    # Get the Voronoit Logs
    # vor_logs = filedialog.askopenfilename(title='Get the vor logs mf')
    sets[density] = {'pdb': pdb, 'pow_logs': pow_logs, 'vor_logs': aw_logs}
    try:
        pow = read_logs2(pow_logs)
        aw = read_logs2(aw_logs)
        old = False
    except:
        pow = read_logs(pow_logs)
        aw = read_logs(aw_logs)
        old = True
    my_sys = System(pdb, simple=True)
    rads, diffs, names = [], [], []

    # Loop through the pdb file
    with open(pdb, 'r') as pdb_reader, open(pdb[:-4] + '_vor_diff_colored.pdb', 'w') as pdb_writer:
        mini, maxi = np.inf, -np.inf
        # Loop through the lines
        for i, line in enumerate(pdb_reader.readlines()):

            if i == 0:
                pdb_writer.write(line)
                continue
            npl = read_pdb_line(line)
            # Check that the number is actually in the dataframe
            try:
                if not old:
                    # Get the pow atom and the vor atom
                    pow_atom, vor_atom = pow['atoms'].loc[pow['atoms']['Index'] == i - 1].to_dict('records')[0], aw['atoms'].loc[aw['atoms']['Index'] == i - 1].to_dict('records')[0]

                    # Calculate the difference in volume
                    vol_diff = (pow_atom['Volume'] - vor_atom['Volume']) / vor_atom['Volume']
                else:
                    # Get the pow atom and the vor atom
                    pow_atom, vor_atom = pow['atoms'].loc[pow['atoms']['num'] == i - 1].to_dict('records')[0], aw['atoms'].loc[aw['atoms']['num'] == i - 1].to_dict('records')[0]

                    # Calculate the difference in volume
                    vol_diff = (pow_atom['volume'] - vor_atom['volume']) / vor_atom['volume']
                # Check for crazy volume difference and trigger a volume difference
                if vol_diff >= 30:
                    vol_diff = 30
                    print(npl['atom_serial_number'], 'Off by {} %'.format(100 * vol_diff))
                diffs.append(vol_diff)
                rads.append(npl['temperature_factor'])
                names.append(npl['atom_name'])
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

    # Write the set code
    with open(pdb[:-4] + '_set_diff.txt', 'w') as set_color:
        # Write the first line
        set_color.write('spectrum q, blue_white_red, minimum={}, maximum={}\n'.format(-maxi, maxi))
        # Select the group to not be colored
        set_color.write('color white, chain Z')


    # def func(x, a, b, c):
    #     return a * np.exp(-b * x) + c
    #
    # # Plot the radius to difference values
    # plt.scatter(rads, [100 * _ for _ in diffs], s=2, alpha=0.3)
    # plt.plot([min(rads), max(rads)], [0, 0], c='k')
    # popt, pcov = curve_fit(func, np.array(rads), np.array(diffs))
    # # plt.text(s='y = {:.2f} * exp(-{:.2f} * x) + {:.2f}'.format(*popt), x=1.5, y=1.5, font=dict(size=10))
    # plt.plot(rads, [100 * _ for _ in func(np.array(rads), *popt)], label=density)
    #
    # plot_another = input('Plot Another?    ')
    # if plot_another == 'n':
    #     break

# plt.xlabel('Ball Radius', fontdict=dict(size=25))
# plt.ylabel('% Difference', fontdict=dict(size=25))
# # plt.ti
# plt.xticks(font=dict(size=30))
# plt.yticks(font=dict(size=30))
# plt.title('% Difference by Radii\n(Non-Overlapping, CV {})'.format(cv), font=dict(size=20))
# plt.tick_params(axis='both', width=2, length=12)
# plt.legend(fontsize=15)
#
# plt.tight_layout()
# plt.show()


# import csv
# import os
# from os import path
# import tkinter as tk
# from tkinter import filedialog
# from scipy import stats
# from matplotlib import pyplot as plt
# from System.system import System
# from Data.Analyze.tools.compare.read_logs2 import read_logs2
# from Data.Analyze.tools.compare.read_logs import read_logs
# from System.sys_funcs.input.pdb import read_pdb_line
# from System.sys_funcs.output.atoms import make_pdb_line
# import numpy as np
# from scipy.optimize import curve_fit
#
#
# root = tk.Tk()
# root.withdraw()
# root.wm_attributes('-topmost', 1)
# sets = {}
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
# while True:
#     folder = filedialog.askdirectory(title='Choose the folder mf')
#     for rroot, directys, files in os.walk(folder):
#         for file in files:
#             if file[-3:] == 'pdb' and 'atoms' not in file and 'diff' not in file:
#                 pdb = rroot + '/' + file
#             if file[-3:] == 'csv' and 'aw' in file:
#                 aw_logs = rroot + '/' + file
#             if file[-3:] == 'csv' and 'pow' in file:
#                 pow_logs = rroot + '/' + file
#
#
#     # Take in the PDB
#     # pdb = filedialog.askopenfilename(title='Get the pdb file mf')
#     my_split_pdb = pdb.split('/')[-1]
#     density = my_split_pdb.split('_')[3]
#     cv = my_split_pdb.split('_')[1]
#     # Get the power logs
#     # pow_logs = filedialog.askopenfilename(title='Get the pow logs mf')
#     # Get the Voronoit Logs
#     # vor_logs = filedialog.askopenfilename(title='Get the vor logs mf')
#     sets[density] = {'pdb': pdb, 'pow_logs': pow_logs, 'vor_logs': aw_logs}
#     try:
#         pow = read_logs2(pow_logs)
#         aw = read_logs2(aw_logs)
#         old = False
#     except:
#         pow = read_logs(pow_logs)
#         aw = read_logs(aw_logs)
#         old = True
#     my_sys = System(pdb, simple=True)
#     rads, diffs, log_diffs, names = [], [], [], []
#
#     # Loop through the pdb file
#     with (open(pdb, 'r') as pdb_reader, open(pdb[:-4] + '_vor_diff_colored.pdb', 'w') as pdb_writer, open(pdb[:-4] + '_vor_diff_colored_log.pdb', 'w') as pdb_log_writer):
#         mini, maxi = np.inf, -np.inf
#         # Loop through the lines
#         for i, line in enumerate(pdb_reader.readlines()):
#
#             if i == 0:
#                 pdb_writer.write(line)
#                 continue
#             npl = read_pdb_line(line)
#             # Check that the number is actually in the dataframe
#             try:
#                 if not old:
#                     # Get the pow atom and the vor atom
#                     pow_atom, vor_atom = pow['atoms'].loc[pow['atoms']['Index'] == i - 1].to_dict('records')[0], aw['atoms'].loc[aw['atoms']['Index'] == i - 1].to_dict('records')[0]
#
#                     # Calculate the difference in volume
#                     vol_diff = (pow_atom['Volume'] - vor_atom['Volume']) / vor_atom['Volume']
#                 else:
#                     # Get the pow atom and the vor atom
#                     pow_atom, vor_atom = pow['atoms'].loc[pow['atoms']['num'] == i - 1].to_dict('records')[0], aw['atoms'].loc[aw['atoms']['num'] == i - 1].to_dict('records')[0]
#
#                     # Calculate the difference in volume
#                     vol_diff = (pow_atom['volume'] - vor_atom['volume']) / vor_atom['volume']
#                 # Check for crazy volume difference and trigger a volume difference
#                 if vol_diff >= 30:
#                     print(npl['atom_serial_number'], 'Off by {} %'.format(100 * vol_diff))
#                     _ = float('')
#                 log_diffs.append(np.log(vol_diff))
#                 diffs.append(vol_diff)
#                 rads.append(npl['temperature_factor'])
#                 names.append(npl['atom_name'])
#                 if vol_diff < mini:
#                     mini = vol_diff
#                 if vol_diff > maxi:
#                     maxi = vol_diff
#
#                 new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'],
#                                              res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
#                                              x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
#                                              z=float(npl['z_coordinate']), occ=vol_diff,
#                                              tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
#                 new_log_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'],
#                                                  res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
#                                                  x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
#                                                  z=float(npl['z_coordinate']), occ=np.log(vol_diff),
#                                                  tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
#                 pdb_writer.write(new_pdb_line)
#                 pdb_log_writer.write(new_log_pdb_line)
#
#             except IndexError:
#                 new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'], chain='Z',
#                                              res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
#                                              x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
#                                              z=float(npl['z_coordinate']), occ=0.0,
#                                              tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
#                 pdb_writer.write(new_pdb_line)
#                 pdb_log_writer.write(new_pdb_line)
#             except ValueError:
#                 new_pdb_line = make_pdb_line(ser_num=int(npl['atom_serial_number']), name=npl['atom_name'], chain='Z',
#                                              res_name=npl['residue_name'], res_seq=int(npl['residue_sequence_number']),
#                                              x=float(npl['x_coordinate']), y=float(npl['y_coordinate']),
#                                              z=float(npl['z_coordinate']), occ=0.0,
#                                              tfact=float(npl['temperature_factor']), elem=npl['element_symbol'])
#                 pdb_writer.write(new_pdb_line)
#                 pdb_log_writer.write(new_pdb_line)
#     # Write the set code
#     with open(pdb[:-4] + '_set_diff.txt', 'w') as set_color:
#         # Write the first line
#         set_color.write('spectrum q, blue_white_red, minimum={}, maximum={}\n'.format(-maxi, maxi))
#         # Select the group to not be colored
#         set_color.write('color white, chain Z')
#     # Write the set code
#     with open(pdb[:-4] + '_set_log_diff.txt', 'w') as set_color:
#         # Write the first line
#         set_color.write('spectrum q, blue_white_red, minimum={}, maximum={}\n'.format(-np.log(maxi), np.log(maxi)))
#         # Select the group to not be colored
#         set_color.write('color white, chain Z')
#
#
#     def func(x, a, b, c):
#         return a * np.exp(-b * x) + c
#
#     # Plot the radius to difference values
#     ax1.scatter(rads, [100 * _ for _ in diffs], s=2, alpha=0.3)
#     ax1.plot([min(rads), max(rads)], [0, 0], c='k')
#     popt, pcov = curve_fit(func, np.array(rads), np.array(diffs))
#     # plt.text(s='y = {:.2f} * exp(-{:.2f} * x) + {:.2f}'.format(*popt), x=1.5, y=1.5, font=dict(size=10))
#     ax1.plot(rads, [100 * _ for _ in func(np.array(rads), *popt)], label=density)
#     # Plot the radius to difference values
#     ax2.scatter(rads, [100 * _ for _ in log_diffs], s=2, alpha=0.3)
#     ax2.plot([min(rads), max(rads)], [0, 0], c='k')
#     popt, pcov = curve_fit(func, np.array(rads), np.array(log_diffs))
#     # plt.text(s='y = {:.2f} * exp(-{:.2f} * x) + {:.2f}'.format(*popt), x=1.5, y=1.5, font=dict(size=10))
#     ax2.plot(rads, [100 * _ for _ in func(np.array(rads), *popt)], label=density)
#
#     plot_another = input('Plot Another?    ')
#     if plot_another == 'n':
#         break
#
# plt.xlabel('Ball Radius', fontdict=dict(size=25))
# plt.ylabel('% Difference', fontdict=dict(size=25))
# # plt.ti
# plt.xticks(font=dict(size=30))
# plt.yticks(font=dict(size=30))
# plt.title('% Difference by Radii\n(Non-Overlapping, CV {})'.format(cv), font=dict(size=20))
# plt.tick_params(axis='both', width=2, length=12)
# plt.legend(fontsize=15)
#
# plt.tight_layout()
# plt.show()


