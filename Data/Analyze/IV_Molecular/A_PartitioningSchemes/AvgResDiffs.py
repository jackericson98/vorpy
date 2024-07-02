import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.get_res_data import residue_data


if __name__ == '__main__':
    # Get the dropbox folder
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    folder = filedialog.askdirectory()
    # Get the systems in the designated folder
    systems = []
    for root, directory, files in os.walk(folder):
        for file in files:
            if file[-3:] == 'pdb':
                my_sys = System(file=folder + '/' + file)
                my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
                systems.append(my_sys)

    # Sort atoms by number of atoms
    num_atoms = [len(_.atoms) for _ in systems]
    systems = [x for _, x in sorted(zip(num_atoms, systems))]
    # Create the logs dictionary
    my_sys_names = [__.name for __ in systems]

    # Create the log file name dictionary
    my_log_files = {_: {__: folder + '/' + _ + '_{}_logs.csv'.format(__) for __ in {'vor', 'pow', 'del'}}
                    for _ in my_sys_names}

    # Create the log dictionary
    my_logs = {}
    # print('39: my_log_files = {}'.format(my_log_files))

    # Get the log information
    (pow_vol_avg_diff, del_vol_avg_diff, pow_vol_se, del_vol_se, pow_sa_avg_diff, del_sa_avg_diff, pow_sa_se,
     del_sa_se) = [], [], [], [], [], [], [], []
    for system in systems:
        # print('45: system = {}'.format(system.name))
        pow_vols, del_vols, pow_sas, del_sas = [], [], [], []
        # Get the values from the residue function
        vor_reses = residue_data(system, read_logs(my_log_files[system.name]['vor']), get_all=True)
        # print('49: vor_reses = {}'.format(vor_reses))
        pow_reses = residue_data(system, read_logs(my_log_files[system.name]['pow']), get_all=True)
        # print('51: pow_reses = {}'.format(pow_reses))
        del_reses = residue_data(system, read_logs(my_log_files[system.name]['del']), get_all=True)
        # print('53: del_reses = {}'.format(del_reses))
        # Find the percent differences by residue
        # Classification level
        for _ in vor_reses:
            # Sub class level
            for __ in vor_reses[_]:
                # Res_seq level
                for ___ in vor_reses[_][__]:
                    if vor_reses[_][__][___] == {}:
                        continue
                    if vor_reses[_][__][___]['vol'] == 0 or vor_reses[_][__][___]['sa'] == 0:
                        continue
                    pow_vols.append(abs(vor_reses[_][__][___]['vol'] - pow_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol'])
                    del_vols.append(abs(vor_reses[_][__][___]['vol'] - del_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol'])
                    pow_sas.append(
                        abs(vor_reses[_][__][___]['sa'] - pow_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa'])
                    del_sas.append(
                        abs(vor_reses[_][__][___]['sa'] - del_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa'])
        # Get the averages
        print('AvgResDiffs 73: ', system.name)
        print('AvgResDiffs 74: ', 100 * sum(pow_vols)/len(pow_vols), 100 * np.std(pow_vols)/np.sqrt(len(pow_vols)))
        print('AvgResDiffs 75: ', 100 * sum(del_vols)/len(del_vols), 100 * np.std(pow_vols)/np.sqrt(len(pow_vols)))
        print('AvgResDiffs 76: ', 100 * sum(pow_sas)/len(pow_sas), 100 * np.std(pow_sas)/np.sqrt(len(pow_sas)))
        print('AvgResDiffs 77: ', 100 * sum(del_sas)/len(del_sas), 100 * (np.std(del_sas)/np.sqrt(len(del_sas))))
        # Get the standard Errors
        pow_vol_avg_diff.append(100 * sum(pow_vols)/len(pow_vols))
        del_vol_avg_diff.append(100 * sum(del_vols)/len(del_vols))
        pow_sa_avg_diff.append(100 * sum(pow_sas)/len(pow_sas))
        del_sa_avg_diff.append(100 * sum(del_sas)/len(del_sas))
        # Get the standard Errors
        pow_vol_se.append(np.std(pow_vols)/np.sqrt(len(pow_vols)))
        del_vol_se.append(np.std(del_vols)/np.sqrt(len(del_vols)))
        pow_sa_se.append(np.std(pow_sas)/np.sqrt(len(pow_sas)))
        del_sa_se.append(np.std(del_sas)/np.sqrt(len(del_sas)))

    # Create the dictionary for converting the labels
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP', 'BSA': 'BSA', '1BNA': '1BNA',
                     'DB1976': 'DB1976'}
                    [_] for _ in my_sys_names]
    # Set the label codes
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    # Create the bar graph
    bar([pow_vol_avg_diff, del_vol_avg_diff], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Volume Difference',
        errors=[pow_vol_se, del_vol_se])
    # Create the bar graph
    bar([pow_sa_avg_diff, del_sa_avg_diff], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Surface Area Difference',
        errors=[pow_sa_se, del_sa_se])

