import os
import numpy as np
import tkinter as tk
from tkinter import filedialog

import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.system.system import System
from vorpy.src.group.group import Group
from vorpy.src.analyze.tools.plot_templates.bar import bar
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def get_res_sa(logs):
    # Create the dictionary for the residue data
    res_data = {}
    # Create a dictionary to find the residue key for a given surface
    res_key_dict = {}
    # Loop through and get all of the residues
    for i, atom in logs['atoms'].iterrows():
        # Create the residue key
        res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
        # Check to see if the residue key is in the dictionary
        if res_key not in res_data:
            res_data[res_key] = [atom['Index']]
        else:
            res_data[res_key].append(atom['Index'])
        # Add the residue key to the dictionary
        res_key_dict[atom['Index']] = res_key
    # Create the residue surface area dictionary
    res_sa_dict = {}
    # Loop through the surfaces and add the surface area if they are outer surfaces
    for i, surf in logs['surfs'].iterrows():
        # Get the first ball in the surface
        ball_1, ball_2 = surf['Balls']
        # Get the residue key for the first ball
        if ball_1 in res_key_dict:
            res_key = res_key_dict[ball_1]
        elif ball_2 in res_key_dict:
            res_key = res_key_dict[ball_2]
        else:
            print(f"Ball {ball_1} or {ball_2} not in the residue key dictionary")
            continue
        # Check if the res_sa_dict has the residue key
        if res_key not in res_sa_dict:
            res_sa_dict[res_key] = 0
        # Now make sure that only one of the balls is in the res_data dictionary
        if ball_1 in res_data[res_key] and ball_2 in res_data[res_key]:
            continue
        else:
            res_sa_dict[res_key] += surf['Surface Area']


    # Return the dictionary
    return res_sa_dict


def get_res_data(folder=None, exclude_keys=[], get_sa=False):
    """
    Function that gets the residue data for the given systems
    """
    if folder is None:
        # Get the dropbox folder
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory()
    # Create the dictionary for the system residue averages
    sys_res_data = {}
    # Go through the folder and get the systems
    for subfolder in os.listdir(folder):
        my_key = subfolder.split('_')[0]
        if my_key in exclude_keys:
            continue
        # get the logs files addresses
        aw_logs = os.path.join(folder, subfolder, 'aw_logs.csv')
        pow_logs = os.path.join(folder, subfolder, 'pow_logs.csv')
        prm_logs = os.path.join(folder, subfolder, 'prm_logs.csv')
        # Read the logs to get the residue data
        aw_logs = read_logs2(aw_logs, all_=False, balls=True, surfs=True)
        pow_logs = read_logs2(pow_logs, all_=False, balls=True, surfs=True)
        prm_logs = read_logs2(prm_logs, all_=False, balls=True, surfs=True)
        # Check to see if we want to get the surface area
        if get_sa:
            aw_res_sa = get_res_sa(aw_logs)
            pow_res_sa = get_res_sa(pow_logs)
            prm_res_sa = get_res_sa(prm_logs)
        # Create the dictionary for the residue data
        res_data = {}
        aw_vols, pow_vols, prm_vols = [], [], []
        # Get the residue data
        for i, atom in aw_logs['atoms'].iterrows():
            # Get the pow and prm atoms
            pow_atom = pow_logs['atoms'].loc[pow_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
            prm_atom = prm_logs['atoms'].loc[prm_logs['atoms']['Index'] == atom['Index']].to_dict(orient='records')[0]
            # Create the residue key
            res_key = f"{atom['Residue']}_{atom['Residue Sequence']}"
            # Check to see if the residue key is in the dictionary
            if res_key not in res_data:
                res_data[res_key] = {'aw': {'vol': atom['Volume'], 'sa': 0}, 'pow': {'vol': pow_atom['Volume'], 'sa': 0}, 'prm': {'vol': prm_atom['Volume'], 'sa': 0}}
            else:
                res_data[res_key]['aw']['vol'] += atom['Volume']
                res_data[res_key]['pow']['vol'] += pow_atom['Volume']
                res_data[res_key]['prm']['vol'] += prm_atom['Volume']
                aw_vols.append(res_data[res_key]['aw']['vol'])
                pow_vols.append(res_data[res_key]['pow']['vol'])
                prm_vols.append(res_data[res_key]['prm']['vol'])

        aw_sas, pow_sas, prm_sas = [], [], []
        # Add the surface areas to the residues
        for res_key in res_data:
            res_data[res_key]['aw']['sa'] = aw_res_sa[res_key]
            res_data[res_key]['pow']['sa'] = pow_res_sa[res_key]
            res_data[res_key]['prm']['sa'] = prm_res_sa[res_key]
            aw_sas.append(res_data[res_key]['aw']['sa'])
            pow_sas.append(res_data[res_key]['pow']['sa'])
            prm_sas.append(res_data[res_key]['prm']['sa'])

        # Now get the averages and standard deviations for thew systems sa and vol
        sys_res_data[my_key] = {'aw': {'avg vol': sum(aw_vols)/len(aw_vols), 'avg sa': sum(aw_sas)/len(aw_sas), 'se vol': np.std(aw_vols)/np.sqrt(len(aw_vols)), 'se sa': np.std(aw_sas)/np.sqrt(len(aw_sas))}, 
                                'pow': {'avg vol': sum(pow_vols)/len(pow_vols), 'avg sa': sum(pow_sas)/len(pow_sas), 'se vol': np.std(pow_vols)/np.sqrt(len(pow_vols)), 'se sa': np.std(pow_sas)/np.sqrt(len(pow_sas))}, 
                                'prm': {'avg vol': sum(prm_vols)/len(prm_vols), 'avg sa': sum(prm_sas)/len(prm_sas), 'se vol': np.std(prm_vols)/np.sqrt(len(prm_vols)), 'se sa': np.std(prm_sas)/np.sqrt(len(prm_sas))}}

    # Return the dictionary
    return sys_res_data

def plot_data(sys_res_data):
    #     # Create the bar graph
    bar([[sys_res_data[key]['aw']['avg vol'] for key in sys_res_data], [sys_res_data[key]['pow']['avg vol'] for key in sys_res_data]], x_names=list(sys_res_data.keys()), legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Volume Difference',
        errors=[[sys_res_data[key]['aw']['se vol'] for key in sys_res_data], [sys_res_data[key]['pow']['se vol'] for key in sys_res_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)
    # Create the bar graph
    bar([[sys_res_data[key]['aw']['avg sa'] for key in sys_res_data], [sys_res_data[key]['pow']['avg sa'] for key in sys_res_data]], x_names=list(sys_res_data.keys()), legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Surface Area Difference',
        errors=[[sys_res_data[key]['aw']['se sa'] for key in sys_res_data], [sys_res_data[key]['pow']['se sa'] for key in sys_res_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)


if __name__ == '__main__':
    # Get the data
    sys_res_data = get_res_data(get_sa=True, exclude_keys=['A', 'B', 'C'])
    # Plot the data
    plot_data(sys_res_data)






# # pre_made_data = [
# #     ['cambrin', 2.657358947567952, 0.1701522540304929, 2.59465391164789, 0.1701522540304929, 0.7071402007950848, 0.059791991731881135, 2.3096433988252345, 0.21297727309940268],
# #     ['hairpin', 2.120015342765748, 0.17989523306316388, 4.207094792681904, 0.17989523306316388, 0.3151423480069726, 0.050484980940870186, 0.587710293326916, 0.09756079607665612],
# #     ['p53tet', 2.9331662712893984, 0.18342810870548817, 2.7747918309068043, 0.18342810870548817, 0.6593449473781604, 0.10422962591710631, 2.370253686067341, 0.25677219904214393],
# #     ['pl_complex', 2.7747879360584684, 0.07901575226030394, 2.8010158048936318, 0.07901575226030394, 0.62219526460963, 0.03121395618438131, 2.422886871774085, 0.11430094978326413],
# #     ['streptavidin', 3.0522725810134834, 0.1238965156383565, 2.174625260864225, 0.1238965156383565, 0.4774545212283423, 0.0397510499787414, 2.072412116883273, 0.1310601195844641],
# #     ['hammerhead', 6.576136641455463, 1.3355495236483024, 8.36539089536273, 1.3355495236483024, 7.027173386404598, 1.337333613171973, 5.30256130664321, 2.34512218840759],
# #     ['NCP', 8.826521657702823, 0.9350257210684129, 8.618442644565837, 0.9350257210684129, 1.7301587766519828, 0.19436185835523279, 7.435837930712498, 1.0235150388051224],
# #     ['BSA', 2.8614722131257735, 0.5170327576040869, 8.53612886017547, 1.5170327576040869, 0.6362188314154772, 0.2564681508906418, 1.47362166167768, 0.16446425197757493]
# # ]


# if __name__ == '__main__':
#     # Get the dropbox folder
#     root = tk.Tk()
#     root.withdraw()
#     root.wm_attributes('-topmost', 1)
#     folder = filedialog.askdirectory()
#     # Get the systems in the designated folder
#     systems = []
#     for root, directory, files in os.walk(folder):
#         for file in files:
#             if file[-3:] == 'pdb':
#                 my_sys = System(file=folder + '/' + file)
#                 my_sys.groups = [Group(sys=my_sys, residues=my_sys.residues)]
#                 systems.append(my_sys)

#     # Sort atoms by number of atoms
#     num_atoms = [len(_.atoms) for _ in systems]
#     systems = [x for _, x in sorted(zip(num_atoms, systems))]
#     # Create the logs dictionary
#     my_sys_names = [__.name for __ in systems]

#     # Create the log file name dictionary
#     my_log_files = {_: {__: folder + '/' + _ + '_{}_logs.csv'.format(__) for __ in {'vor', 'pow', 'del'}}
#                     for _ in my_sys_names}

#     # Create the log dictionary
#     my_logs = {}
#     # print('39: my_log_files = {}'.format(my_log_files))

#     # Get the log information
#     (pow_vol_avg_diff, del_vol_avg_diff, pow_vol_se, del_vol_se, pow_sa_avg_diff, del_sa_avg_diff, pow_sa_se,
#      del_sa_se) = [], [], [], [], [], [], [], []
#     for system in systems:
#         print(system.name)
#         vor_out, vor_in = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
#         if os.path.exists(folder + '/res_data/{}_vor_res.csv'.format(system.name)):
#             vor_in, vor_out = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
#         pow_out, pow_in = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
#         if os.path.exists(folder + '/res_data/{}_pow_res.csv'.format(system.name)):
#             pow_in, pow_out = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
#         del_out, del_in = folder + '/res_data/{}_del_res.csv'.format(system.name), None
#         if os.path.exists(folder + '/res_data/{}_del_res.csv'.format(system.name)):
#             del_in, del_out = folder + '/res_data/{}_del_res.csv'.format(system.name), None
#         # print('45: system = {}'.format(system.name))
#         pow_vols, del_vols, pow_sas, del_sas = [], [], [], []
#         # Get the values from the residue function
#         vor_reses = residue_data(system, read_logs(my_log_files[system.name]['vor']), get_all=True, read_file=vor_in, output_file=vor_out)
#         # print('49: vor_reses = {}'.format(vor_reses))
#         pow_reses = residue_data(system, read_logs(my_log_files[system.name]['pow']), get_all=True, read_file=pow_in, output_file=pow_out)
#         # print('51: pow_reses = {}'.format(pow_reses))
#         del_reses = residue_data(system, read_logs(my_log_files[system.name]['del']), get_all=True, read_file=del_in, output_file=del_out)
#         # print('53: del_reses = {}'.format(del_reses))
#         # Find the percent differences by residue
#         # Classification level
#         for _ in vor_reses:
#             # Sub class level
#             for __ in vor_reses[_]:
#                 # Res_seq level
#                 for ___ in vor_reses[_][__]:
#                     if vor_reses[_][__][___] == {}:
#                         continue
#                     if vor_reses[_][__][___]['vol'] == 0 or vor_reses[_][__][___]['sa'] == 0:
#                         continue
#                     pow_vol_diff = (pow_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
#                     del_vol_diff = (pow_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
#                     pow_vols.append(pow_vol_diff)
#                     del_vols.append(del_vol_diff)
#                     pow_sa_diff = (vor_reses[_][__][___]['sa'] - pow_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa']
#                     del_sa_diff = (vor_reses[_][__][___]['sa'] - del_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa']
#                     pow_sas.append(pow_sa_diff)
#                     del_sas.append(del_sa_diff)
#                     if abs(del_vol_diff) > 100:
#                         print(__, ___, del_vol_diff, del_sa_diff)
#         # # Get the averages
#         # print('AvgResDiffs 73: ', system.name)
#         # print('AvgResDiffs 74: ', 100 * sum(pow_vols)/len(pow_vols), 100 * np.std(pow_vols)/np.sqrt(len(pow_vols)))
#         # print('AvgResDiffs 75: ', 100 * sum(del_vols)/len(del_vols), 100 * np.std(del_vols)/np.sqrt(len(del_vols)))
#         # print('AvgResDiffs 76: ', 100 * sum(pow_sas)/len(pow_sas), 100 * np.std(pow_sas)/np.sqrt(len(pow_sas)))
#         # print('AvgResDiffs 77: ', 100 * sum(del_sas)/len(del_sas), 100 * (np.std(del_sas)/np.sqrt(len(del_sas))))


#         # Get the standard Errors
#         pow_vol_avg_diff.append(100 * sum(pow_vols)/len(pow_vols))
#         del_vol_avg_diff.append(100 * sum(del_vols)/len(del_vols))
#         pow_sa_avg_diff.append(100 * sum(pow_sas)/len(pow_sas))
#         del_sa_avg_diff.append(100 * sum(del_sas)/len(del_sas))
#         # Get the standard Errors
#         pow_vol_se.append(np.std(pow_vols)/np.sqrt(len(pow_vols)))
#         del_vol_se.append(np.std(del_vols)/np.sqrt(len(del_vols)))
#         pow_sa_se.append(np.std(pow_sas)/np.sqrt(len(pow_sas)))
#         del_sa_se.append(np.std(del_sas)/np.sqrt(len(del_sas)))

#     # Create the dictionary for converting the labels
#     graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet', 'pl_complex': 'Prot-Lig',
#                      'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'BSA': 'BSA', '1BNA': '1BNA',
#                      'DB1976': 'DB1976'}[_] for _ in my_sys_names]
#     # Set the label codes
#     code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
#                  'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J'}
#     new_graph_labels = [code_dict[_] for _ in graph_labels]

#     def sort_3_lists(lista, listb):
#         # Zipping lists together and sorting by the first list
#         sorted_lists = sorted(zip(lista, listb), key=lambda x: x[0])

#         # Unpacking the sorted lists
#         lista, listb = zip(*sorted_lists)

#         # Converting tuples back to lists if needed
#         lista = list(lista)
#         listb = list(listb)

#         # Return the lists
#         return lista, listb

#     new_graph_labels, pre_made_data = sort_3_lists(new_graph_labels, new_graph_labels)

#     # Create the bar graph
#     bar([[_[1] for _ in pre_made_data], [_[3] for _ in pre_made_data]], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
#         Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Volume Difference',
#         errors=[[_[2] for _ in pre_made_data], [_[4] for _ in pre_made_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)
#     # Create the bar graph
#     bar([[_[5] for _ in pre_made_data], [_[7] for _ in pre_made_data]], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
#         Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Surface Area Difference',
#         errors=[[_[6] for _ in pre_made_data], [_[8] for _ in pre_made_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)

