import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
from System.system import System
from System.Group.group import Group
from Data.Analyze.tools.plot_templates.bar import bar
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.get_res_data import residue_data

# pre_made_data = [
#     ['cambrin', 2.657358947567952, 0.1701522540304929, 2.59465391164789, 0.1701522540304929, 0.7071402007950848, 0.059791991731881135, 2.3096433988252345, 0.21297727309940268],
#     ['hairpin', 2.120015342765748, 0.17989523306316388, 4.207094792681904, 0.17989523306316388, 0.3151423480069726, 0.050484980940870186, 0.587710293326916, 0.09756079607665612],
#     ['p53tet', 2.9331662712893984, 0.18342810870548817, 2.7747918309068043, 0.18342810870548817, 0.6593449473781604, 0.10422962591710631, 2.370253686067341, 0.25677219904214393],
#     ['pl_complex', 2.7747879360584684, 0.07901575226030394, 2.8010158048936318, 0.07901575226030394, 0.62219526460963, 0.03121395618438131, 2.422886871774085, 0.11430094978326413],
#     ['streptavidin', 3.0522725810134834, 0.1238965156383565, 2.174625260864225, 0.1238965156383565, 0.4774545212283423, 0.0397510499787414, 2.072412116883273, 0.1310601195844641],
#     ['hammerhead', 6.576136641455463, 1.3355495236483024, 8.36539089536273, 1.3355495236483024, 7.027173386404598, 1.337333613171973, 5.30256130664321, 2.34512218840759],
#     ['NCP', 8.826521657702823, 0.9350257210684129, 8.618442644565837, 0.9350257210684129, 1.7301587766519828, 0.19436185835523279, 7.435837930712498, 1.0235150388051224],
#     ['BSA', 2.8614722131257735, 0.5170327576040869, 8.53612886017547, 1.5170327576040869, 0.6362188314154772, 0.2564681508906418, 1.47362166167768, 0.16446425197757493]
# ]


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
        print(system.name)
        vor_out, vor_in = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_vor_res.csv'.format(system.name)):
            vor_in, vor_out = folder + '/res_data/{}_vor_res.csv'.format(system.name), None
        pow_out, pow_in = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_pow_res.csv'.format(system.name)):
            pow_in, pow_out = folder + '/res_data/{}_pow_res.csv'.format(system.name), None
        del_out, del_in = folder + '/res_data/{}_del_res.csv'.format(system.name), None
        if os.path.exists(folder + '/res_data/{}_del_res.csv'.format(system.name)):
            del_in, del_out = folder + '/res_data/{}_del_res.csv'.format(system.name), None
        # print('45: system = {}'.format(system.name))
        pow_vols, del_vols, pow_sas, del_sas = [], [], [], []
        # Get the values from the residue function
        vor_reses = residue_data(system, read_logs(my_log_files[system.name]['vor']), get_all=True, read_file=vor_in, output_file=vor_out)
        # print('49: vor_reses = {}'.format(vor_reses))
        pow_reses = residue_data(system, read_logs(my_log_files[system.name]['pow']), get_all=True, read_file=pow_in, output_file=pow_out)
        # print('51: pow_reses = {}'.format(pow_reses))
        del_reses = residue_data(system, read_logs(my_log_files[system.name]['del']), get_all=True, read_file=del_in, output_file=del_out)
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
                    pow_vol_diff = (pow_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
                    del_vol_diff = (pow_reses[_][__][___]['vol'] - vor_reses[_][__][___]['vol']) / vor_reses[_][__][___]['vol']
                    pow_vols.append(pow_vol_diff)
                    del_vols.append(del_vol_diff)
                    pow_sa_diff = (vor_reses[_][__][___]['sa'] - pow_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa']
                    del_sa_diff = (vor_reses[_][__][___]['sa'] - del_reses[_][__][___]['sa']) / vor_reses[_][__][___]['sa']
                    pow_sas.append(pow_sa_diff)
                    del_sas.append(del_sa_diff)
                    if abs(del_vol_diff) > 100:
                        print(__, ___, del_vol_diff, del_sa_diff)
        # Get the averages
        print('AvgResDiffs 73: ', system.name)
        print('AvgResDiffs 74: ', 100 * sum(pow_vols)/len(pow_vols), 100 * np.std(pow_vols)/np.sqrt(len(pow_vols)))
        print('AvgResDiffs 75: ', 100 * sum(del_vols)/len(del_vols), 100 * np.std(del_vols)/np.sqrt(len(del_vols)))
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
    graph_labels = [{'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet', 'pl_complex': 'Prot-Lig',
                     'streptavidin': 'STVDN', 'hammerhead': 'H-Head', 'NCP': 'NCP', 'BSA': 'BSA', '1BNA': '1BNA',
                     'DB1976': 'DB1976'}[_[0]] for _ in pre_made_data]
    # Set the label codes
    code_dict = {'Na5': 'A', 'EDTA': 'B', 'Hairpin': 'C', 'Cambrin': 'D', 'H-Head': 'E', 'p53tet': 'F',
                 'Prot-Lig': 'G', 'STVDN': 'H', 'NCP': 'I', 'BSA': 'J'}
    new_graph_labels = [code_dict[_] for _ in graph_labels]

    def sort_3_lists(lista, listb):
        # Zipping lists together and sorting by the first list
        sorted_lists = sorted(zip(lista, listb), key=lambda x: x[0])

        # Unpacking the sorted lists
        lista, listb = zip(*sorted_lists)

        # Converting tuples back to lists if needed
        lista = list(lista)
        listb = list(listb)

        # Return the lists
        return lista, listb

    new_graph_labels, pre_made_data = sort_3_lists(new_graph_labels, pre_made_data)

    # Create the bar graph
    bar([[_[1] for _ in pre_made_data], [_[3] for _ in pre_made_data]], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Volume Difference',
        errors=[[_[2] for _ in pre_made_data], [_[4] for _ in pre_made_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)
    # Create the bar graph
    bar([[_[5] for _ in pre_made_data], [_[7] for _ in pre_made_data]], x_names=new_graph_labels, legend_names=['Power', 'Primitive'],
        Show=True, y_axis_title='% Difference', x_axis_title='Model', title='Average Residue Surface Area Difference',
        errors=[[_[6] for _ in pre_made_data], [_[8] for _ in pre_made_data]], y_range=[0, None], xtick_label_size=25, ytick_label_size=25, ylabel_size=30, xlabel_size=30, tick_length=12, tick_width=2)

