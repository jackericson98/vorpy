import os
import tkinter as tk
from tkinter import filedialog

import matplotlib.pyplot as plt

from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.plot_templates.line import line_plot
from vorpy.src.system.system import System


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
folder = filedialog.askdirectory()

aw_logs, pow_logs, prm_logs, pdb_paths = [], [], [], []
for file in os.listdir(folder):
    if file[-4:] == '.pdb':
        pdb_paths.append(folder + '/' + file)
    elif file[-6:] == 'aw.csv':
        aw_logs.append(folder + '/' + file)
    elif file[-7:] == 'prm.csv':
        prm_logs.append(folder + '/' + file)
    elif file[-7:] == 'pow.csv':
        pow_logs.append(folder + '/' + file)
    else:
        print('Bad File Read: ', file)

vols = [[], [], []]
sas = [[], [], []]

vol_devs = [[], []]
sa_devs = [[], []]

for i in range(len(pdb_paths)):
    my_system = System(file=pdb_paths[i])
    aw_info = read_logs(aw_logs[i])
    prm_info = read_logs(prm_logs[i])
    pow_info = read_logs(pow_logs[i])
    aw_vol, pow_vol, prm_vol = aw_info['group data']['volume'], pow_info['group data']['volume'], prm_info['group data']['volume']
    aw_sa, pow_sa, prm_sa = aw_info['group data']['sa'], pow_info['group data']['sa'], prm_info['group data']['sa']

    vols[0].append(aw_vol)
    vols[1].append(pow_vol)
    vols[2].append(prm_vol)
    sas[0].append(aw_sa)
    sas[1].append(pow_sa)
    sas[2].append(prm_sa)
    vol_devs[0].append(100 * (pow_vol - aw_vol) / aw_vol)
    vol_devs[1].append(100 * (prm_vol - aw_vol) / aw_vol)
    sa_devs[0].append(100 * (pow_sa - aw_sa) / aw_sa)
    sa_devs[1].append(100 * (prm_sa - aw_sa) / aw_sa)


frm_vol_devs = [[100 * (vols[j][i] - vols[j][0]) / vols[j][0] for i in range(len(vols[0]) - 1)] for j in range(3)]
frm_sa_devs = [[100 * (sas[j][i] - sas[j][0]) / sas[j][0] for i in range(len(sas[0]) - 1)] for j in range(3)]

x_data = [_ + 1 for _ in range(len(pdb_paths) - 1)]

x_data1 = [_ + 1 for _ in range(len(pdb_paths))]

# Plot 1: Frame volume deviation
line_plot([x_data for _ in range(3)], frm_vol_devs, Show=False, colors=['k', 'skyblue', 'orange'],
          title='Cambrin Frame Volume Deviation', y_label='% Difference', x_label='Frame', labels=['Power', 'Primitive', 'Additively Weighted'], legend_orientation='Horizontal', x_label_size=35, y_label_size=35, tick_val_size=30, tick_width=2, tick_length=12, linewidth=4)
# plot 2: Frame SA deviation
line_plot([x_data for _ in range(3)], frm_sa_devs, Show=False, colors=['k', 'skyblue', 'orange'],
          title='Cambrin Frame SA Deviation', y_label='% Difference', x_label='Frame', labels=['Power', 'Primitive', 'Additively Weighted'], x_label_size=35, y_label_size=35, tick_val_size=30, tick_width=2, tick_length=12, linewidth=4)
# Plot 3: Deviation from AW volume
line_plot([x_data1 for _ in range(2)], vol_devs, Show=False, title='Cambrin Volume Variation from AWVd', x_label='Frame', colors=['orange', 'green'],
          y_label='% Difference', x_label_size=35, y_label_size=35, tick_val_size=30, tick_width=2, tick_length=12, linewidth=4)
# Plot 4: Deviation from AW SA
line_plot([x_data1 for _ in range(2)], sa_devs, title='Cambrin SA Variation from AWVd', x_label='Frame', colors=['orange', 'green'],
          y_label='% Difference', legend_orientation='Horizontal', x_label_size=35, y_label_size=35, tick_val_size=30, tick_width=2, tick_length=12, linewidth=4)
# plt.plot(x_data, data)
# plt.show()


