import csv
import os
import time

import numpy as np
import tkinter as tk
from tkinter import filedialog
from vorpy.src.system.system import System
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.sys_funcs.calcs.calcs import calc_sphericity, calc_dist, get_time


def get_logs_and_pdbs(make_file=True, output_file_name=None, cv=None, density=None):
    logs_pdbs = {}

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    user_data = filedialog.askdirectory(title='Get User Data')
    for rroot, directories, files in os.walk(user_data):
        for directory in directories:
            if 'aw' in directory or 'pow' in directory:
                continue
            logs_pdbs[directory] = {}
            for rrooot, dircs, filese in os.walk(rroot + '/' + directory):
                for file in filese:
                    if file[-3:] == 'pdb' and 'atoms' not in file:
                        logs_pdbs[directory]['pdb'] = rrooot + '/' + file
                # print(filese)
                for dircy_dirc in dircs:
                    if dircy_dirc[-2:] == 'aw':
                        for rootytooty, dincretories, flies in os.walk(rrooot):
                            for file in flies:
                                if file[-3:] == 'csv' and rootytooty[-2:] == 'aw':
                                    logs_pdbs[directory]['aw'] = rootytooty + '/' + file
                                    # print(rootytooty + '/' + file)
                    elif dircy_dirc[-3:] == 'pow':
                        for rootytooty, dincretories, flies in os.walk(rrooot):
                            for file in flies:
                                if file[-3:] == 'csv' and rootytooty[-3:] == 'pow':
                                    logs_pdbs[directory]['pow'] = rootytooty + '/' + file
                # print(logs_pdbs[directory])
    if make_file:
        if output_file_name is None:
            output_file_name = 'logs_pdbs.txt'
        with open(output_file_name, 'w') as loggy_woggys:
            for _ in logs_pdbs:
                # print(logs_pdbs[_])
                if 'pdb' in logs_pdbs[_] and 'aw' in logs_pdbs[_] and 'pow' in logs_pdbs[_]:
                    loggy_woggys.write(logs_pdbs[_]['pdb'] + '\n')
                    loggy_woggys.write(logs_pdbs[_]['aw'] + '\n')
                    loggy_woggys.write(logs_pdbs[_]['pow'] + '\n')

    return logs_pdbs


def make_new_logs(logs_pdbs_dict=None, logs_pdbs_file=None):

    if logs_pdbs_file is not None:
        logs_pdbs_dict = {}
        with open(logs_pdbs_file, 'r') as loggy_woggys:
            key = None
            for i, line in enumerate(loggy_woggys.readlines()):

                if i % 3 == 0:
                    split_line = line.split('/')
                    key = split_line[-2]
                    logs_pdbs_dict[key] = {'pdb': line[:-1]}
                elif i % 3 == 1:
                    logs_pdbs_dict[key]['aw'] = line[:-1]
                elif i % 3 == 2:
                    logs_pdbs_dict[key]['pow'] = line[:-1]

    new_info_dict = {}
    length_of_dic = len(logs_pdbs_dict)
    counter = 0
    start = time.perf_counter()
    with open('Isoq.csv', 'w') as new_logs:
        n_logs = csv.writer(new_logs)
        n_logs.writerow(['file', "avg isoq aw", 'avg isoq pow', 'avg abs diff isoq (aw base)', 'avg abs diff isoq (aw base)'])
        for _ in logs_pdbs_dict:
            counter += 1
            my_sys = System(logs_pdbs_dict[_]['pdb'], simple=True)
            aw_logs = read_logs2(logs_pdbs_dict[_]['aw'])
            pow_logs = read_logs2(logs_pdbs_dict[_]['pow'])
            new_info_dict[_] = {'isoq_aw': [], 'isoq_pow': [], 'rads': []}

            for i, ball in my_sys.balls.iterrows():
                # Get the aw_ball and the pow ball
                try:
                    aw_ball = aw_logs['atoms'][aw_logs['atoms']['Index'] == ball['num']].iloc[0].to_dict()
                    pow_ball = pow_logs['atoms'][pow_logs['atoms']['Index'] == ball['num']].iloc[0].to_dict()
                except ValueError:
                    continue
                except IndexError:
                    continue
                # if len(aw_ball['num']) == 0:
                #     continue
                # aw_avg_sa = aw_ball['sa'] / aw_num_neighbors
                # print(i, aw_ball)
                # Get neighbors
                new_info_dict[_]['isoq_aw'].append(aw_ball['Isometric Quotient'])
                new_info_dict[_]['isoq_pow'].append(pow_ball['Isometric Quotient'])

                # Add the radius of the ball to the information dictionary
                new_info_dict[_]['rads'].append(ball['rad'])

            # Get the dictionary
            dicty = new_info_dict[_]
            # Get the
            diffs1 = [abs(dicty['isoq_aw'][i] - dicty['isoq_pow'][i]) / dicty['isoq_aw'][i]
                      for i in range(len(dicty['isoq_aw']))]
            diffs2 = [abs(dicty['isoq_aw'][i] - dicty['isoq_pow'][i]) / dicty['isoq_pow'][i]
                      for i in range(len(dicty['isoq_aw']))]
            try:
                line = [_, np.mean(dicty['isoq_aw']), np.mean(dicty['isoq_pow']), np.mean(diffs1), np.mean(diffs2)]
            except:
                continue
            timey_wimey = get_time(time.perf_counter() - start)
            print('{}/{} Done, {} %, Time elapsed = {}:{}:{}  -   Data ---->>>     '
                  .format(counter, length_of_dic, round(100 * counter / length_of_dic, 3), round(timey_wimey[0]),
                          round(timey_wimey[1]), round(timey_wimey[2], 2)), line)
            n_logs.writerow(line)

    return new_info_dict


if __name__ == '__main__':
    get_logs_and_pdbs(True, output_file_name='loggy_woggys1.txt')
    os.chdir('../../../..')
    noopy_ploopster = make_new_logs(logs_pdbs_file=os.getcwd() + '/Data/Analyze/tools/batch/loggy_woggys1.txt')
