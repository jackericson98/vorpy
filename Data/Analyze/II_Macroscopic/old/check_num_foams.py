import os
from os import path
import tkinter as tk
from tkinter import filedialog

# root = tk.Tk()
# root.withdraw()
# root.wm_attributes('-topmost', 1)
# rootdir = filedialog.askdirectory()

rootdir = '/home/jack/PycharmProjects/foam_gen/Data/user_data'

my_file_types = {}
for i in range(20):
    for j in range(17):
        i_num = str(float(round((i + 1) * 0.025, 3)))
        j_num = str(float(round((j + 4) * 0.5, 3)))
        if j_num in my_file_types:
            if i_num not in my_file_types[j_num]:
                my_file_types[j_num][i_num] = 0
        else:
            my_file_types[j_num] = {i_num: 0}

for dir_ in os.listdir(rootdir):
    dir_info = dir_.split('_')
    try:
        my_file_types[dir_info[1]][dir_info[3]] += 1
    except KeyError:
        continue
count = 0
reverse_list, my_list = [], []
for _ in my_file_types:
    for __ in my_file_types[_]:
        if my_file_types[_][__] < 20:
            print(_, __, my_file_types[_][__])
        for i in range(20 - my_file_types[_][__]):
            count += 1
            if my_file_types[_][__] == 0 and i == 0:
                my_list.insert(0, 'python3 foam_gen.py 10.0 {} 1000 {} False gamma\n'.format(_, __))
            else:
                my_list.append('python3 foam_gen.py 10.0 {} 1000 {} False gamma\n'.format(_, __))


        # print('python3 foam_gen.py', *poopy)
        # reverse_list.insert(0, ['python3 foam_gen.py', *poopy])
# for _ in reverse_list:
#     print(*_)
# my_strs = [x for x, _ in sorted(my_list, key=lambda x: x[1])]
#
# for srt in my_strs:
#     print(*srt)

# for i, _ in enumerate(my_list):
#     i_div = i // 200
    # if not path.exists('foamify_{}.sh'.format(i_div)):
        # with open('foamify_{}.sh'.format(i_div), 'a') as gay:
        #     gay.write('#!/bin/sh\n')
    # print(_)
    # with open('foamify_{}.sh'.format(i_div), 'a') as gay:
    #     # gay.write(_ + '\n')

# num_files = 12
# for i, _ in enumerate(my_list):
#     run_file_num = i % 12
#     if not path.exists('foam_runs_' + str(run_file_num) + '.sh'):
#         with open('foam_runs_' + str(run_file_num) + '.sh', 'a') as write_file:
#             write_file.write('#!/bin/bash\n')
#
#     with open('foam_runs_' + str(run_file_num) + '.sh', 'a') as write_file:
#         write_file.write(_)


print("{} Runs Left".format(count))


# cvs = {}
# for _ in my_file_types:
#     name = _.split('_')
#     cv = name[1]
#     dens = name[3]
#     if cv in cvs:
#         if dens in cvs[cv]:
#             cvs[cv][dens] += my_file_types[_]
#         else:
#             cvs[cv][dens] = my_file_types[_]
#     else:
#         cvs[cv] = {dens: my_file_types[_]}
#
# print()
