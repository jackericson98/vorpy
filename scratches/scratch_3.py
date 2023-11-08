import csv

import numpy as np
import plotly.graph_objects as go

with open('C:/Users/jacke/PycharmProjects/vorpy/Data/user_data/foam_data.csv', 'r') as foam_data:
    f_dat = csv.reader(foam_data)
    my_data = []
    x, y, z1, z2 = [], [], [], []
    for line in foam_data.readlines()[1:]:
        my_line = line.split(',')
        if len(my_line) == 1:
            continue
        # print(foam_params)
        print(my_line)
        my_data.append({'num': my_line[0], 'avg rad size': float(my_line[2]), 'box size': float(my_line[1]),
                        'rad std': float(my_line[3]), 'num balls': int(my_line[4]),
                        'density': float(my_line[5]), 'vol diff vor': my_line[6], 'sa diff vor': my_line[7],
                        'vol diff pow': my_line[8], 'sa diff pow': my_line[9], 'num cells': my_line[10]})

    lists = {}
    nums = []
    under_100 = []
    i = 0
    while i < len(my_data):
        nums.append(int(my_data[i]['num']))
        # if int(my_data[i]['num']) != i + 1:
        #     print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/foam_{}/lognormal.pdb -s nt compare -s mv 1000'.format(i + 1))

        if int(my_data[i]['num cells']) < 100:
            under_100.append(int(my_data[i]['num']))
            # print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/foam_{}/lognormal.pdb -s nt compare -s mv 1000'.format(i + 1))
        i += 1
    for i in range(1, 1301):
        if i not in nums:
            under_100.append(i)
    under_100.sort()
    for _ in under_100:
        print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/foam_{}/lognormal.pdb -s nt compare -s mv 1000'.format(_))

        # dp = my_data[i]
        # # Check if the data has been added before
        # if dp['rad std'] in lists:
        #     if dp['density'] in lists[dp['rad std']]:
        #         lists[dp['rad std']][dp['density']][0].append(dp['vol diff vor'])
        #         lists[dp['rad std']][dp['density']][1].append(dp['sa diff vor'])
        #         lists[dp['rad std']][dp['density']][2].append(dp['vol diff pow'])
        #         lists[dp['rad std']][dp['density']][3].append(dp['sa diff pow'])
        #     else:
        #         lists[dp['rad std']][dp['density']] = [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']], [dp['sa diff pow']]]
        # else:
        #     lists[dp['rad std']] = {dp['density']: [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']], [dp['sa diff pow']]]}
    # xi = np.linspace(0.05, 0.5, 10)
    # yi = np.linspace(0.05, 0.5, 10)
    # my_grid = np.meshgrid(xi, yi)
    
    # fig = go.Figure(data=[go.Surface(x=x, y=y, z=z1)])
    # fig.show()
        