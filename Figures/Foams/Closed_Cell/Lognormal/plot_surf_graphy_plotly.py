import csv
import os

import pandas as pd
import numpy as np
import plotly.graph_objects as go


plot_type = 'lognormal'
for i in range(4):
    os.chdir('..')
with open(os.getcwd() + '/Data/user_data/foam_data.csv', 'r') as my_foam_data:

    my_data = []
    # x, y, z1, z2 = [], [], [], []
    foam_data = csv.reader(my_foam_data)
    for my_line in foam_data:
        if len(my_line) <= 1 or plot_type not in my_line[0] or int(my_line[4]) < 100:
            continue
        my_data.append({'num': my_line[0], 'avg rad size': float(my_line[2]), 'box size': float(my_line[1]),
                        'rad std': float(my_line[3]), 'num balls': int(my_line[4]),
                        'density': float(my_line[5]), 'vol diff vor': float(my_line[6]),
                        'sa diff vor': float(my_line[7]),
                        'vol diff pow': float(my_line[8]), 'sa diff pow': float(my_line[9]),
                        'num cells': int(my_line[10])})


lists = {}

for dp in my_data:
    # Check if the data has been added before
    if dp['rad std'] in lists:
        if dp['density'] in lists[dp['rad std']]:
            lists[dp['rad std']][dp['density']][0].append(dp['vol diff vor'])
            lists[dp['rad std']][dp['density']][1].append(dp['sa diff vor'])
            lists[dp['rad std']][dp['density']][2].append(dp['vol diff pow'])
            lists[dp['rad std']][dp['density']][3].append(dp['sa diff pow'])
        else:
            lists[dp['rad std']][dp['density']] = [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']],
                                                   [dp['sa diff pow']]]
    else:
        lists[dp['rad std']] = {
            dp['density']: [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']], [dp['sa diff pow']]]}


if plot_type == 'gamma':
    my_densities = np.arange(0.025, 0.45, 0.45)
    my_sds = np.arange(1.5, 6.5, 0.25)

if plot_type == 'lognormal':
    my_densities = np.arange(0.025, 0.475, 0.025)
    my_sds = np.arange(0.1, 0.525, 0.025)

my_densities = [round(_, 3) for _ in my_densities]
my_sds = [round(_, 3) for _ in my_sds]

# Initialize lists using list comprehensions
datavvm, datavsm, datapvm, datapsm = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvms, datavsms, datapvms, datapsms = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvps, datavsps, datapvps, datapsps = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]

# Iterate over my_sds and my_densities using nested loops
for i, num in enumerate(my_sds):
    for j, num2 in enumerate(my_densities):
        means = [np.mean(lst) * 100 for lst in lists[num][num2]]
        sds = [np.std([_ * 100 for _ in lst]) / np.sqrt(len(lst)) for lst in lists[num][num2]]

        # Append means and mean +/- SD to respective lists
        datavvm[i].append(means[0]); datavvms[i].append(means[0] - sds[0]); datavvps[i].append(means[0] + sds[0])
        datavsm[i].append(means[1]); datavsms[i].append(means[1] - sds[1]); datavsps[i].append(means[1] + sds[1])
        datapvm[i].append(means[2]); datapvms[i].append(means[2] - sds[2]); datapvps[i].append(means[2] + sds[2])
        datapsm[i].append(means[3]); datapsms[i].append(means[3] - sds[3]); datapsps[i].append(means[3] + sds[3])


xi = my_sds
yi = my_densities


# Convert the matrix to a NumPy array
matrix_array = np.array(datavvm)

# Create x and y coordinates
x_coords, y_coords = np.meshgrid(my_sds, my_densities)

# Create a surface plot
fig = go.Figure(data=[go.Surface(z=matrix_array, x=x_coords, y=y_coords)])

# Update layout
fig.update_layout(scene=dict(
                    xaxis_title='Coefficient of Variation',
                    yaxis_title='Density',
                    zaxis_title='% difference'))

# Show the plot
fig.show()
#
# df_vvm = pd.DataFrame(datavvm, xi, xi)
# df_vsm = pd.DataFrame(datavsm, xi, xi)
#
# fig = go.Figure(
#     data=[go.Surface(x=xi, y=yi, z=datavvm), go.Surface(x=xi, y=yi, z=datavvms, showscale=False, opacity=0.5),
#           go.Surface(x=xi, y=yi, z=datavvps, showscale=False, opacity=0.5)])
# fig1 = go.Figure(
#     data=[go.Surface(x=xi, y=yi, z=datavsm), go.Surface(x=xi, y=yi, z=datavsms, showscale=False, opacity=0.5),
#           go.Surface(x=xi, y=yi, z=datavsps, showscale=False, opacity=0.5)])
# fig2 = go.Figure(data=[go.Surface(x=xi, y=yi, z=datapvm), go.Surface(x=xi, y=yi, z=datapvms),
#                        go.Surface(x=xi, y=yi, z=datapvps)])
# fig3 = go.Figure(data=[go.Surface(x=xi, y=yi, z=datapsm), go.Surface(x=xi, y=yi, z=datapsms, opacity=0.5),
#                        go.Surface(x=xi, y=yi, z=datapsps, opacity=0.5)])
# # fig3.add_table(datapsm)
# fig.update_layout(title="Average Percent Difference - Volume")
# fig.update_scenes(xaxis_title_text="Density",
#                   yaxis_title_text="Coefficient of Variation", zaxis_title_text="% Difference")
# fig.add_table(cells=dict(values=df_vvm.values.tolist()))
#
# fig1.update_layout(title='Average Percent Difference - Surface Area')
# fig1.update_scenes(xaxis_title_text="Density",
#                    yaxis_title_text="Coefficient of Variation", zaxis_title_text="% Difference")
# fig.show()
# fig1.show()
# # fig2.show()
# # fig3.show()
