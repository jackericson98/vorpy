import csv
import pandas as pd
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
        my_data.append({'num': my_line[0], 'avg rad size': float(my_line[2]), 'box size': float(my_line[1]),
                        'rad std': float(my_line[3]), 'num balls': int(my_line[4]),
                        'density': float(my_line[5]), 'vol diff vor': float(my_line[6]),
                        'sa diff vor': float(my_line[7]),
                        'vol diff pow': float(my_line[8]), 'sa diff pow': float(my_line[9]),
                        'num cells': int(my_line[10])})

lists = {}
for dp in my_data:
    if float(dp['vol diff vor']) > 1:
        continue
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

my_nums = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
datavvm, datavsm, datapvm, datapsm = [], [], [], []
datavvms, datavsms, datapvms, datapsms = [], [], [], []
datavvps, datavsps, datapvps, datapsps = [], [], [], []
for num in my_nums:
    datavvm.append([])
    datavsm.append([])
    datapvm.append([])
    datapsm.append([])
    datavvms.append([])
    datavsms.append([])
    datapvms.append([])
    datapsms.append([])
    datavvps.append([])
    datavsps.append([])
    datapvps.append([])
    datapsps.append([])
    for num2 in my_nums:
        means, sds = [], []
        for i in range(4):
            means.append(np.mean(lists[num][num2][i]) * 100)
            sds.append(np.std([_ * 100 for _ in lists[num][num2][i]]) / np.sqrt(len(lists[num][num2])))
        datavvm[-1].append(means[0])
        datavsm[-1].append(means[1])
        datapvm[-1].append(means[2])
        datapsm[-1].append(means[3])
        datavvps[-1].append(means[0] + sds[0])
        datavsps[-1].append(means[1] + sds[1])
        datapvps[-1].append(means[2] + sds[2])
        datapsps[-1].append(means[3] + sds[3])
        datavvms[-1].append(means[0] - sds[0])
        datavsms[-1].append(means[1] - sds[1])
        datapvms[-1].append(means[2] - sds[2])
        datapsms[-1].append(means[3] - sds[3])

xi = np.linspace(0.05, 0.5, 10)
yi = np.linspace(0.05, 0.5, 10)
df_vvm = pd.DataFrame(datavvm, xi, xi)
df_vsm = pd.DataFrame(datavsm, xi, xi)

fig = go.Figure(
    data=[go.Surface(x=xi, y=yi, z=datavvm), go.Surface(x=xi, y=yi, z=datavvms, showscale=False, opacity=0.5),
          go.Surface(x=xi, y=yi, z=datavvps, showscale=False, opacity=0.5)])
fig1 = go.Figure(
    data=[go.Surface(x=xi, y=yi, z=datavsm), go.Surface(x=xi, y=yi, z=datavsms, showscale=False, opacity=0.5),
          go.Surface(x=xi, y=yi, z=datavsps, showscale=False, opacity=0.5)])
fig2 = go.Figure(data=[go.Surface(x=xi, y=yi, z=datapvm), go.Surface(x=xi, y=yi, z=datapvms),
                       go.Surface(x=xi, y=yi, z=datapvps)])
fig3 = go.Figure(data=[go.Surface(x=xi, y=yi, z=datapsm), go.Surface(x=xi, y=yi, z=datapsms, opacity=0.5),
                       go.Surface(x=xi, y=yi, z=datapsps, opacity=0.5)])
# fig3.add_table(datapsm)
fig.update_layout(title="Average Percent Difference - Volume")
fig.update_scenes(xaxis_title_text="Density",
                  yaxis_title_text="Coefficient of Variation", zaxis_title_text="% Difference")
fig.add_table(cells=dict(values=df_vvm.values.tolist()))

fig1.update_layout(title='Average Percent Difference - Surface Area')
fig1.update_scenes(xaxis_title_text="Density",
                   yaxis_title_text="Coefficient of Variation", zaxis_title_text="% Difference")
fig.show()
fig1.show()
# fig2.show()
# fig3.show()
