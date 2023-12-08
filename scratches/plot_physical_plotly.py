import csv
import pandas as pd
import numpy as np
import plotly.graph_objects as go


plot_type = 'physical3'

with open('C:/Users/jacke/PycharmProjects/vorpy/Data/user_data/foam_data.csv', 'r') as my_foam_data:
    my_data = []
    # x, y, z1, z2 = [], [], [], []
    foam_data = csv.reader(my_foam_data)
    for my_line in foam_data:
        if len(my_line) <= 1 or plot_type not in my_line[0] or int(my_line[4]) < 100:
            continue
        my_data.append({'num': my_line[0], 'avg rad size': float(my_line[2]), 'box size': float(my_line[1]),
                        'rad std': float(my_line[3]), 'num balls': int(my_line[4]),
                        'density': float(my_line[5]), 'vol diff vor': 100 * float(my_line[6]),
                        'sa diff vor': 100 * float(my_line[7]),
                        'vol diff pow': 100 * float(my_line[8]), 'sa diff pow': 100 * float(my_line[9]),
                        'num cells': int(my_line[10])})


lists = {}
for dp in my_data:
    # Check if the data has been added before
    if dp['density'] in lists:
        lists[dp['density']][0].append(dp['vol diff vor'])
        lists[dp['density']][1].append(dp['sa diff vor'])
        # lists[dp['density']][2].append(dp['vol diff pow'])
        # lists[dp['density']][3].append(dp['sa diff pow'])
    else:
        # lists[dp['density']] = [[dp['vol diff vor']], [dp['sa diff vor']], [dp['vol diff pow']],
        #                                        [dp['sa diff pow']]]
        lists[dp['density']] = [[dp['vol diff vor']], [dp['sa diff vor']]]
# for _ in lists:
#     print(_, [len(lists[_][i]) for i in range(4)])

# Calculate averages across sets at each index
averages = [
    [sum(set_values[i]) / len(set_values[i]) for set_values in lists.values()]
    for i in range(len(lists[0.025]))
]
std_errs = [[np.std(set_values[i]) / np.sqrt(len(set_values[i])) for set_values in lists.values()]
            for i in range(len(lists[0.025]))]


# Extract x values for plotting
x_values = list(lists.keys())

# Create traces for each line
trace_names = ['Volume', 'Surface Area']
traces = [go.Scatter(x=x_values, y=avg, mode='lines', name=trace_names[i], error_y=dict(type='data', array=std_errs[i], visible=True)) for i, avg in enumerate(averages)]

# Create layout for the plot
layout = go.Layout(title=dict(text='Gal-Or & Hoelsher AWVd vs. Power for Densities', font=dict(size=40)),
                   xaxis=dict(title='Density', tickfont=dict(size=16), titlefont=dict(size=30)),
                   yaxis=dict(title='Percent Difference', tickfont=dict(size=16), titlefont=dict(size=30)),
                   legend=dict(font=dict(size=30)))

# Create figure
fig = go.Figure(data=traces, layout=layout)

# Show the plot
fig.show()