# import os
#
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.optimize import curve_fit
# from Data.Analyze.tools.plot_templates.line import line_plot
#
#
# xs, ys, lognormals, num_balls = [], [], [], []
# data = {}
# with open(os.getcwd() + '/density_adjustments.txt', 'r') as density_adjustments:
#     for line in density_adjustments.readlines():
#         split_line = line.split()
#         num_balls = split_line[3]
#         # if split_line[5] == 'lognormal':
#         #     continue
#         if num_balls in data:
#
#             if float(split_line[4]) in data[num_balls]['data']:
#                 data[num_balls]['data'][float(split_line[4])].append(float(split_line[0]))
#             else:
#
#                 data[num_balls]['data'][float(split_line[4])] = [float(split_line[0])]
#             data[num_balls]['cv'].append(float(split_line[2]))
#         else:
#             data[num_balls] = {'data': {float(split_line[4]): [float(split_line[0])]}, 'cv': [float(split_line[2])]}
#
#
# # Define the model function
# def sqrt_model(x, a, b, c):
#     return a * x ** 2 + b * x + c
#
#
# colors = ['skyblue', 'orange', 'red', 'blue', 'green', 'pink']
# sorted_data = {key: data[key] for key in sorted(data)}
# cmap = plt.cm.get_cmap('rainbow')
# x_data, y_data, err_data = [], [], []
# for i, _ in enumerate(sorted_data):
#     xs = [__ for __ in sorted_data[_]['data']]
#     ys = [np.mean(sorted_data[_]['data'][__]) for __ in sorted_data[_]['data']]
#     std_errs = [np.std(sorted_data[_]['data'][__]) / np.sqrt(len(sorted_data[_]['data'][__])) for __ in sorted_data[_]['data']]
#     # Perform the curve fitting
#     params, cov = curve_fit(sqrt_model, xs, ys)
#
#     # Extract the parameters
#     a, b, c = params
#     x_data.append(xs)
#     y_data.append([a * x ** 2 + b * x + c for x in xs])
#     err_data.append(std_errs)
#
# line_plot(x_data, y_data, errors=err_data, labels=[_ for _ in sorted_data])

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os


xs, ys, lognormals, num_balls = [], [], [], []
data = {}
with open(os.getcwd() + '/density_adjustments.txt', 'r') as density_adjustments:
    for line in density_adjustments.readlines():
        split_line = line.split()
        num_balls = split_line[3]
        # if split_line[5] == 'lognormal':
        #     continue
        if num_balls in data:
            data[num_balls]['xs'].append(float(split_line[4]))
            data[num_balls]['ys'].append(float(split_line[0]))
            data[num_balls]['cv'].append(float(split_line[2]))
        else:
            data[num_balls] = {'xs': [float(split_line[4])], 'ys': [float(split_line[0])], 'cv': [float(split_line[2])]}


# Define the model function
def sqrt_model(x, a, b, c):
    return a * x ** 2 + b * x + c


colors = ['skyblue', 'orange', 'red', 'blue', 'green', 'pink']

sorted_data = {key: data[key] for key in sorted(data)}
cmap = plt.cm.get_cmap('rainbow')

for i, _ in enumerate(sorted_data):
    xs, ys = sorted_data[_]['xs'], sorted_data[_]['ys']
    colors1 = [cmap(_) for _ in sorted_data[_]['cv']]
    # Perform the curve fitting
    params, cov = curve_fit(sqrt_model, xs, ys)

    # Extract the parameters
    a, b, c = params
    print('{} Balls - y = {}x^2 + {}x + {}'.format(_, a, b, c))
    x_fit = np.linspace(min(xs), max(xs), 100)
    y_fit = [a * x ** 2 + b * x + c for x in x_fit]

    plt.plot(x_fit, y_fit, c=colors[i], label=_)
    plt.scatter(xs, ys, c=colors[i], alpha=0.1)

plt.ylabel("Closed Density", fontdict=dict(size=25))
plt.xlabel('Open Density', fontdict=dict(size=25))
# plt.xticks(rotation=45, ha='right', font=dict(size=xtick_label_size))
plt.yticks([0.1, 0.3, 0.5], font=dict(size=20))
plt.xticks([0.10, 0.30, 0.50, 0.70], font=dict(size=20))
plt.tick_params(axis='both', width=2, length=12)
legend = plt.legend(title='# of Balls', loc='upper left', shadow=True, ncol=1, prop={'size': 12})
legend.get_title().set_fontsize(str(15))
plt.title('Open vs Closed Density', fontsize=25)
plt.tight_layout()
plt.show()
