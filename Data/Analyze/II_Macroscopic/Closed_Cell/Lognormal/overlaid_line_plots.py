import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_foams_file = askopenfilename()


plot_type = 'lognormal'



with open(my_foams_file, 'r') as my_foam_data:

    my_data = []
    # x, y, z1, z2 = [], [], [], []
    foam_data = csv.reader(my_foam_data)
    for my_line in foam_data:
        if len(my_line) <= 1:
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
    my_densities = np.arange(0.025, 0.45, 0.025)
    my_sds = np.arange(1.5, 6.5, 0.25)

if plot_type == 'lognormal':
    my_densities = np.arange(0.025, 0.525, 0.025)
    my_sds = np.arange(0.1, 0.5, 0.025)

my_densities = [round(_, 3) for _ in my_densities]
my_sds = [round(_, 3) for _ in my_sds]

# Initialize lists using list comprehensions
datavvm, datavsm, datapvm, datapsm = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvms, datavsms, datapvms, datapsms = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]
datavvps, datavsps, datapvps, datapsps = [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds], [[] for _ in my_sds]


# Iterate over my_sds and my_densities using nested loops
for i, num in enumerate(my_sds):
    for j, num2 in enumerate(my_densities):
        means = []
        sds = []
        # Get the current Data
        try:
            curr_data = lists[num][num2]
        except KeyError:
            print(num, num2)
            continue

        for data1 in curr_data:
            # Calculate the Z-scores
            z_scores = np.abs((data1 - np.mean(data1)) / np.std(data1))

            # Set a Z-score threshold (e.g., 3)
            z_score_threshold = 1

            # Exclude outliers based on the Z-score
            filtered_data = np.array(data1)[z_scores < z_score_threshold]
            means.append(100*np.mean(filtered_data))
            sds.append(np.std([_ * 100 for _ in filtered_data]) / np.sqrt(len(filtered_data)))

        # Append means and mean +/- SD to respective lists
        datavvm[i].append(means[0]); datavvms[i].append(means[0] - sds[0]); datavvps[i].append(means[0] + sds[0])
        datavsm[i].append(means[1]); datavsms[i].append(means[1] - sds[1]); datavsps[i].append(means[1] + sds[1])
        datapvm[i].append(means[2]); datapvms[i].append(means[2] - sds[2]); datapvps[i].append(means[2] + sds[2])
        datapsm[i].append(means[3]); datapsms[i].append(means[3] - sds[3]); datapsps[i].append(means[3] + sds[3])


# datavvmlog = [[9.514, 9.613, 9.424, 9.488, 8.606, 8.857, 8.926, 8.653, 8.862, 8.746, 8.461, 8.396, 8.184, 8.333, 8.74, 7.671, 8.26, 8.501],
# [12.068, 12.087, 12.515, 12.528, 11.759, 12.806, 11.03, 11.569, 11.278, 11.491, 11.467, 11.556, 10.997, 11.684, 11.044, 10.929, 10.995, 10.55],
# [14.56, 14.68, 14.618, 14.697, 14.484, 15.301, 14.717, 15.014, 14.158, 14.3, 14.537, 14.118, 15.033, 13.748, 14.26, 14.665, 13.576, 14.433],
# [16.576, 17.969, 17.906, 18.444, 18.107, 17.971, 17.855, 18.99, 18.69, 19.662, 18.29, 18.549, 18.887, 17.918, 18.199, 17.742, 19.46, 19.177],
# [19.4, 20.621, 21.278, 20.699, 23.003, 23.035, 22.177, 22.415, 23.891, 24.836, 24.612, 24.257, 22.736, 24.682, 22.932, 23.43, 22.964, 22.37],
# [22.751, 24.378, 25.718, 25.479, 27.282, 26.433, 26.293, 29.072, 27.936, 27.686, 29.754, 28.97, 28.379, 28.997, 28.423, 30.629, 30.602, 28.336],
# [24.13, 26.412, 28.109, 30.077, 29.495, 34.316, 32.356, 33.548, 36.653, 37.862, 34.131, 36.568, 35.976, 37.878, 37.005, 36.528, 38.166, 38.744],
# [26.282, 31.517, 33.797, 39.083, 38.334, 38.882, 40.984, 41.259, 42.497, 43.552, 46.539, 44.291, 47.563, 44.823, 47.189, 49.768, 48.315, 48.197],
# [27.584, 36.294, 36.315, 38.581, 38.392, 39.693, 45.036, 45.572, 51.458, 49.203, 55.718, 52.182, 57.392, 56.729, 53.913, 58.574, 59.015, 60.586],
# [31.096, 37.991, 40.483, 42.63, 48.566, 52.0, 58.42, 61.689, 61.154, 62.93, 65.834, 65.854, 68.374, 66.131, 67.239, 69.02, 70.044, 72.937],
# [32.526, 39.705, 44.198, 47.102, 49.282, 52.982, 55.366, 61.796, 62.594, 62.827, 66.14, 68.023, 70.068, 70.514, 73.095, 75.729, 76.073, 80.276],
# [35.119, 43.485, 46.352, 48.093, 54.015, 57.208, 60.187, 69.62, 74.752, 75.869, 77.032, 79.856, 81.643, 85.003, 90.629, 93.995, 102.253, 101.066],
# [34.773, 45.995, 53.169, 55.023, 61.667, 66.056, 69.702, 71.806, 79.406, 83.234, 83.241, 83.033, 88.117, 92.361, 101.48, 104.308, 106.533, 108.552],
# [38.968, 46.883, 51.942, 59.307, 62.42, 68.35, 70.554, 74.252, 80.744, 90.618, 97.26, 102.934, 107.285, 108.3, 108.426, 114.437, 117.147, 123.434],
# [36.184, 47.667, 57.161, 65.001, 70.04, 77.403, 80.636, 89.792, 90.303, 97.509, 96.533, 108.898, 116.116, 118.33, 125.916, 121.019, 131.29, 134.975],
# [40.917, 49.667, 56.51, 66.5, 71.979, 76.959, 83.165, 87.197, 100.348, 101.89, 103.763, 117.151, 119.656, 122.836, 133.052, 132.914, 140.714, 140.71],
# [41.297, 53.603, 65.555, 72.232, 74.375, 77.438, 92.761, 96.589, 101.159, 114.974, 119.927, 129.13, 133.618, 137.795, 135.935, 143.729, 144.138, 143.515]]
#
# datavvm1 = [
#             [28.230, 31.694, 36.784, 37.442, 38.810, 41.817, 42.029, 43.246, 44.686, 45.668, 46.043, 46.069, 46.237, 46.609, 47.974, 48.198, 49.147],
#             [29.090, 35.767, 36.679, 39.707, 40.490, 41.390, 46.687, 47.503, 48.245, 49.960, 51.580, 51.538, 52.543, 52.724, 52.764, 53.150, 54.083],
#             [29.524, 34.857, 37.530, 40.640, 41.934, 43.616, 45.852, 48.809, 49.081, 51.579, 51.877, 52.244, 54.278, 58.090, 59.699, 60.470, 60.751],
#             [30.801, 36.067, 41.639, 42.293, 43.504, 45.561, 46.587, 46.102, 48.934, 50.302, 53.333, 57.856, 56.230, 60.186, 58.221, 59.833, 62.081],
#             [30.919, 36.814, 39.135, 42.054, 44.796, 46.740, 50.781, 52.829, 54.922, 56.034, 57.723, 56.145, 58.764, 59.721, 60.282, 61.803, 64.005],
#             [31.337, 36.602, 38.311, 42.260, 48.441, 47.283, 46.119, 48.107, 51.410, 49.955, 54.916, 55.034, 60.177, 60.015, 61.856, 63.525, 66.532],
#             [29.373, 36.156, 38.398, 42.187, 45.771, 45.519, 49.012, 52.567, 54.740, 54.008, 54.277, 60.576, 63.733, 65.459, 63.612, 66.333, 67.574],
#             [32.264, 37.715, 40.759, 41.695, 47.494, 46.368, 50.756, 52.663, 54.234, 55.353, 55.910, 57.716, 60.400, 59.054, 61.473, 68.808, 68.159],
#             [30.645, 36.985, 39.525, 45.245, 46.976, 47.909, 47.001, 48.277, 52.084, 57.300, 56.078, 56.659, 60.702, 59.554, 67.842, 68.556, 69.535],
#             [30.913, 37.097, 40.297, 43.398, 46.392, 50.615, 50.588, 53.212, 53.017, 56.112, 57.515, 60.575, 58.062, 60.126, 61.591, 66.308, 71.554],
#             [30.366, 36.206, 41.303, 43.053, 46.163, 52.378, 49.983, 54.898, 56.108, 60.931, 61.543, 59.279, 59.964, 64.065, 67.119, 71.734, 71.963],
#             [29.821, 36.431, 41.585, 42.642, 46.340, 45.696, 55.629, 55.394, 57.396, 59.069, 60.229, 60.526, 64.046, 62.307, 71.413, 71.313, 71.885],
#             [32.375, 36.851, 40.707, 43.817, 46.184, 50.492, 55.076, 55.062, 57.653, 58.697, 56.638, 60.218, 64.235, 64.453, 69.446, 70.206, 71.571],
#             [33.366, 40.229, 41.253, 45.984, 46.396, 47.221, 50.093, 55.031, 58.116, 59.858, 57.242, 60.701, 62.696, 67.621, 73.058, 71.255, 73.955],
#             [30.980, 35.382, 42.945, 46.980, 48.878, 52.718, 51.777, 56.971, 59.681, 56.419, 63.981, 65.707, 71.951, 73.609, 72.174, 74.305, 77.027],
#             [32.303, 37.698, 40.563, 45.321, 48.788, 50.678, 50.546, 58.392, 59.879, 60.837, 65.209, 64.619, 67.468, 64.462, 69.798, 71.794, 71.922],
#             [31.524, 36.498, 42.502, 43.430, 48.392, 46.413, 54.825, 59.254, 53.912, 60.556, 59.604, 64.158, 68.601, 69.827, 69.308, 74.223, 75.878],
#             [32.438, 39.126, 42.513, 44.267, 44.849, 50.834, 55.750, 62.524, 62.035, 66.136, 68.720, 69.423, 68.506, 70.539, 75.949, 75.052, 77.995],
#             [32.477, 38.538, 43.118, 47.872, 47.821, 53.741, 59.020, 61.596, 63.077, 66.054, 68.850, 70.072, 69.438, 72.142, 77.708, 75.040, 79.066],
#             [31.606, 37.301, 43.429, 48.339, 49.143, 49.349, 63.606, 63.777, 65.906, 66.273, 69.162, 71.956, 72.210, 76.017, 78.303, 85.450, 85.231]
#     ]


print(len(my_sds), len(my_densities))

# Coefficient of Variation (CV) and Density values

# Create a single plot
fig, ax = plt.subplots(figsize=(8, 6))

for i in range(len(datavvm)):
    ax.plot(my_densities[1:], datavvm[i][1:], label=str(my_sds[i]))

# Set plot title and legend
ax.set_xticks(np.arange(my_densities[1], my_densities[-1] + 0.05, 0.05))
ax.set_title('AW vs Power % Difference Volume', fontsize=30)
ax.set_xlabel('Density', fontsize=20)
ax.set_ylabel('% Difference', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
legend = ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1))
legend.set_title('CV')

# Adjust the right margin to make room for the legend
plt.subplots_adjust(right=0.8)

# Show the plot
plt.show()
