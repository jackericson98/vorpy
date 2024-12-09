import csv
import numpy as np

import matplotlib.pyplot as plt


def convert_point(string_point):
    vals = string_point.split(',')
    return float(vals[0][1:]), float(vals[1]), float(vals[2][:-1])


tracker = {}
with open('location_data.csv', 'r') as my_data:
    reader = csv.reader(my_data)
    for line in reader:
        if len(line) == 0:
            continue
        cv, density = line[:2]
        if (cv, density) not in tracker:
            tracker[(cv, density)] = []

        for point in line[2:]:
            tracker[(cv, density)].append(np.array(convert_point(point)))

for cv, density in tracker:
    if cv != '0.05' or density != '0.05':
        continue
    fig, ax = plt.figure(figsize=(8, 6)), plt.gca()
    points = np.array(tracker[(cv, density)])
    scatter = ax.scatter(points[:, 0], points[:, 1], c=points[:, 2], cmap='gray', s=0.8, marker='x')
    ax.set_xlabel('X value', fontsize=25)
    ax.set_ylabel('Y value', fontsize=25)
    ax.set_title('Distribution of Coordinates\n(cv={}, d={})'.format(cv, density), fontsize=30)

    # Adjusting color bar
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label('Z value', fontsize=25)
    cbar.ax.tick_params(labelsize=20)  # Adjust color bar tick font size

    # Adjusting x and y tick font sizes
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.tight_layout()
    plt.show()
