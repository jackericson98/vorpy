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
    if cv != '0.5' or density != '0.25':
        continue
    fig, ax = plt.figure(figsize=(8, 6)), plt.gca()
    points = np.array(tracker[(cv, density)])
    scatter = ax.scatter(points[:, 0], points[:, 1], c=points[:, 2], cmap='gray', s=0.8, marker='x')
    ax.set_xlabel('X value')
    ax.set_ylabel('Y value')
    ax.set_title('CV 0.5, Density 0.25')
    fig.colorbar(scatter, ax=ax, label='Z value')

    plt.show()
