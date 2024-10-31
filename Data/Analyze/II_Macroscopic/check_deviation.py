import csv
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import tkinter as tk
from tkinter import filedialog
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_file = filedialog.askopenfilename()

with open(my_file, 'r') as read_file:
    csv_read = csv.reader(read_file)
    data = {}
    my_sds = []
    for line in csv_read:
        cv, den = line[3], line[5]
        my_sds.append(cv)
        if (cv, den) in data:
            data[(cv, den)].append(float(line[6]))
        else:
            data[(cv, den)] = [float(line[6])]


# Function to calculate standard deviation
def calculate_standard_error(subset):
    return np.std(subset) / np.sqrt(len(subset))


def calculate_standard_deviation(subset):
    return np.std(subset)


def coefficient_of_variation(subset):
    return np.std(subset) / np.mean(subset)


# Plotting the results
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
ax1, ax2, ax3 = axes

cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
norm = Normalize(vmin=min(my_sds), vmax=max(my_sds))
sm = ScalarMappable(norm=norm, cmap=cmap)

for _ in data:
    my_data = np.array(data[_])

    # Calculate the overall mean
    overall_mean = np.mean(my_data)

    # Store average standard deviations for each subset size
    average_std_devs, average_std_errs, average_cvs = [], [], []

    # Calculate deviations for each combination of subset sizes
    for subset_size in range(2, len(my_data) + 1):
        std_devs, std_errs, cvs = [], [], []
        meannys = []
        # Generate all combinations of the current subset size
        for subset in combinations(my_data, subset_size):
            # Calculate standard deviation of the subset
            std_dev = calculate_standard_deviation(subset)
            std_devs.append(std_dev)
            # Calculate the standard error
            std_err = calculate_standard_error(subset)
            std_errs.append(std_err)
            # Calculate the coeffient of variation
            coef_var = coefficient_of_variation(subset)
            cvs.append(coef_var)
            meannys.append(np.mean(subset))

        # Calculate the average standard deviation for the current subset size
        average_std_devs.append(np.mean(std_devs))
        average_std_errs.append(np.mean(std_errs))
        average_cvs.append(np.mean(cvs))
        print([(meannys[i], std_devs[i]) for i in range(len(meannys))])
        ax3.plot([subset_size for _ in range(len(meannys))], meannys)
    ax1.plot(range(2, len(my_data) + 1), average_std_devs, c=cmap(norm(float(_[0]))))
    ax2.plot(range(2, len(my_data) + 1), average_std_errs, c=cmap(norm(float(_[0]))))
    # ax3.plot(range(2, len(my_data) + 1), average_cvs, c=cmap(norm(float(_[0]))))
    # Adding a color bar that uses the created ScalarMappable
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))
ax1.set_title('Average Standard Deviation vs. Number of Data Points')
ax2.set_title('Average Standard Error vs. Number of Data Points')
ax3.set_title('Average Coefficient of Variation vs. Number of Data Points')
ax1.set_xlabel('Number of Data Points in Subset')
ax2.set_xlabel('Number of Data Points in Subset')
ax3.set_xlabel('Number of Data Points in Subset')
ax1.set_ylabel('Average Standard Deviation')
ax2.set_ylabel('Average Standard Error')
ax3.set_ylabel('Average Coefficient of Variation')
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
plt.show()
