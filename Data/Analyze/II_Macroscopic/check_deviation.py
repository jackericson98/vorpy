import csv
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
my_file = filedialog.askopenfilename()

with open(my_file, 'r') as read_file:
    csv_read = csv.reader(read_file)
    data = []
    for line in csv_read:
        data.append(float(line[6]))

data = np.array(data)

# Calculate the overall mean
overall_mean = np.mean(data)


# Function to calculate standard deviation
def calculate_standard_deviation(subset):
    return np.std(subset) / np.sqrt(len(subset))


# Store average standard deviations for each subset size
average_std_devs = []

# Calculate deviations for each combination of subset sizes
for subset_size in range(2, 21):
    std_devs = []
    # Generate all combinations of the current subset size
    for subset in combinations(data, subset_size):
        # Calculate standard deviation of the subset
        std_dev = calculate_standard_deviation(subset)
        std_devs.append(std_dev)

    # Calculate the average standard deviation for the current subset size
    average_std_devs.append(np.mean(std_devs))

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(range(2, 21), average_std_devs, marker='o')
plt.title('Average Standard Deviation vs. Number of Data Points')
plt.xlabel('Number of Data Points in Subset')
plt.ylabel('Average Standard Deviation')
plt.grid(True)
plt.show()