import csv
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.plot_templates.histogram import histogram
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory as blend
from matplotlib.ticker import LogLocator


def distribution_of_overlaps(file=None, output_folder=None, bins=10):
    # Check if the file is nothing
    if file is None:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        file = filedialog.askopenfilename(title="Overlap Data")
    # if output_folder is None:
    #     output_folder = filedialog.askdirectory(title="Output Directory")
    # Open the file
    with open(file, 'r') as my_file:
        # Create the csv file
        my_c = csv.reader(my_file)
        # Create the dictionary to store the data
        my_dict = {}
        # Set the maximum value
        my_max_val = 0
        # Loop through the lines
        for line in my_c:
            # Print the data progress


            # Skip the zero lines
            if len(line) == 0:
                continue
            # Split up the file name so that it is in terms of density and cv
            file_name_split = line[0].split('/')
            file_name_split_further = file_name_split[-1].split('_')
            if len(file_name_split_further) < 4:
                continue
            cv, density, number = file_name_split_further[1], file_name_split_further[3], file_name_split_further[-1]
            olap_value = float(file_name_split_further[4])
            if number == 'True' or number == 'False':
                number = '0'
            # Add the new file to the dictionary if not in there
            if (cv, density, number) not in my_dict:
                # Create the file dictionary in the data dictionary
                vals =  [float(_) for _ in line[2:]]
                my_dict[(cv, density, number)] = {line[1]: vals}
                if max(vals) > my_max_val:
                    my_max_val = max(vals)
            # If the file is in the dictionary keep adding the other ball values
            else:
                vals = [float(_) for _ in line[2:]]
                my_dict[(cv, density, number)][line[1]] = vals
                if max(vals) > my_max_val:
                    my_max_val = max(vals)
    # Create the data storing dictionary
    new_data_dict = {}
    # Create a list of cv vals and density vals
    cv_vals, density_vals = [], []
    # Go through each file name and grab the density and cv data
    for cv, density, number in my_dict:
        # Go through the balls in the file
        for ball in my_dict[(cv, density, number)]:
            # # Separate the non zeros
            # non_zeros = [_ for _ in my_dict[(cv, density, number)][ball] if _ != 0]
            # # Count the zeros
            # num_zeros = len(my_dict[(cv, density, number)][ball]) - len(non_zeros)
            # Add the data
            if cv in new_data_dict:
                if density in new_data_dict[cv]:
                    new_data_dict[cv][density] += my_dict[(cv, density, number)][ball]
                # new_data_dict[(cv, density)]['num zeros'] += num_zeros
                else:
                    new_data_dict[cv][density] = my_dict[(cv, density, number)][ball]
            else:
                new_data_dict[cv] = {density: my_dict[(cv, density, number)][ball]}
                if cv not in cv_vals:
                    cv_vals.append(cv)
                if density not in density_vals:
                    density_vals.append(density)

    # Sort the density values
    density_vals.sort()
    cv_vals.sort()
    # Create a 10x11 figure
    fig, axes = plt.subplots(10, 11, figsize=(20, 18), sharex='all', sharey='all')
    # axes = axes.flatten()  # Flatten for easy indexing

    # Function to create a histogram with a y-axis break
    def plot_histogram_with_break(ax, data, bins, ybreak=10000):
        hist, edges = np.histogram(data, bins=bins)
        heights = hist.copy()

        # Handle y-axis break
        heights = np.where(heights > ybreak, ybreak, heights)

        ax.bar(edges[:-1], heights, width=np.diff(edges), align="edge", color="blue", alpha=0.7, edgecolor="black")

        # Add a break indicator
        if max(hist) > ybreak:
            transform = blend(ax.transData, ax.transAxes)
            ax.set_yscale('log')
            ax.plot(
                [edges[0] + np.diff(edges)[0] * 0.2, edges[0] + np.diff(edges)[0] * 0.8],
                [ybreak, ybreak],
                linestyle="--",
                color="black",
                transform=transform,
                clip_on=False,
            )
            # ax.yaxis.set_major_locator(LogLocator(base=10.0))
    # Plot each set
    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            # Grab the axes
            ax = axes[i, j]
            # Get the data
            data = new_data_dict[cv][density]
            # Bin the data
            plot_histogram_with_break(ax, data, np.linspace(0.0, olap_value, bins))
            ax.set_xticks([])
            ax.set_yticks([])
            # ax.hist(data, bins=bins, color='blue', alpha=0.7, edgecolor='black')

    # Add figure-wide x and y axes labels
    for ax, col in zip(axes[-1], cv_vals):
        ax.set_xlabel("")  # Remove direct subplot labels, handled below

    for ax, row in zip(axes[:, 0], density_vals):
        ax.set_ylabel("")  # Remove direct subplot labels, handled below

    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.15 + i * (0.8 / len(cv_vals)), 0.07, col, ha='center', fontsize=15)

    # Add density values to the left of the figure
    for i, row in enumerate(density_vals):
        fig.text(0.07, 0.15 + i * (0.8 / len(density_vals)), row, va='center', rotation='horizontal', fontsize=15)

    fig.text(0.5, 0.04, 'CV Values', ha='center', fontsize=15)  # X-axis label for the figure
    fig.text(0.04, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=15)  # Y-axis label for the figure

    # Add a main title for the entire figure
    fig.suptitle("Distribution of Overlaps by CV and Density", fontsize=20)

    # Adjust layout to prevent overlapping labels
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.95)
    # plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    # Run the code
    distribution_of_overlaps()