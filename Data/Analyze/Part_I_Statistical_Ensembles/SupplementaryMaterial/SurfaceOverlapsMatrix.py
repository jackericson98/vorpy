import csv
import os
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import LogLocator, FuncFormatter


from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs
from Data.Analyze.tools.compare.read_logs2 import read_logs2
import pandas as pd
from System.sys_funcs.calcs.calcs import calc_dist
from System.system import System


def get_overlap_data(file_name):
    # First gather the log files
    my_logs = get_logs_and_pdbs(False)
    # Create the overlaps directory
    olaps = {}
    log_length = len(my_logs)
    # Loop through each of the logs files looking for the overlap data
    for j, loggy in enumerate(my_logs):
        # print the progress
        print("\rReading Log {}/{} - {:.2f}%".format(j + 1, log_length, 100 * (j + 1) / log_length), end="")
        # split the loggy value
        loggy_list = loggy.split('_')
        # Ge the cv and density
        cv, density = float(loggy_list[1]), float(loggy_list[3])
        # Get the files
        try:
            pdb_file, aw_logs_file, pow_logs_file = [my_logs[loggy][_] for _ in ['pdb', 'aw', 'pow']]
        except KeyError:
            continue
        # Get the logs dictionaries
        aw_logs = read_logs2(aw_logs_file, True, all_=False, balls=True)
        # System of balls
        my_sys = System(pdb_file, simple=True)
        if 'rad' not in my_sys.balls:
            continue
        # Create the dataframe for the logs
        aw_dataframe = pd.DataFrame(aw_logs['atoms'])
        # pow_logs = read_logs2(pow_logs_file, True, all_=False, balls=True)
        for i, ball in aw_dataframe.iterrows():
            # Loop through the balls
            overlaps = []
            # Loop through each of the neighbors looking for size
            for neighbor in ball['Neighbors']:
                # Get the neighbors rad and loc
                nrad, nloc = my_sys.balls['rad'][neighbor], my_sys.balls['loc'][neighbor]
                # Check the size
                if nrad > ball['rad']:
                    # Add the overlap percentage
                    overlaps.append(round(max(nrad + ball['rad'] - calc_dist(ball['loc'], nloc), 0) / ball['rad'], 5))
            # # Write the overlap data to match the overlaps file we already have
            with open(file_name, 'a') as my_overlap_file:
                # Create the csv_writer
                olap_csv = csv.writer(my_overlap_file)
                # Write the new line
                olap_csv.writerow([loggy, i] + overlaps)


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
            if len(line) <= 2:
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
                vals = [float(_) for _ in line[2:]]
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
            if cv not in new_data_dict:
                new_data_dict[cv] = {}
                if cv not in cv_vals:
                    cv_vals.append(cv)
            if density not in new_data_dict[cv]:
                new_data_dict[cv][density] = []
                if density not in density_vals:
                    density_vals.append(density)
            new_data_dict[cv][density] += my_dict[(cv, density, number)][ball]

    density_vals.sort(reverse=True)
    cv_vals.sort()
    print(density_vals, cv_vals)

    fig, axes = plt.subplots(10, 11, figsize=(20, 18), sharex='all', sharey='all')

    def format_ticks(val, pos):
        """Format tick labels as powers of ten."""
        if val == 0:
            return "$10^0$"
        else:
            exponent = int(np.log10(val))
            return f"$10^{{{exponent}}}$"

    for i, density in enumerate(density_vals):
        for j, cv in enumerate(cv_vals):
            ax = axes[i, j]
            data = new_data_dict[cv][density]

            my_data = ax.hist(data, bins, range=(0, olap_value), density=False, log=True, color='cyan', edgecolor='k')
            ax.set_yscale('log')
            ax.set_ylim([1, 1000000])
            if j == 0:  # Far-left subplots
                ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))  # Half number of ticks
                ax.yaxis.set_minor_formatter(FuncFormatter(format_ticks))  # Remove scientific notation
                ax.yaxis.set_major_formatter(FuncFormatter(format_ticks))  # Remove scientific notation
            else:
                ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=[], numticks=6))
                ax.yaxis.set_major_formatter(FuncFormatter(format_ticks))


    for ax, col in zip(axes[-1], cv_vals):
        ax.set_xlabel("")  # Remove direct subplot labels, handled below

    for ax, row in zip(axes[:, 0], density_vals):
        ax.set_ylabel("")  # Remove direct subplot labels, handled below

    # Add CV values to the bottom of the figure
    for i, col in enumerate(cv_vals):
        fig.text(0.13 + i * (0.815 / len(cv_vals)), 0.05, col, ha='center', fontsize=15)

    for i, row in enumerate(density_vals[::-1]):
        fig.text(0.05, 0.125 + i * (0.825 / len(density_vals)), row, va='center', rotation='horizontal', fontsize=15)

    fig.text(0.5, 0.02, 'CV Values', ha='center', fontsize=20)  # X-axis label for the figure
    fig.text(0.02, 0.5, 'Density Values', va='center', rotation='vertical', fontsize=20)  # Y-axis label for the figure

    # Add a main title for the entire figure
    fig.suptitle("Distribution of Overlaps by CV and Density", fontsize=20)

    # Adjust layout to prevent overlapping labels
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9)
    # plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    os.chdir('../../../..')
    file_naaame = 'overlaps_1.csv'
    get_overlap_data(file_naaame)
    # Run the code
    distribution_of_overlaps(file_naaame)
