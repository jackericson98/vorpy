import csv
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, FuncFormatter


def format_ticks(val, pos):
    """Format tick labels as powers of ten."""
    if val == 0:
        return "$10^0$"
    else:
        exponent = int(np.log10(val))
        return f"$10^{{{exponent}}}$"


def distribution_of_overlaps(cv, den, bins=30):
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    file_path = filedialog.askopenfilename(title="Select Overlap Data File")
    data_dict = {}

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        max_value = 0
        for line in reader:
            if len(line) > 2:
                parts = line[0].split('/')
                filename_parts = parts[-1].split('_')
                if len(filename_parts) >= 4:
                    file_cv, file_den, number = float(filename_parts[1]), float(filename_parts[3]), filename_parts[-1]
                    overlap_value = float(filename_parts[4])

                    if file_cv == cv and file_den == den:
                        number = number if number not in ['True', 'False'] else '0'
                        key = (cv, file_den, number)

                        if key not in data_dict:
                            data_dict[key] = {}

                        values = [float(x) for x in line[2:]]
                        data_dict[key][line[1]] = values
                        if max(values) > max_value:
                            max_value = max(values)

    # Prepare data for plotting
    aggregated_data = {}
    density_vals, cv_vals = set(), set()

    for (cv_key, den_key, _), balls in data_dict.items():
        for ball_id, values in balls.items():
            if cv_key not in aggregated_data:
                aggregated_data[cv_key] = {}
            if den_key not in aggregated_data[cv_key]:
                aggregated_data[cv_key][den_key] = []
            aggregated_data[cv_key][den_key].extend(values)
            density_vals.add(den_key)
            cv_vals.add(cv_key)

    density_vals = sorted(density_vals, reverse=True)
    cv_vals = sorted(cv_vals)

    fig, ax = plt.subplots()
    for density in density_vals:
        for cv in cv_vals:
            data = aggregated_data[cv][density]
            my_hist = ax.hist(data, bins, range=(0, overlap_value), density=False, log=True, color='cyan', edgecolor='black')
            ax.set_yscale('log')
            ax.set_ylim([1, 1000000])
            # Enhanced LogLocator for major and minor ticks
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=10))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=100))
            ax.yaxis.set_major_formatter(FuncFormatter(format_ticks))
            ax.tick_params(axis='both', which='major', labelsize=20)  # Set tick label size
            ax.tick_params(axis='both', which='minor', labelsize=10)  # Set minor tick label size
            ax.set_ylabel("Surface Count", fontsize=25)
            ax.set_xlabel("Small Radius Overlap", fontsize=25)
            ax.set_title("Distribution of Overlaps", fontsize=25)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    distribution_of_overlaps(cv=1.0, den=0.05)
