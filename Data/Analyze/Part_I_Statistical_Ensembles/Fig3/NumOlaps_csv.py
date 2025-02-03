import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import tkinter as tk
from tkinter import filedialog
import csv


def parse_data(filepath):
    # Dictionary to hold the parsed data
    data_dict = {}

    # Open and read the file
    with open(filepath, 'r') as file1:
        file_csv = csv.reader(file1)
        for line in file_csv:

            # Extract the filename from the path
            try:
                path_parts = line[0].split('/')
            except IndexError:
                continue
            filename = path_parts[-1]

            # Extract details from the filename
            try:
                details = filename.split('_')

                mean = float(details[0])
                cv = float(details[1])
                num_balls = int(details[2])
                density = float(details[3].replace('_', '.'))
            except IndexError:
                print(line)
                continue
            except ValueError:
                print(line)
                continue
            try:
                file_number = int(details[-1])
            except ValueError:
                file_number = 0


            # Key for the main dictionary
            key = (cv, density, file_number)

            # Second column as key for the subdictionary
            sub_key = int(line[1])

            # Remaining columns as list of values
            values = [float(col) if col.isdigit() or '.' in col else int(col) for col in line[2:]]

            # Insert data into the dictionary
            if key not in data_dict:
                data_dict[key] = {}
            data_dict[key][sub_key] = values

    return data_dict


def plot_overlapping_line_plots(data):
    perdic, cvs, dens = {}, [], []
    for cv, density, file_number in data:
        curdic = data[(cv, density, file_number)]
        percent = 100 * len([_ for _ in curdic if sum(curdic[_]) != 0]) / len(curdic)
        if (cv, density) in perdic:
            perdic[(cv, density)].append(percent)
        else:
            perdic[(cv, density)] = [percent]
        if cv not in cvs:
            cvs.append(cv)
        if density not in dens:
            dens.append(density)
    cvs.sort()
    dens.sort()
    data_means, data_std_up, data_std_dn = [[] for _ in range(len(cvs))], [[] for _ in range(len(cvs))], [[] for _ in range(len(cvs))]
    for i, num in enumerate(cvs):
        for j, num2 in enumerate(dens):
            currdata = perdic[(num, num2)]
            mean = np.mean(currdata)
            std = np.std(currdata)
            data_means[i].append(mean)
            data_std_dn[i].append(mean - std)
            data_std_up[i].append(mean + std)

    # Coefficient of Variation (CV) and Density values
    cmap = plt.cm.rainbow  # Choose a colormap that does not have yellow and works well in grayscale
    norm = Normalize(vmin=min(cvs), vmax=max(cvs))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots(figsize=(8, 12))

    for i, sd in enumerate(cvs):
        if sd == 0.05:
            continue
        try:
            # Colors for each line based on 'sd' which is used as an index into the colormap
            color = cmap(norm(sd))
            # if round(sd, 4) > 0.5:
            #     datavvm[i] = [0.9 * _ for _ in datavvm[i]]
            #     datavvms[i], datavvps[i] = [0.9 * _ for _ in datavvms[i]], [0.9 * _ for _ in datavvps[i]]
            # elif 0.35 < sd < 0.45:
            #     datavvm[i] = [1.1 * _ for _ in datavvm[i]]
            #     datavvms[i], datavvps[i] = [1.1 * _ for _ in datavvms[i]], [1.1 * _ for _ in datavvps[i]]
            # if round(sd, 4) <= 0.5:
            #     datavsm[i] = [0.78 * _ for _ in datavsm[i]]
            #     datavsms[i], datavsps[i] = [0.78 * _ for _ in datavsms[i]], [0.78 * _ for _ in datavsps[i]]
            # # elif 0.35 < sd < 0.45:
            # #     datavvm[i] = [1.1 * _ for _ in datavvm[i]]
            # #     datavvms[i], datavvps[i] = [1.1 * _ for _ in datavvms[i]], [1.1 * _ for _ in datavvps[i]]
            # if value == 'vol':
            ax.plot(dens, data_means[i], color=color)
            ax.fill_between(dens, data_std_dn[i], data_std_up[i], color=color, alpha=0.2)
            # elif value == 'sa':
            #     ax.plot(my_densities, datavsm[i], color=color)
            #     ax.fill_between(my_densities, datavsms[i], datavsps[i], color=color, alpha=0.2)
        except ValueError:
            continue

    # Adding a color bar that uses the created ScalarMappable
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Coefficient of Variation (CV)', fontdict=dict(size=25))

    # Set plot titles and labels
    ax.set_xticks(np.arange(dens[0] + 0.05, dens[-1] + 0.05, 0.1))
    # if value == 'vol':
    #     ax.set_ylim([0, 200])
    # elif value == 'sa':
    #     ax.set_ylim([0, 50])
    # ax.set_title('{} Power {}\nAbsolute % Difference'
    #              .format('Overlapping' if cell_type == 'Open' else 'Non-Overlapping',
    #                      {'sa': 'Surface Area', 'vol': 'Volume'}[value]), fontsize=20)
    ax.set_xlabel('Density', fontsize=25)
    ax.set_ylabel('% Balls Overlapping', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=12)

    cbar.ax.tick_params(labelsize=20, size=10, width=2, length=12)
    plt.tight_layout()

        # Show the plot
    plt.show()

# Example usage:
# data = parse_data('path_to_your_file.txt')
# print(data)


if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    file = filedialog.askopenfilename()

    my_data = parse_data(file)
    plot_overlapping_line_plots(my_data)

