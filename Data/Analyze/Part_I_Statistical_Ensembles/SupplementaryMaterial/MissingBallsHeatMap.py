import os
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.system import System


def num_missing_balls(folder=None):
    # If the folder option isnt chosen prompt the user to choose a folder
    if folder is None:
        # Get the folder with all the logs
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes()
        folder = filedialog.askdirectory(title="Choose a data folder")

    # Create the densities data
    missing_balls = {}
    # Loop through the folders
    for subfolder in os.listdir(folder):

        # Get the cv and density values
        split_subfolder = subfolder.split("_")
        try:
            sf_cv, sf_den = float(split_subfolder[1]), float(split_subfolder[3])
        except:
            continue

        # Get the pdb, aw, and pow
        pdb, aw, pow = get_files(folder + '/' + subfolder)
        # Get the logs and make a system
        try:
            my_aw = read_logs2(aw)
        except TypeError:
            print(pdb, aw, pow)
            continue
        ther_balls_counter = 0
        # Loop through the balls
        for i, ball in my_aw['atoms'].iterrows():

            # Make sure the ball is complete
            if ball['Complete Cell?']:
                ther_balls_counter += 1
        # Create the missing balls counter
        if (sf_cv, sf_den) in missing_balls:
            missing_balls[(sf_cv, sf_den)][0] += ther_balls_counter
            missing_balls[(sf_cv, sf_den)][1] += 1000
        else:
            missing_balls[(sf_cv, sf_den)] = [ther_balls_counter, 1000]

    return missing_balls


def plot_heatmap(missing_balls, data=None, cv_values=None, density_values=None):
    if data is None:
        # Extract unique CV and density values from missing_balls
        cv_values = sorted({key[0] for key in missing_balls.keys()})
        density_values = sorted({key[1] for key in missing_balls.keys()}, reverse=True)

        # Create a 2D array for percentages
        data = np.zeros((len(density_values), len(cv_values)))

        # Populate the data array with percentages from missing_balls
        for (cv, density), (ther_balls_counter, total_balls) in missing_balls.items():
            i = density_values.index(density)
            j = cv_values.index(cv)
            percentage = (ther_balls_counter / total_balls) * 100
            data[i, j] = percentage
    print(cv_values)
    print(density_values)
    print(data)
    # Mask for values of 0
    masked_data = np.ma.masked_where(data == 0, data)

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    heatmap = ax.imshow(masked_data, cmap="RdYlGn", aspect="auto", vmin=0, vmax=100)

    # Add color bar
    cbar = plt.colorbar(heatmap)
    cbar.set_label("Percentage (%)", fontsize=25)
    cbar.ax.tick_params(labelsize=20, size=10, width=2, length=12)

    # Set axis labels and title
    ax.set_xlabel("CV", fontsize=30)
    ax.set_ylabel("Density", fontsize=30)
    ax.set_title("Percentage of Complete Cells", fontsize=30)

    # Set axis tick labels
    ax.set_xticks(np.arange(len(cv_values)))
    ax.set_yticks(np.arange(len(density_values)))
    ax.set_xticklabels([f"{v:0.1f}" for v in cv_values], fontsize=20)
    ax.set_yticklabels([f"{v:.2f}" for v in density_values], fontsize=20)
    ax.tick_params('both', length=12, width=2)

    # Rotate x-axis tick labels for better readability
    # plt.xticks(rotation=45)

    # Annotate each cell with the percentage value
    for i in range(len(density_values)):
        for j in range(len(cv_values)):
            if data[i, j] == 0:
                ax.text(j, i, "N/A", ha="center", va="center", color="black", fontsize=15)
            else:
                ax.text(j, i, f"{data[i, j]:.0f}", ha="center", va="center", color="black", fontsize=15)

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    os.chdir('../../../..')
    m_blizzys = num_missing_balls()
    cv_values, density_values, data = None, None, None
    # cv_values = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    # density_values = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
    # data = np.array([[79.09473684, 78.95263158, 79.52105263, 79.24, 79.82, 79.355, 79.335, 79.17, 78.965   ,   78.76  ,     78.84      ], [78.12105263 ,78.05789474 ,78.3     ,   78.115   ,   78.17 ,      78.245 , 78.01    ,   77.97     ,  77.64   ,    77.22631579, 76.6       ] ,[76.83684211 ,76.77894737, 76.83157895, 76.76   ,    77.185  ,    77.035 , 76.89   ,    76.2 ,       75.85789474 ,75.7 ,       75.115     ], [74.96842105, 75.06842105, 75.07894737 ,74.965    ,  75.085   ,   75.335 , 75.245  ,    74.775  ,    74.17894737 ,73.43684211 ,73.465     ], [73.48421053 ,73.6   ,     73.42631579, 73.575  ,    73.485 ,     73.8 , 73.495    ,  72.61    ,   72.54   ,    71.82631579 ,71.535     ], [71.25789474 ,71.24736842, 70.91052632 ,71.055  ,    71.435     , 71.345  ,71.455  ,    71.28    ,   70.41    ,   69.87368421, 69.22      ], [69.58888889, 69.20526316 ,68.62105263 ,68.495 ,     68.99     ,  69.455,  68.9    ,    68.53    ,   67.945   ,   67.91052632 ,66.57      ] ,[ 0.   ,      66.48947368, 65.685    ,  65.645    ,  65.88    ,   66.125 , 66.11052632, 65.865   ,   65.42     ,  64.72105263 ,63.84      ] ,[ 0.    ,      0.  ,       63.08     ,  62.635    ,  62.84    ,   63.295 , 62.68     ,  62.7     ,   62.58   ,    61.93157895 ,60.405     ] ,[ 0.   ,       0.    ,      0.    ,     59.26666667 ,59.41,       59.595,  59.535,      59.29,       59.03,       58.1,        57.46315789]])
    # data =np.array([[100.,    100.,    100.,    100.,     99.955,  99.695,  99.02,   97.965,  96.61,  95.075,  93.095],
    #                 [100.,    100.,    100.,     99.995,  99.985,  99.68,   99.225,  97.96,   96.61,  95.03,   93.645],
    #                 [100.,    100.,    100.,    100.,     99.955,  99.69,   99.14,   98.235,  97.015,  95.695,  94.05 ],
    #                 [100.,    100.,    100.,     99.985,  99.93,   99.69,   99.075,  98.05,   97.42,   96.24,   94.88 ],
    #                 [100.,    100.,    100.,    100.,     99.905,  99.705,  99.18,   98.28,   97.505,  96.61,   95.05 ],
    #                 [100.,    100.,    100.,     99.995,  99.93,   99.785,  99.225,  98.42,   97.645,  96.695,  95.79 ],
    #                 [100.,    100.,    100.,    99.995,  99.94,   99.67,   99.345,  98.635,  98.025,  97.23,   96.87 ],
    #                 [100.,    100.,    100.,    100.,     99.955,  99.75,   99.42,   98.945,  98.375,  97.82,   97.575],
    #                 [100.,    100.,     99.99,  100.,     99.94,   99.895,  99.495,  99.195,  98.94,   98.55,   98.235],
    #                 [99.99,   99.975,  99.97,   99.965,  99.96,   99.875,  99.79,   99.49,   99.415,  99.31,   99.055]])
    plot_heatmap(m_blizzys, data=data, cv_values=cv_values, density_values=density_values)
