import os
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from Data.Analyze.tools.batch.get_files import get_files
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from Data.Analyze.tools.compare.read_logs import read_logs
from System.system import System


def num_missing_balls(folder=None):
    # If the folder option isnt chosen prompt the user to choose a folder
    if folder is None:
        # Get the folder with all the logs
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes()
        folder = filedialog.askdirectory(title="Choose a data folder")
    print(folder)
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
        if pdb is None or aw is None or pow is None:
            continue
        silly = False
        try:

            my_aw = read_logs2(aw)
        except:
            my_aw = read_logs(aw)
            silly = True
        ther_balls_counter = 0
        # Loop through the balls
        for i, ball in my_aw['atoms'].iterrows():

            # Make sure the ball is complete
            if not silly and ball['Complete Cell?']:
                ther_balls_counter += 1
            if silly and ball['complete']:
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
            print(cv, density, ther_balls_counter, total_balls)
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
    # m_blizzys = num_missing_balls()
    m_blizzys = None
    # cv_values, density_values, data = None, None, None
    cv_values = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    density_values = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
    # nonPBC 0.5
    # data = np.array([[79.09473684, 78.95263158, 79.52105263, 79.24, 79.82, 79.355, 79.335, 79.17, 78.965   ,   78.76  ,     78.84      ], [78.12105263 ,78.05789474 ,78.3     ,   78.115   ,   78.17 ,      78.245 , 78.01    ,   77.97     ,  77.64   ,    77.22631579, 76.6       ] ,[76.83684211 ,76.77894737, 76.83157895, 76.76   ,    77.185  ,    77.035 , 76.89   ,    76.2 ,       75.85789474 ,75.7 ,       75.115     ], [74.96842105, 75.06842105, 75.07894737 ,74.965    ,  75.085   ,   75.335 , 75.245  ,    74.775  ,    74.17894737 ,73.43684211 ,73.465     ], [73.48421053 ,73.6   ,     73.42631579, 73.575  ,    73.485 ,     73.8 , 73.495    ,  72.61    ,   72.54   ,    71.82631579 ,71.535     ], [71.25789474 ,71.24736842, 70.91052632 ,71.055  ,    71.435     , 71.345  ,71.455  ,    71.28    ,   70.41    ,   69.87368421, 69.22      ], [69.58888889, 69.20526316 ,68.62105263 ,68.495 ,     68.99     ,  69.455,  68.9    ,    68.53    ,   67.945   ,   67.91052632 ,66.57      ] ,[ 0.   ,      66.48947368, 65.685    ,  65.645    ,  65.88    ,   66.125 , 66.11052632, 65.865   ,   65.42     ,  64.72105263 ,63.84      ] ,[ 0.    ,      0.  ,       63.08     ,  62.635    ,  62.84    ,   63.295 , 62.68     ,  62.7     ,   62.58   ,    61.93157895 ,60.405     ] ,[ 0.   ,       0.    ,      0.    ,     59.26666667 ,59.41,       59.595,  59.535,      59.29,       59.03,       58.1,        57.46315789]])
    # PBC 0.0
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
    # PBC 0.5
    # data = np.array([[100.0, 100.0, 100.0, 99.97, 99.85, 99.545, 98.84, 97.555, 96.305, 94.95, 93.39],
    #                  [100.0, 100.0, 100.0, 100.0, 99.88, 99.34, 98.955, 97.685, 96.305, 95.39, 93.805],
    #                  [100.0, 100.0, 100.0, 99.965, 99.86, 99.6, 98.825, 97.96, 96.59, 96.06, 94.205],
    #                  [100.0, 100.0, 100.0, 99.98, 99.795, 99.375, 98.805, 97.94, 96.995, 96.04, 95.27],
    #                  [100.0, 100.0, 100.0, 99.955, 99.765, 99.42, 98.955, 97.955, 97.05, 96.275, 95.315],
    #                  [100.0, 100.0, 100.0, 100.0, 99.82, 99.39, 99.13, 98.395, 97.81, 96.53, 95.99],
    #                  [100.0, 100.0, 100.0, 99.93, 99.745, 99.56, 99.19, 98.65, 97.59, 97.575, 96.865],
    #                  [100.0, 100.0, 100.0, 99.925, 99.875, 99.7, 99.385, 99.015, 98.59, 98.24, 97.84],
    #                  [100.0, 100.0, 100.0, 99.97, 99.87, 99.755, 99.48, 99.26, 99.155, 98.805, 98.735],
    #                  [100.0, 100.0, 100.0, 99.885, 99.85, 99.76, 99.615, 98.885, 98.845, 98.725, 98.885]])
#     print('nonpbc 1.0')
    # nonPBC 1.0
    # data = np.array([
    #                 [61.11764706, 60.49375, 59.57692308, 58.4125, 56.86666667, 55.89166667, 54.18, 52.54285714, 50.40666667, 49.87058824, 46.09],
    #                 [61.09285714, 61.0125, 59.77142857, 58.69333333, 58.00714286, 57.03571429, 55.76666667, 53.82352941, 53.32142857, 51.15, 48.675],
    #                 [61.05625, 61.1, 60.35, 59.43333333, 58.40909091, 57.93846154, 56.91538462, 55.95882353, 54.03846154, 52.76111111, 50.305],
    #                 [61.02142857, 61.17692308, 60.24166667, 59.88666667, 59.13076923, 58.40909091, 57.08823529, 56.38461538, 55.9875, 54.66666667, 52.96],
    #                 [61.225, 61.05384615, 60.8625, 59.61538462, 59.54666667, 58.81333333, 57.8, 57.54285714, 56.36666667, 56.47142857, 55.15],
    #                 [61.13, 60.42, 60.13571429, 59.86, 59.7, 58.84615385, 58.57777778, 58.49333333, 57.47222222, 57.36666667, 56.375],
    #                 [60.59230769, 59.84545455, 60.01818182, 60.17333333, 59.64375, 59.03333333, 59.16363636, 58.43529412, 58.59285714, 59.0, 58.535],
    #                 [60.24, 60.43076923, 60.14, 59.16923077, 60.06, 59.84375, 59.26153846, 59.82666667, 59.52352941, 59.75625, 60.21],
    #                 [60.1, 59.73125, 59.45714286, 59.7, 59.55384615, 59.89411765, 60.62142857, 61.17857143, 62.14, 62.79230769, 63.47058824],
    #                 [59.67142857, 60.26666667, 59.84666667, 61.12666667, 61.575, 63.37333333, 64.1, 65.2125, 66.775, 66.3875, 68.255]])
    # PBC 1.0
    # data = np.array([
    #     [100.000, 100.000, 100.000, 99.895, 99.540, 98.920, 98.240, 97.185, 95.490, 95.005, 92.855],
    #     [100.000, 100.000, 100.000, 99.870, 99.525, 98.955, 97.995, 96.925, 95.480, 94.870, 93.300],
    #     [100.000, 100.000, 99.985, 99.925, 99.675, 98.900, 98.400, 97.200, 95.960, 95.155, 93.655],
    #     [100.000, 100.000, 99.985, 99.845, 99.600, 98.850, 98.190, 97.200, 96.115, 95.160, 94.765],
    #     [100.000, 100.000, 99.995, 99.840, 99.695, 98.945, 98.255, 97.420, 96.920, 95.735, 94.935],
    #     [100.000, 100.000, 99.995, 99.915, 99.600, 99.050, 98.630, 97.820, 96.665, 96.075, 95.450],
    #     [100.000, 100.000, 99.985, 99.860, 99.545, 99.250, 98.645, 98.190, 97.335, 96.635, 95.980],
    #     [100.000, 100.000, 99.970, 99.830, 99.630, 99.175, 98.730, 98.375, 97.750, 97.375, 96.750],
    #     [100.000, 100.000, 99.985, 99.870, 99.705, 99.395, 98.995, 98.650, 98.295, 97.990, 97.680],
    #     [100.000, 100.000, 99.960, 99.930, 99.790, 99.670, 99.400, 99.160, 99.015, 98.860, 98.710]
    # ])

    plot_heatmap(m_blizzys, data, cv_values, density_values)
