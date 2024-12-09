import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from Data.Analyze.tools.batch.compile_new_logs import get_logs_and_pdbs
from Data.Analyze.tools.compare.read_logs import read_logs
from Data.Analyze.tools.compare.read_logs2 import read_logs2
from System.system import System
from scipy import stats
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d


def calculate_cdf(pdf, x_values):
    cdf_values = np.array([quad(pdf, 0, x)[0] for x in x_values])
    cdf_values /= cdf_values[-1]  # Normalize to [0, 1]
    return cdf_values


# Generate random samples using inverse transform sampling
def inverse_transform_sampling(pdf, x_values, n_samples):
    cdf_values = calculate_cdf(pdf, x_values)
    inverse_cdf = interp1d(cdf_values, x_values, kind='linear', fill_value='extrapolate')
    u = np.random.rand(n_samples)
    return inverse_cdf(u)


def get_hist_of_radii(density=0.5, cv=0.5, bins=100):
    # First pick the folder to find the cv and density
    # root = tk.Tk()
    # root.withdraw()
    # root.wm_attributes('-topmost', 1)
    # complete_folder = filedialog.askdirectory(title="New Data")
    # old_folder = filedialog.askdirectory(title="Old Data")
    def gamma(r, mu=1):
        # Gamma parameters
        alpha = 1 / cv ** 2
        beta = alpha / mu  # To keep mean = 1
        gamma_dist = stats.gamma(a=alpha, scale=1 / beta)
        # Compute PDFs
        return gamma_dist.pdf(r)
    # # Look through the different folders to get all logs and pdbs that match the set cv and density
    # new_logs_pdbs = get_logs_and_pdbs(True, output_file_name='loggy_woggys_new.txt')
    old_logs_pdbs = get_logs_and_pdbs(make_file=True, output_file_name='../../tools/batch/loggy_woggys_old.txt')

    # Loop through the new_logs
    # Record the radii
    radii = []
    real_radii = []
    for _ in old_logs_pdbs:
        values = _.split('_')
        if values[1] != str(cv) or values[3] != str(density):
            continue
        # Create the system
        my_sys = System(old_logs_pdbs[_]['pdb'], simple=True)
        # Get the logs
        aw_logs = read_logs(old_logs_pdbs[_]['aw'])
        pow_logs = read_logs(old_logs_pdbs[_]['pow'])

        for i, ball in my_sys.balls.iterrows():
            if i > 999:
                continue
            real_radii.append(ball['rad'])
            # Get the aw_ball and the pow ball
            try:
                aw_ball = aw_logs['atoms'][aw_logs['atoms']['num'] == ball['num']].iloc[0].to_dict()
                pow_ball = pow_logs['atoms'][pow_logs['atoms']['num'] == ball['num']].iloc[0].to_dict()
            except ValueError:
                continue
            except IndexError:
                continue
            radii.append(ball['rad'])


    # Create a figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Real Radii

    data1 = ax1.hist(real_radii, bins=bins, alpha=0.5, color='red', edgecolor='k')

    data2 = ax1.hist(radii, bins=bins, alpha=0.5, color='blue', edgecolor='k')
    diffs1, diffs2 = [], []
    x_values = np.linspace(0, 5, bins)
    for i, val in enumerate(x_values):
        gamma_val = gamma(val)

        diffs1.append(100 * abs(len(real_radii) * gamma(val) - data1[0][i]) / len(real_radii) * gamma(val))

        diffs2.append(100 * abs(len(radii) * gamma(val) - data2[0][i]) / len(radii) * gamma(val))
    print(diffs1)
    print(diffs2)
    ax1.set_xlabel('Bubble Radius')
    ax1.set_xticks(np.arange(0, 5, 0.5))

    # ax1.set_ylabel('Probability Density Function', color='blue')
    # ax1.tick_params('y', colors='blue')

    # Set the limits of the primary y-axis
    ax1.set_ylim(bottom=0)

    # Create a secondary y-axis for the histogram
    ax2 = ax1.twinx()

    ax2.plot(x_values, len(real_radii) * gamma(x_values), label='PDF', color='red', alpha=0.5)
    ax2.plot(x_values, len(radii) * gamma(x_values), label='PDF', color='blue', alpha=0.5)
    # Plot the histogram on the secondary y-axis
    # ax2.set_ylabel('Number of Bubbles', color='red')
    # ax2.tick_params('y', colors='red', )

    # Set the limits of the secondary y-axis
    ax2.set_ylim(bottom=0)

    # Display the plot
    # plt.title(title)
    plt.show()

    # plt.plot(x_values, [diffs2[i] - diffs2])

    plt.plot(x_values, diffs2)
    plt.plot(x_values, diffs1)
    plt.show()


if __name__ == '__main__':
    get_hist_of_radii()
