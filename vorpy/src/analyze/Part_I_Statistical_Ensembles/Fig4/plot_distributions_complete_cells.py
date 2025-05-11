import os

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.special import gamma as gamma_func
from scipy.optimize import fsolve
import tkinter as tk
from tkinter import filedialog
from vorpy.src.analyze.tools.compare.read_logs import read_logs
from vorpy.src.inputs import read_pdb_line
from vorpy.src.system.system import System


cv = '0.5'
density = '0.5'
open_cell = True


def lognormal(r, cv, mu=1):
    sigma = np.sqrt(np.log(cv**2 + 1))
    mu_log = np.log(mu / np.sqrt(1 + cv**2))  # Adjusted to incorporate mu
    lognormal_dist = stats.lognorm(s=sigma, scale=np.exp(mu_log))
    return lognormal_dist.pdf(r)


def gamma(r, cv, mu=1):
    # Gamma parameters
    alpha = 1 / cv ** 2
    beta = alpha / mu  # To keep mean = 1
    gamma_dist = stats.gamma(a=alpha, scale=1 / beta)
    # Compute PDFs
    return gamma_dist.pdf(r)


def weibull(r, cv, mu=1):


    # Define equations to solve for kappa and lambda
    def equations(p):
        kappa, lambda_ = p
        mean_eq = lambda_ * gamma_func(1 + 1/kappa) - mu
        var_eq = lambda_**2 * (gamma_func(1 + 2/kappa) - gamma_func(1 + 1/kappa)**2) - (cv * mu)**2
        return (mean_eq, var_eq)

    # Initial guesses for kappa and lambda
    kappa_initial = 0.75
    lambda_initial = mu

    # Solve for kappa and lambda
    kappa, lambda_ = fsolve(equations, (kappa_initial, lambda_initial))

    # Create Weibull distribution
    weibull_dist = stats.weibull_min(c=kappa, scale=lambda_)
    return weibull_dist.pdf(r)


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
# Get the root directory
root_dir = filedialog.askdirectory()

directories_we_want = []
# Get the pdb file, and the two logs files
for sub_dir, dirs, files in os.walk(root_dir):
    for direc in dirs:
        my_cv, my_den = direc[4:8], direc[14:18]
        if my_cv[-1] == '_':
            my_cv = my_cv[:-1]
            my_den = direc[13:17]
        if my_den[-1] == '_':
            my_den = my_den[:-1]
        if my_den == density and my_cv == cv:
            directories_we_want.append(root_dir + '/' + direc)

pdb_files, pow_logs, vor_logs = [], [], []
# Go through the directories we want and add the logs
for direct in directories_we_want:
    for sub_dir, dirs, files in os.walk(direct):
        # get the logs
        if sub_dir[-2:] == 'aw':
            vor_logs.append(sub_dir + '/' + files[0])
        elif sub_dir[-3:] == 'pow':
            pow_logs.append(sub_dir + '/' + files[0])
        else:
            pdb_files.append(sub_dir + '/' + files[0])


print(len(pdb_files), len(vor_logs), len(pow_logs))

# Change the directory so the system doesn't crap
os.chdir('../../../..')

rads, complete_rads = [], []
for i, pdb in enumerate(pdb_files):
    # Interpret the logs
    try:
        pow_ = read_logs(pow_logs[i])
        vor_ = read_logs(vor_logs[i])
        sys_ = System(pdb)
    except IndexError:
        continue
    # Loop through the lines
    for j, atom in sys_.balls.iterrows():
        # Get the power atom line
        try:
            pow_atom = pow_['atoms'].loc[pow_['atoms']['num'] == j, 'complete'].to_list()[0]
            vor_atom = vor_['atoms'].loc[vor_['atoms']['num'] == j, 'complete'].to_list()[0]
            rads.append(atom['rad'])
        except IndexError:
            pass
        complete_rads.append(atom['rad'])





# Get the full folder
#
# # Take in the PDB
# pdb = filedialog.askopenfilename(title='Get the pdb file mf')
# # Get the power logs
# pow_logs = filedialog.askopenfilename(title='Get the pow logs mf')
# # Get the Voronoit Logs
# vor_logs = filedialog.askopenfilename(title='Get the vor logs mf')







# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(10, 5))

# Get the x values for the pdf
x_values = np.linspace(0, 5, 100)

# Plot the PDF and the histogram on the primary y-axis
ax1.plot(x_values, gamma(x_values, 0.5), label='PDF', color='blue')
ax1.set_xlabel('Bubble Radius', fontsize=20)
ax1.set_xticks(np.arange(0, 5, 0.5))
ax1.set_ylabel('Probability Density Function', color='blue', fontsize=20)
ax1.tick_params('y', colors='blue', labelsize=15)
ax1.tick_params('x', labelsize=15)

# Set the limits of the primary y-axis
ax1.set_ylim(bottom=0)

# Create a secondary y-axis for the histogram
ax2 = ax1.twinx()

# Plot the histogram on the secondary y-axis
ax2.hist(complete_rads, bins=50, alpha=0.4, color='blue')
ax2.hist(rads, bins=50, alpha=0.4, color='red')
ax2.set_ylabel('Number of Bubbles', color='red', fontsize=20)
ax2.tick_params('y', colors='red', labelsize=15)

# Set the limits of the secondary y-axis
ax2.set_ylim(bottom=0)

cell_type = 'Non-Overlapping'
if open_cell:
    cell_type = 'Overlapping'

# Display the plot
plt.title('Gamma {}, Density = {}, CV = {} rads'.format(cell_type, density, cv), fontsize=25)
plt.show()
