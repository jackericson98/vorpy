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


def gamma(r, cv, mu=1):
    # Gamma parameters
    alpha = 1 / cv ** 2
    beta = alpha / mu  # To keep mean = 1
    gamma_dist = stats.gamma(a=alpha, scale=1 / beta)
    # Compute PDFs
    return gamma_dist.pdf(r)


root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
file = filedialog.askopenfilename()

my_sys = System(file, simple=True)

file_attributes = file.split('/')
file_attributes = file_attributes[-1].split('_')

cv, den = float(file_attributes[1]), float(file_attributes[3])

# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(10, 8))

# Real Radii
radii = my_sys.balls['rad']

# Plot histogram
data2 = ax1.hist(radii, bins=30, alpha=0.5, color='blue', edgecolor='k', density=False)

# Calculate bin width
bin_edges = data2[1]
bin_width = bin_edges[1] - bin_edges[0]
total_area = sum(data2[0]) * bin_width

# Generate x values for the PDF
x_values = np.linspace(0, 3.5, 100)

# Scale gamma PDF to match the histogram area
pdf_values = gamma(x_values, cv)
scaled_pdf = pdf_values * total_area

# Update plot formatting
ax1.set_xlabel('Bubble Radius', fontsize=25)
ax1.set_ylabel('Count', fontsize=25, color='blue')
ax1.set_xticks(np.arange(0, 5, 0.5))
ax1.tick_params(axis='both', labelsize=20)
ax1.set_ylim(bottom=0)

# Create a secondary y-axis for the histogram
ax2 = ax1.twinx()

# Plot the scaled gamma PDF
ax1.plot(x_values, scaled_pdf, label='Scaled PDF', color='red', linewidth=2)
ax2.set_ylabel('Probability', fontsize=25, color='red')
ax2.tick_params(axis='both', labelsize=20)
ax2.set_ylim(bottom=0, top=100)

# Add title
plt.title(f"Distribution of Radii vs. PDF\nCV = {cv}, Density = {den}", fontsize=30)
plt.tight_layout()
# Display the plot
plt.show()
