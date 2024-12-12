import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

# Normalizing function for PDFs
def normalize_pdf(pdf, range):
    integral, _ = quad(pdf, *range)
    return lambda x: pdf(x) / integral

# Unnormalized PDFs
def p1_unnormalized(r):
    return 2.082 * r / (1 + 0.387 * r**2)**4

def p2_unnormalized(r):
    return (32 / np.pi**2) * r**2 * np.exp(-(4 / np.pi) * r**2)

def p3_unnormalized(x):
    return (16 / np.pi) * x**2 * np.exp(-np.sqrt(16 / np.pi) * x**2)

# Normalized PDFs
p1 = normalize_pdf(p1_unnormalized, (0, np.inf))
p2 = normalize_pdf(p2_unnormalized, (0, np.inf))
p3 = normalize_pdf(p3_unnormalized, (0, np.inf))

# Calculate the cumulative distribution function (CDF)
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

# Parameters
num_samples = 1000
x_values = np.linspace(0, 3, num_samples)

# Define common bins
common_bins = np.linspace(0, 3, 31)  # 30 equally spaced bins

# Plot setup
fig, ax1 = plt.subplots(figsize=(7, 4))
colors = ['purple', 'orange', 'green']  # Distinguishable colors
titles = ['Dev', 'R & L', 'G & H']
pdfs = [p1, p2, p3]
alpha = 0.8

for i, (pdf, title, color) in enumerate(zip(pdfs, titles, colors), 1):
    alpha -= 0.2
    # Generate random samples
    random_numbers = inverse_transform_sampling(pdf, x_values, num_samples)

    # Plot normalized PDF
    ax1.plot(x_values, pdf(x_values), label=title, color=color, linewidth=4)

    # Plot histogram with common bins
    plt.hist(random_numbers, bins=common_bins, density=True, alpha=alpha, color=color, edgecolor='k', label=f'{title} Hist')

    # Finalize plot
ax1.set_xlabel('Bubble Radius', fontsize=25)
ax1.set_ylabel('Probability', fontsize=25)
ax1.tick_params('both', labelsize=25)
ax1.set_yticks([])
plt.title('Physics-Based PDFs', fontsize=30)
plt.legend(fontsize=15)
plt.tight_layout()
plt.show()
