import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy import stats
from plotly import graph_objects as go


# Define the probability distribution function
def p1(r):
    return 2.082 * r / (1 + 0.387 * r**2)**4


# Define the probability distribution function
def p2(r):
    return (32/np.pi**2) * r**2 * np.exp(-(4/np.pi) * r**2)


# Define the probability density function (PDF)
def p3(x):
    # Example PDF (you should replace this with your own PDF)
    return (16/np.pi) * x**2 * np.exp(-np.sqrt(16/np.pi) * x**2)


def gamma(r):
    return 0.5 * stats.gamma.pdf(r, 4, scale=1 / 6.5)


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


plot_num = 3
num_samples = 5000

x_values = np.linspace(0, 5, num_samples)


# Generate random numbers
random_numbers1 = inverse_transform_sampling(p1, x_values, num_samples)
random_numbers2 = inverse_transform_sampling(p2, x_values, num_samples)
random_numbers3 = inverse_transform_sampling(p3, x_values, num_samples)
random_numbers_gamma = inverse_transform_sampling(gamma, x_values, num_samples)


if plot_num == 1:
    samples = random_numbers1
    pdf = p1
    title = 'Devries Bubble Distribution - N = {}'.format(num_samples)
elif plot_num == 2:
    samples = random_numbers2
    pdf = p2
    title = 'Ranadive & Lemilch Bubble Distribution - N = {}'.format(num_samples)
elif plot_num == 3:
    samples = random_numbers3
    pdf = p3
    title = 'Gal-Or & Hoelsher Bubble Distribution - N = {}'.format(num_samples)
elif plot_num == 'gamma':
    samples = random_numbers_gamma
    pdf = gamma
    title = 'gamma'

# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(10, 5))

# Plot the PDF and the histogram on the primary y-axis
ax1.plot(x_values, pdf(x_values), label='PDF', color='blue')
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
ax2.hist(samples, bins=30, alpha=0.5, color='red')
ax2.set_ylabel('Number of Bubbles', color='red', fontsize=20)
ax2.tick_params('y', colors='red', labelsize=15)

# Set the limits of the secondary y-axis
ax2.set_ylim(bottom=0)

# Display the plot
plt.title(title, fontsize=25)
plt.show()
#
# # Plot the histogram of generated random numbers
# plt.hist(random_numbers3, bins=50, density=True, alpha=0.5, label='Generated Data')
#
#
# plt.plot(x_data, [p3(_) for _ in x_data])
#
# plt.xlabel('r')
# plt.ylabel('Probability Density')
# plt.title('Random Numbers from Given Distribution')
# plt.legend()
# plt.show()
