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


num_samples = 1000
x_values = np.linspace(0, 3, num_samples)
# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(8, 5))

for plot_num in (1, 2, 3):


    # Generate random numbers
    random_numbers1 = inverse_transform_sampling(p1, x_values, num_samples)
    random_numbers2 = inverse_transform_sampling(p2, x_values, num_samples)
    random_numbers3 = inverse_transform_sampling(p3, x_values, num_samples)


    if plot_num == 1:
        samples = random_numbers1
        pdf = p1
        title = 'Devries'.format(num_samples)
    elif plot_num == 2:
        samples = random_numbers2
        pdf = p2
        title = 'Ranadive & Lemilch'.format(num_samples)
    elif plot_num == 3:
        samples = random_numbers3
        pdf = p3
        title = 'Gal-Or & Hoelsher'.format(num_samples)




    # Plot the PDF and the histogram on the primary y-axis
    ax1.plot(x_values, pdf(x_values), label=title)


    # Set the limits of the primary y-axis
    # ax1.set_ylim()



# Display the plot
ax1.set_xlabel('Bubble Radius', fontsize=20)
ax1.set_xticks(np.arange(0, 3, 0.5))
ax1.set_ylabel('Probability Density', fontsize=20)
ax1.tick_params('x', labelsize=20)
ax1.set_yticks([])
plt.title('Physics-Based PDFs', fontsize=25)
plt.tight_layout()
plt.legend(fontsize=20)
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
