import matplotlib.pyplot as plt
import numpy as np


def histogram(samples, num_bins=30, pdf=None, pdf_res=1000, fig=None, ax=None, bar_color=None, pdf_color=None, title='', xlabel='',
              ylabel='', pdf_ylabel=''):
    """
    Takes in a list of values and bins them based on the num_bins variable, plotting a Histogram.
    """
    if fig is None:
        # Create a figure and axis
        fig = plt.figure(figsize=(10, 5))
    if bar_color is None:
        bar_color = 'blue'
    if pdf_color is None:
        pdf_color = 'red'

    # Create the axis
    if ax is None:
        ax = fig.add_subplot()
    # Allow for different histogram plots
    if type(samples[0]) is not list:
        samples = [samples]
        bar_color = [bar_color]

    for i, sample_set in enumerate(samples):
        # Plot the histogram on the secondary y-axis
        ax.hist(samples, bins=num_bins, alpha=0.5, color='red')
        ax.set_ylabel(ylabel, color='red', fontsize=20)
        ax.tick_params('y', colors='red', labelsize=15)

    # Set the limits of the primary y-axis
    ax.set_ylim(bottom=0)

    if pdf is not None:
        # create the x values
        x_values = np.linspace(start=min(samples), stop=max(samples), num=pdf_res)
        # Create a secondary y-axis for the histogram
        ax2 = ax.twinx()
        # Plot the PDF and the histogram on the primary y-axis
        ax2.plot(x_values, pdf(x_values), label='PDF', color=pdf_color)
        ax2.set_xlabel(xlabel, fontsize=20)
        ax2.set_ylabel(pdf_ylabel, color=pdf_color, fontsize=20)
        ax2.tick_params('y', colors=pdf_color, labelsize=15)
        ax2.tick_params('x', labelsize=15)
        # Set the limits of the secondary y-axis
        ax2.set_ylim(bottom=0)

    # Display the plot
    plt.title(title, fontsize=25)
    plt.show()