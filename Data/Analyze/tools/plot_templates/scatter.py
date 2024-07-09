import matplotlib.pyplot as plt
import numpy as np


def scatter(xs, ys, labels=None, errors=None, title='', x_axis_title='', y_axis_title='', Show=False, save=None,
            legend_title=None, unit='', title_size=25, tick_width=2, tick_length=12, xlabel_size=20, ylabel_size=20,
            xtick_label_size=20, ytick_label_size=20, legend_entry_size=20, x_range=None, y_range=None,
            legend_orientation='Horizontal', alpha=1, legend_title_size=20):

    # Set the colors
    colors = ['skyblue', 'orange', 'red', 'blue', 'green', 'pink', 'y', 'purple']

    # Set the markers
    markers = ['.', ',', 'o', 's', 'v', '^', '<', '>', '1', '2', '3', '4', '*', '+', 'x', 'D', 'd', 'p', 'h', 'H', '|',
               '_']

    # Create a single plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the actual data
    for i in range(len(xs)):
        if labels is not None:
            ax.scatter(xs[i], ys[i], marker=markers[i], label=labels[i], c=colors[i], alpha=alpha)
        else:
            ax.scatter(xs[i], ys[i], marker=markers[i], c=colors[i], alpha=alpha)
    # Check how the data is set up and make sure it is a list of lists
    if type(xs[0]) is not list:
        xs = [xs]
        ys = [ys]

    # Get the total maximum for the list of lists
    err_max = 0
    if errors is not None:
        err_max = max([max(_) for _ in errors])
    ymax = max([max(_) for _ in ys]) + err_max
    ymin = min([min(_) for _ in ys]) - err_max

    # Get the number of bar groups to plot
    num_groups = range(len(ys[0]))

    # Plot the title, ylabel and xlabel
    # ax.title(title, fontdict=dict(size=title_size))
    plt.ylabel(y_axis_title, fontdict=dict(size=ylabel_size))
    plt.xlabel(x_axis_title, fontdict=dict(size=xlabel_size))

    # Label the bar groups
    plt.xticks(rotation=45, ha='right', font=dict(size=xtick_label_size))
    plt.yticks(font=dict(size=ytick_label_size))
    plt.tick_params(axis='both', width=tick_width, length=tick_length)

    # Add the legend
    leg_col = 1
    if legend_orientation == 'Horizontal':
        leg_col = len(xs)

    if legend_title is not None:
        legend = ax.legend(title=legend_title, loc='upper right', bbox_to_anchor=(1.25, 0.97), shadow=True, ncol=leg_col,
                  prop={'size': legend_entry_size})
        legend.get_title().set_fontsize(str(legend_title_size))
    elif len(xs) > 1:
        legend = ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=leg_col,
                  prop={'size': legend_entry_size})

        legend.get_title().set_fontsize(str(legend_title_size))
    # Set the y limit
    multiplier = 1.3
    if ymin < 0:
        multiplier = 1.5
    if y_range is None:
        plt.ylim(ymin * 1.1, multiplier * ymax)
    else:
        if y_range[0] is not None:
            ymin = y_range[0]
        if y_range[1] is not None:
            ymax = y_range[1]
        plt.ylim(ymin, multiplier * ymax)

    # Set the x limits
    if x_range is not None:
        plt.xlim(*x_range)

    # Set the figure size
    plt.tight_layout()

    if Show:
        plt.show()
