import matplotlib.pyplot as plt
import numpy as np


def bar(data, errors=None, x_names=None, legend_names=None, title='', x_axis_title='', y_axis_title='', bar_width=0.35,
        Show=False, save=None, legend_title=None, print_vals_on_bars=False, unit='', title_size=25, tick_width=2,
        tick_length=12, xlabel_size=20, ylabel_size=20, xtick_label_size=20, ytick_label_size=20, legend_entry_size=20,
        x_range=None, y_range=None, legend_orientation='Horizontal', x_tick_rotation=0, y_tick_rotation=0):

    # Check how the data is set up and make sure it is a list of lists
    if type(data[0]) is not list:
        data = [data]

    # Get the total maximum for the list of lists
    err_max = 0
    if errors is not None:
        err_max = max([max(_) for _ in errors])
    ymax = max([max(_) for _ in data]) + err_max
    ymin = min([min(_) for _ in data]) - err_max

    # Get the number of bars to plot
    num_bars = len(data)

    # Get the number of bar groups to plot
    num_groups = range(len(data[0]))

    # Get the names for the individual bars
    if legend_names is None or len(legend_names) != len(data):
        legend_names = ['' for _ in range(len(data))]

    # Get the names for each group of bars
    if x_names is None or len(x_names) != len(data[0]):
        x_names = ['' for _ in range(len(data[0]))]

    # Set the colors
    colors = ['skyblue', 'orange', 'lavender', 'red', 'goldenrod', 'slategray', 'rose', 'coral', 'periwinkle',
              'turquoise']

    # Plot the bars
    xlocs = [[i * bar_width + j * bar_width * (num_bars + 1) for j in num_groups] for i in range(len(data))]
    for i in range(len(data)):
        my_bars = plt.bar(xlocs[i], data[i], width=bar_width, label=legend_names[i], color=colors[i], edgecolor='black')

        # Plot the error bars
        if errors is not None:
            for j, my_bar in enumerate(my_bars):
                plt.errorbar(my_bar.get_x() + my_bar.get_width() / 2, data[i][j], yerr=errors[i][j], capsize=5, capthick=2,
                             color='black', alpha=0.8)

    # Plot the title, ylabel and xlabel
    plt.title(title, fontdict=dict(size=title_size))
    plt.ylabel(y_axis_title, fontdict=dict(size=ylabel_size))
    plt.xlabel(x_axis_title, fontdict=dict(size=xlabel_size))

    # Label the bar groups
    x_tick_locs = [np.mean([xlocs[j][i] for j in range(len(data))]) for i in range(len(data[0]))]
    plt.xticks(x_tick_locs, x_names, font=dict(size=xtick_label_size), rotation=x_tick_rotation)
    plt.yticks(font=dict(size=ytick_label_size), rotation=y_tick_rotation)
    plt.tick_params(axis='both', width=tick_width, length=tick_length)

    # Plot the data on the bars
    if print_vals_on_bars:
        for j in range(len(data)):
            for i, v in enumerate(data[j]):
                if v < 0.3 * ymax:
                    height = 0.5 * ymax
                else:
                    height = v / 2
                plt.text(xlocs[j][i] + 0.03, height, str(v) + unit, ha='center', va='center', rotation=90)

    # Add the legend
    leg_col = 1
    if legend_orientation == 'Horizontal':
        leg_col = len(data)
    try:
        if legend_title is not None:
            plt.legend(title=legend_title, loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=leg_col)
        elif len(data) > 1:
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=leg_col,
                       prop={'size': legend_entry_size})
    except Exception:
        pass

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

    # Show the plot if chosen to
    if Show:
        plt.show()

    # Save the graph
    if save is not None:
        plt.savefig(save)
