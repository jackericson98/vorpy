import matplotlib.pyplot as plt


def bar(data, errors=None, x_names=None, legend_names=None, title='', x_axis_title='', y_axis_title='', bar_width=0.35, Show=False,
        save=None, legend_title=None):

    # Check how the data is set up and make sure it is a list of lists
    if type(data[0]) is not list:
        data = [data]

    # Get the total maximum for the list of lists
    ymax = max([max(_) for _ in data])

    # Set the bar width
    bar_width = 0.35

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
    colors = ['skyblue', 'orange', 'lavender', 'mintgreen', 'peach', 'goldenrod', 'slategray', 'rose', 'coral', 'periwinkle',
              'turquoise']

    # Plot the bars
    for i in range(len(data)):
        x_locs = [i * bar_width + j * bar_width * (num_bars + 1) for j in num_groups]
        my_bars = plt.bar(x_locs, data[i], width=bar_width, label=legend_names[i], color=colors[i], edgecolor='black')

        # Plot the error bars
        if errors is not None:
            for j, my_bar in enumerate(my_bars):
                plt.errorbar(my_bar.get_x() + my_bar.get_width() / 2, data[i][j], yerr=errors[i][j], capsize=5, capthick=2,
                             color='black', alpha=0.8)

    # Plot the title, ylabel and xlabel
    plt.title(title, fontdict=dict(size=20))
    plt.ylabel(y_axis_title, fontdict=dict(size=15))
    plt.xlabel(x_axis_title, fontdict=dict(size=15))

    # Label the bar groups
    x_locs = [j * bar_width * (num_bars + 1) + num_bars * bar_width / 2 for j in num_groups]
    plt.xticks(x_locs, x_names, rotation=45, ha='right', font=dict(size=10))

    # Add the legend
    if legend_title is not None:
        plt.legend(title=legend_title)
    elif len(data) > 1:
        plt.legend()

    # Set the y limit
    plt.ylim(0, 1.3 * ymax)

    # Show the plot if chosen to
    if Show:
        plt.show()

    # Save the graph
    if save is not None:
        plt.savefig(save)
