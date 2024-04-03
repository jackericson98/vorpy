import matplotlib.pyplot as plt


def box_whisker(data, x_names=None, title='', x_axis_title='', y_axis_title='', Show=False, save=None):

    plt.boxplot(data, labels=x_names)

    # Plot the title, ylabel and xlabel
    plt.title(title, fontdict=dict(size=20))
    plt.ylabel(y_axis_title, fontdict=dict(size=15))
    plt.xlabel(x_axis_title, fontdict=dict(size=15))

    # Label the bar groups
    # x_locs = [j * bar_width * (num_bars + 1) + num_bars * bar_width / 2 for j in num_groups]
    # plt.xticks(x_locs, x_names, rotation=45, ha='right', font=dict(size=10))

    # Show the plot if chosen to
    if Show:
        plt.show()

    # Save the graph
    if save is not None:
        plt.savefig(save)

