import matplotlib.pyplot as plt
import numpy as np


def create_graph(title='', xlabel='', ylabel='', xtick_labels=[], ytick_labels=[], title_size=20, xlabel_size=20,
                 ylabel_size=20, xtick_label_size=20, y_tick_label_size=20, tick_width=2, tick_length=12):
    # Create a blank figure and axis
    fig, ax = plt.subplots()

    # Setting the ticks on the y-axis
    ax.set_yticks(ytick_labels)
    ax.set_yticklabels(np.arange(0, 1, len(ytick_labels)), fontsize=y_tick_label_size)
    ax.tick_params(axis='y', length=tick_length, width=tick_width)

    # Setting the ticks on the x-axis
    ax.set_xticks(xtick_labels)
    ax.set_xticklabels(np.arange(0, 1, len(xtick_labels)), fontsize=xtick_label_size)
    ax.tick_params(axis='x', length=tick_length, width=tick_width)

    # Setting labels for axes
    ax.set_xlabel(xlabel, fontsize=xlabel_size)
    ax.set_ylabel(ylabel, fontsize=ylabel_size)

    # Set the title
    ax.title(title, fontsize=title_size)

    # Setting the limit for x and y axis for visibility
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Turn off right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    # Display the plot
    plt.show()


if __name__ == "__main__":
    create_graph()
