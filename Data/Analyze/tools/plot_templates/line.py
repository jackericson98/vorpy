import matplotlib.pyplot as plt
import numpy as np


def line_plot(xs, ys, errors=None, labels=None, error_alpha=0.2, title=None, x_label=None, y_label=None,
              legend_title=None, Show=True):
    # Create a single plot
    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(len(xs)):
        ax.plot(xs[i], ys[i], label=labels[i])
        if errors is not None:
            ax.fill_between(xs[i], [ys[i][j] - errors[i][j] for j in range(len(ys[i]))],
                            [ys[i][j] + errors[i][j] for j in range(len(ys[i]))], alpha=error_alpha)

    # Set plot title and legend
    # ax.set_xticks(np.arange(xs[0][1], xs[0][-1] + 0.05, 0.05))
    if title is not None:
        ax.set_title(title, fontsize=30)
    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=20)
    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=15)
    if labels is not None:
        legend = ax.legend(loc='upper right', bbox_to_anchor=(1.3, 0.8))
    if legend_title is not None:
        legend.set_title(legend_title)

    # Adjust the right margin to make room for the legend
    plt.subplots_adjust(right=0.8)

    # Show the plot
    if Show:
        plt.show()
