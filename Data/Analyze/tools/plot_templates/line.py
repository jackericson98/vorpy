import matplotlib.pyplot as plt
import numpy as np


def line_plot(xs, ys, errors=None, labels=None, error_alpha=0.2, title=None, x_label=None, y_label=None,
              legend_title=None, Show=True, title_size=25, x_label_size=20, y_label_size=20, legend_title_size=20,
              legend_label_size=20, tick_width=2, tick_length=12, legend_orientation='vertical'):
    # Create a single plot
    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(len(xs)):
        if labels is not None:
            ax.plot(xs[i], ys[i], label=labels[i])
        else:
            ax.plot(xs[i], ys[i])
        if errors is not None:
            ax.fill_between(xs[i], [ys[i][j] - errors[i][j] for j in range(len(ys[i]))],
                            [ys[i][j] + errors[i][j] for j in range(len(ys[i]))], alpha=error_alpha)

    # Set plot title and legend
    # ax.set_xticks(np.arange(xs[0][1], xs[0][-1] + 0.05, 0.05))
    if title is not None:
        ax.set_title(title, fontsize=title_size)
    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=x_label_size)
    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=y_label_size)
    ax.tick_params(axis='both', which='major', labelsize=legend_label_size, length=tick_length, width=tick_width)
    if labels is not None:
        ncol = 1
        if legend_orientation == 'horizontal':
            ncol = len(labels)
        legend = ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1), prop={'size': legend_label_size}, ncol=ncol,
                           shadow=True)
        if legend_title is not None:
            legend.set_title(legend_title)
            legend.get_title().set_fontsize(str(legend_title_size))

    # Adjust the right margin to make room for the legend
    plt.subplots_adjust(right=0.8)
    plt.tight_layout()
    # Show the plot
    if Show:
        plt.show()
