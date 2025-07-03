import matplotlib.pyplot as plt
import numpy as np


def line_plot(xs, ys, errors=None, labels=None, error_alpha=0.2, title=None, x_label=None, y_label=None,
              legend_title=None, Show=True, title_size=25, x_label_size=20, y_label_size=20, legend_title_size=20,
              legend_label_size=20, tick_width=2, tick_length=12, legend_orientation='vertical', colors=None,
              tick_val_size=15, linewidth=1, x_ticks=None, x_ticks2=None, y_ticks=None, y_ticks2=None, colorbar=None,
              tight_layout=True, figsize=None):
    # Set default colors if none provided
    if colors is None:
        colors = ['r', 'g', 'b', 'y', 'o', 'pink', 'purple', 'lavender']
    
    # Calculate maximum y value for scaling
    y_max = max([max(_) for _ in ys])
    
    # Create a single plot with specified figure size
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot each line with optional labels and error bars
    for i in range(len(xs)):
        # Plot the main line with or without labels
        if labels is not None:
            ax.plot(xs[i], ys[i], label=labels[i], linewidth=linewidth, c=colors[i])
        else:
            ax.plot(xs[i], ys[i], linewidth=linewidth, c=colors[i])
        
        # Add error bars if provided
        if errors is not None:
            ax.fill_between(xs[i], [ys[i][j] - errors[i][j] for j in range(len(ys[i]))],
                            [ys[i][j] + errors[i][j] for j in range(len(ys[i]))], alpha=error_alpha, color=colors[i])

    # Set plot title and axis labels
    if title is not None:
        ax.set_title(title, fontsize=title_size)
    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=x_label_size)
    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=y_label_size)
    
    # Handle custom x-axis ticks and secondary x-axis
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
        if x_ticks2 is not None:
            ax2 = ax.twiny()
            ax2.set_xticks(x_ticks2)
            ax2.tick_params(axis='both', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width)
    
    # Handle custom y-axis ticks and secondary y-axis
    if y_ticks is not None:
        ax.set_ylabel(labels[0] + ' Counts', color=colors[0])
        ax.set_yticks(np.linspace(min(ys[0]), 0.85 * y_max, len(y_ticks)))
        ax.set_yticklabels([str(_) for _ in y_ticks])
        ax.tick_params(axis='y', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width,
                       labelcolor=colors[0])
        if y_ticks2 is not None:
            ax3 = ax.twinx()
            ax3.set_ylabel(labels[1] + ' Counts', color=colors[1], fontsize=y_label_size)
            ax3.set_yticks(np.linspace(0, 0.85, len(y_ticks2)))
            ax3.set_yticklabels([str(_) for _ in y_ticks2])
            ax3.tick_params(axis='both', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width,
                            labelcolor=colors[1])
    
    # Set tick parameters for main axes
    ax.tick_params(axis='both', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width)
    
    # Handle legend or colorbar
    if labels is not None:
        if colorbar is not None:
            # Create colorbar with legend title
            cbar = plt.colorbar(colorbar, ax=ax)
            cbar.set_label(legend_title, fontdict=dict(size=legend_title_size))
            cbar.ax.tick_params(labelsize=legend_label_size, size=12, width=2, length=12)
        else:
            # Create legend with configurable orientation
            ncol = 1
            if legend_orientation.lower() == 'horizontal':
                ncol = len(labels)
            legend = ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1), prop={'size': legend_label_size}, ncol=ncol,
                               shadow=True)
            if legend_title is not None:
                legend.set_title(legend_title)
                legend.get_title().set_fontsize(str(legend_title_size))

    # Adjust layout to accommodate legend
    plt.subplots_adjust(right=0.8)
    if tight_layout:
        plt.tight_layout()
    
    # Display the plot if requested
    if Show:
        plt.show()