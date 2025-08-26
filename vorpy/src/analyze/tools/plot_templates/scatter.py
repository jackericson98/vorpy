import matplotlib.pyplot as plt
import numpy as np


def scatter(xs, ys, labels=None, errors=None, title='', x_axis_title='', y_axis_title='', Show=False, save=None,
            legend_title=None,         # Title for the legend
            unit='',                   # Unit to display in axis labels or legend
            title_size=25,             # Font size for the plot title
            tick_width=2,              # Width of the axis ticks
            tick_length=12,            # Length of the axis ticks
            xlabel_size=20,            # Font size for x-axis label
            ylabel_size=20,            # Font size for y-axis label
            xtick_label_size=20,       # Font size for x-axis tick labels
            ytick_label_size=20,       # Font size for y-axis tick labels
            legend_entry_size=20,      # Font size for legend entries
            x_range=None,              # Range for x-axis (list or tuple: [min, max])
            y_range=None,              # Range for y-axis (list or tuple: [min, max])
            legend_orientation='Horizontal',  # Orientation of the legend
            alpha=None,                # Transparency for markers (float or list of floats)
            legend_title_size=20,      # Font size for legend title
            markers=None,              # Marker style(s) for scatter points (str or list)
            colors=None,               # Color(s) for scatter points (str or list)
            x_tick_labels=None,        # Custom x-tick labels
            x_tick_label_locs=None,    # Locations for custom x-tick labels
            y_tick_labels=None,        # Custom y-tick labels
            y_tick_label_locs=None,    # Locations for custom y-tick labels
            marker_size=100,           # Size(s) for markers (int, float, or list)
            xlabel_rotation=0,         # Rotation angle for x-axis labels
            xtick_anchor='center',     # Alignment for x-tick labels
            fig=None,                  # Matplotlib Figure object (optional)
            ax=None,                    # Matplotlib Axes object (optional)
            legend_labels=None,         # List of labels for the legend
            legend_handles=None,        # List of handles for the legend
            legend_loc='upper right',   # Location of the legend
            legend_bbox_to_anchor=(1.25, 0.97), # Bounding box for the legend
            axis_line_thickness=1
        ):

    # Set the colors
    if type(colors) is str:
        colors = [colors for _ in xs]
    elif colors is None:
        colors = ['skyblue', 'orange', 'red', 'blue', 'green', 'pink', 'y', 'purple']
        # color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # colors = [color_cycle[1 % len(color_cycle)] for i in range(len(xs))]

    # Set the markers
    if type(markers) is str:
        markers = [markers for _ in xs]
    elif markers is None:
        markers = ['.', ',', 'o', 's', 'v', '^', '<', '>', '1', '2', '3', '4', '*', '+', 'x', 'D', 'd', 'p', 'h', 'H', '|',
                   '_']
    if type(marker_size) is int or type(marker_size) is float:
        marker_size = [marker_size for _ in xs]
    elif marker_size is None:
        marker_size = [100 for _ in xs]

    if type(alpha) is int or type(alpha) is float:
        alpha = [alpha for _ in xs]
    elif alpha is None:
        alpha = [1 for _ in xs]

    # Create a single plot
    if fig is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    # Plot the actual data
    for i in range(len(xs)):
        my_label = None
        if labels is not None:
            my_label = labels[i]
        my_markers = markers[i]
        if type(markers[i]) is not list:
            my_markers = [markers[i] for _ in xs[i]]
        my_colors = colors[i]
        if type(colors[i]) is not list:
            my_colors = [colors[i] for _ in xs[i]]
        my_marker_sizes = marker_size[i]
        if type(marker_size[i]) is not list:
            my_marker_sizes = [marker_size[i] for _ in xs[i]]
        my_alphas = alpha[i]
        if type(my_alphas) is not list:
            my_alphas = [alpha[i] for _ in xs[i]]
        for j in range(len(xs[i])):

            ax.scatter(xs[i][j], ys[i][j], marker=my_markers[j], label=my_label, c=my_colors[j], alpha=alpha,
                       s=my_marker_sizes[j])
            my_label = None

    # Check how the data is set up and make sure it is a list of lists
    if type(xs[0]) is not list:
        xs = [xs]
        ys = [ys]

    # Get the total maximum for the list of lists
    err_max = 0
    if errors is not None:
        err_max = max([max(_) for _ in errors])
    # print(ys)
    ymax = max([max(_) for _ in ys]) + err_max
    ymin = min([min(_) for _ in ys]) - err_max

    # Get the number of bar groups to plot
    num_groups = range(len(ys[0]))

    # Plot the title, ylabel and xlabel
    # ax.title(title, fontdict=dict(size=title_size))
    plt.ylabel(y_axis_title, fontdict=dict(size=ylabel_size))
    plt.xlabel(x_axis_title, fontdict=dict(size=xlabel_size))

    # Label the bar groups
    plt.xticks(ticks=x_tick_label_locs, labels=x_tick_labels, font=dict(size=xtick_label_size), rotation=xlabel_rotation)
    plt.yticks(ticks=y_tick_label_locs, labels=y_tick_labels, font=dict(size=ytick_label_size))
    plt.tick_params(axis='both', width=tick_width, length=tick_length)

    ax.spines['top'].set_linewidth(axis_line_thickness)
    ax.spines['bottom'].set_linewidth(axis_line_thickness)
    ax.spines['left'].set_linewidth(axis_line_thickness)
    ax.spines['right'].set_linewidth(axis_line_thickness)

    # Add the legend
    leg_col = 1
    if legend_orientation == 'Horizontal':
        leg_col = len(xs)
    if legend_labels is not None:
        legend = ax.legend(handles=legend_handles, labels=legend_labels, title=legend_title, loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor, shadow=True, ncol=leg_col,
                    prop={'size': legend_entry_size})
        legend.get_title().set_fontsize(str(legend_title_size))
    elif legend_title is not None and legend_loc is not None:
        legend = ax.legend(title=legend_title, loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor, shadow=True, ncol=leg_col,
                  prop={'size': legend_entry_size})
        legend.get_title().set_fontsize(str(legend_title_size))
    elif type(labels) is list and len(labels) > 1:
        legend = ax.legend(loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor, shadow=True, ncol=leg_col,
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
        plt.ylim(ymin, ymax)
    # Add extra space below the title using 'pad'
    plt.title(title, fontdict=dict(size=title_size), pad=10)
    # Set the x limits
    if x_range is not None:
        plt.xlim(*x_range)

    # Set the figure size - use constrained_layout instead of tight_layout to handle legends outside plot area
    try:
        plt.tight_layout()
    except:
        # If tight_layout fails, use constrained_layout or adjust figure size
        fig.set_size_inches(10, 6)  # Increase figure size to accommodate legend
        plt.subplots_adjust(right=0.75)  # Adjust right margin to make room for legend

    if Show:
        plt.show()
