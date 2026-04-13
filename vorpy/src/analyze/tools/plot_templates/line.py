import matplotlib.pyplot as plt
import numpy as np



def line_plot(xs, ys, errors=None, labels=None, error_alpha=0.2, title=None, x_label=None, y_label=None,
              legend_title=None, Show=True, title_size=25, x_label_size=20, y_label_size=20, legend_title_size=20,
              legend_label_size=20, tick_width=2, tick_length=12, legend_orientation='vertical', colors=None,
              tick_val_size=15, linewidth=1, x_ticks=None, x_ticks2=None, y_ticks=None, y_ticks2=None, colorbar=None,
              tight_layout=True, figsize=None, ylim=None, xlim=None, axis_line_thickness=1.5, legend_loc='upper right',
              legend_bbox_to_anchor=(1.2, 1), alpha=1, y_label2=None, ylim2=None):
    """
    xs, ys: lists of equal-length sequences
      - xs[i] are x-values for line i
      - ys[i] are y-values for line i

    Left axis controls:
      - ylim, y_ticks, y_label

    Right axis controls (optional; created if any of these provided):
      - ylim2, y_ticks2, y_label2
    """

    if colors is None:
        # Fixed invalid 'o' color name -> 'orange'
        colors = ['r', 'g', 'b', 'y', 'orange', 'pink', 'purple', 'lavender']

    fig, ax = plt.subplots(figsize=figsize)

    ax_right = None
    ax_top = None

    # ---- plot main lines on left axis ----
    for i in range(len(xs)):
        x_i = np.asarray(xs[i], dtype=float)
        y_i = np.asarray(ys[i], dtype=float)

        if labels is not None:
            ax.plot(x_i, y_i, label=labels[i], linewidth=linewidth, c=colors[i], alpha=alpha)
        else:
            ax.plot(x_i, y_i, linewidth=linewidth, c=colors[i], alpha=alpha)

        if errors is not None:
            e_i = np.asarray(errors[i], dtype=float)
            ax.fill_between(
                x_i,
                y_i - e_i,
                y_i + e_i,
                alpha=error_alpha,
                color=colors[i]
            )

    # ---- titles / labels ----
    if title is not None:
        ax.set_title(title, fontsize=title_size)

    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=x_label_size)

    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=y_label_size)

    # ---- left axis limits / ticks ----

    ax.set_autoscale_on(False)

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    if x_ticks is not None:
        ax.set_xticks(x_ticks)

    if y_ticks is not None:
        ax.set_yticks(y_ticks)

    # ---- optional top x-axis ----
    if x_ticks2 is not None:
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        ax_top.set_xticks(x_ticks2)
        ax_top.tick_params(axis='both', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width)

    # ---- optional right y-axis (INDEPENDENT) ----
    if (y_ticks2 is not None) or (ylim2 is not None) or (y_label2 is not None):
        ax_right = ax.twinx()

        if y_label2 is not None:
            ax_right.set_ylabel(y_label2, fontsize=y_label_size)

        if ylim2 is not None:
            ax_right.set_ylim(ylim2)

        if y_ticks2 is not None:
            ax_right.set_yticks(y_ticks2)

        ax_right.tick_params(axis='y', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width)

    # ---- styling ----
    ax.tick_params(axis='both', which='major', labelsize=tick_val_size, length=tick_length, width=tick_width)

    for spine in ax.spines.values():
        spine.set_linewidth(axis_line_thickness)

    if ax_top is not None:
        for spine in ax_top.spines.values():
            spine.set_linewidth(axis_line_thickness)

    if ax_right is not None:
        for spine in ax_right.spines.values():
            spine.set_linewidth(axis_line_thickness)

    # ---- colorbar / legend ----
    if colorbar is not None:
        cbar = fig.colorbar(colorbar, ax=ax)
        if legend_title is not None:
            cbar.set_label(legend_title, fontdict=dict(size=legend_title_size))
        cbar.ax.tick_params(labelsize=legend_label_size, size=12, width=2, length=12)

    if labels is not None and colorbar is None:
        ncol = 1
        if str(legend_orientation).lower() == 'horizontal':
            ncol = len(labels)

        legend = ax.legend(
            loc=legend_loc,
            bbox_to_anchor=legend_bbox_to_anchor,
            prop={'size': legend_label_size},
            ncol=ncol,
            shadow=True
        )

        if legend_title is not None:
            legend.set_title(legend_title)
            legend.get_title().set_fontsize(str(legend_title_size))

    # ---- layout ----
    # Only reserve space if they intentionally pushed legend outside
    if labels is not None and colorbar is None and legend_bbox_to_anchor is not None and legend_bbox_to_anchor[0] > 1.0:
        fig.subplots_adjust(right=0.8)

    if tight_layout:
        fig.tight_layout()

    if Show:
        plt.show()
