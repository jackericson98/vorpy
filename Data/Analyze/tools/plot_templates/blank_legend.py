import matplotlib.pyplot as plt


def legend_plot(labels, markers=None, label_colors=None, data_type='Line', legend_orientation='Vertical',
                legend_title=None, legend_title_size=20, legend_entry_size=20, legend_loc='upper left', marker_size=None):

    fig, ax = plt.subplots()

    if data_type == 'scatter':
        if markers is None:
            markers = ['.' for _ in labels]
        if label_colors is None:
            label_colors = ['k' for _ in labels]
        if marker_size is None:
            marker_size = [100 for _ in labels]
        for i, label in enumerate(labels):
            ax.scatter([0], [0], label=label, color=label_colors[i], marker=markers[i], s=marker_size[i])

    if data_type == 'line':
        if label_colors is None:
            label_colors = ['k' for _ in labels]


    # Add the legend
    leg_col = 1
    if legend_orientation == 'Horizontal':
        leg_col = len(xs)

    if legend_title is not None:
        legend = ax.legend(title=legend_title, loc=legend_loc, shadow=True,
                           ncol=leg_col,
                           prop={'size': legend_entry_size})
        legend.get_title().set_fontsize(str(legend_title_size))
    else:
        legend = ax.legend(loc=legend_loc, shadow=True, ncol=leg_col,
                           prop={'size': legend_entry_size})

        legend.get_title().set_fontsize(str(legend_title_size))
    plt.show()

if __name__ == '__main__':
    # color_dict = {'NP': 'green', '+': 'red', 'P': 'cyan', '-': 'k'}
    # legend_plot(['Polar', 'NonPolar', 'Positively Charged', 'Negatively Charged'], label_colors=['cyan', 'green', 'red', 'k'],
    #             markers=['o' for _ in range(4)], legend_title='Residue Type', data_type='scatter',
    #             marker_size=[200 for _ in range(4)])
    legend_plot(['Carbon', 'Oxygen', 'Hydrogen', 'Nitrogen', 'Phosphorous', 'Sulfur', 'Selenium'], label_colors=['grey', 'r', 'pink', 'b', 'darkorange', 'y', 'sandybrown'],
                markers=['o' for _ in range(7)], legend_title='Element Type', data_type='scatter',
                marker_size=[200 for _ in range(7)])