import matplotlib.pyplot as plt
import numpy as np


def inline_plot(data, labels=None, axis_label='', orientation='vertical', log_scale=False, types=None, Show=False,
                legend_title=None, title=''):

    # Get the unique labels
    if types is not None:
        type_vals = {}
        for i in range(len(data)):
            if types[i] in type_vals:
                type_vals[types[i]]['data'].append(data[i])
                if labels is not None:
                    type_vals[types[i]]['labels'].append(labels[i])
            else:
                type_vals[types[i]] = {'data': [data[i]]}
                if labels is not None:
                    type_vals[types[i]]['labels'] = [labels[i]]
    else:
        type_vals = {'Data': data}

    # Create a list of colors to reference
    colors = ['skyblue', 'orange', 'green', 'slategray', 'coral', 'periwinkle',
              'turquoise']

    # Create a list of markers to reference
    markers = ['o', 's', 'x', '^', '-']

    # Plot the data
    for i, _ in enumerate(type_vals):
        if orientation == 'horizontal':
            plt.scatter(type_vals[_]['data'], [0.5 for _ in range(len(type_vals[_]['data']))], c=colors[i], label=_,
                        marker=markers[i])
            # Plot the labels
            if labels is not None:
                for j, label in enumerate(type_vals[_]['labels']):
                    plt.text(type_vals[_]['data'][j], 0.55, label, rotation=45)
        else:
            plt.scatter([0.5 for _ in range(len(type_vals[_]['data']))], type_vals[_]['data'], c=colors[i], label=_,
                        marker=markers[i])
            # Plot the labels
            if labels is not None:
                for j, label in enumerate(type_vals[_]['labels']):
                    plt.text(0.55, type_vals[_]['data'][j] - 0.1, label, va='center')

    # If horizontal hide the y axis
    if orientation == 'horizontal':
        # Hide the y-axis
        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(True)

        # Hide the y ticks
        plt.yticks([])

        # Set x-axis to log scale
        if log_scale:
            plt.xscale('log')

        # Set the range for the y axis
        plt.ylim(0, 2)

        # Plot the x axis label
        plt.xlabel(axis_label, fontdict=dict(size=15))

    else:
        # Hide the x-axis
        plt.gca().spines['left'].set_visible(True)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)

        # Hide the x ticks
        plt.xticks([])

        # Set y-axis to log scale
        if log_scale:
            plt.yscale('log')

        # Set the range for the x axis
        plt.xlim(0, 2)

        # Plot the y axis label
        plt.ylabel(axis_label, fontdict=dict(size=15))

    # Insert the legend
    if legend_title is not None:
        plt.legend(title=legend_title)
    elif len([_ for _ in type_vals]) > 1:
        plt.legend()

    # Plot the title
    plt.title(title, fontdict=dict(size=20))

    # Show the plot
    if Show:
        plt.show()


# Test
# inline_plot([1, 10, 100, 1000, 10000, 100000], types=['gay', 'gay', 'straight', 'bi', 'trans', 'bi'],
#             log_scale=True, Show=True, orientation='horizontal',
#             labels=['aaaaaaaaa', 'bbbbbbbb', 'cccccccc', 'dddddddd', 'eeeeeeeeee', 'ffffffff'], title='num gays in sf', axis_label='num gays')
