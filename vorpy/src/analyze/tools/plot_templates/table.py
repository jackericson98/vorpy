import matplotlib.pyplot as plt
import matplotlib as mpl


def color_cells(rows, header_color, color_cols=None, color_map='RdYlGn'):
    # Create the colormap
    my_cmap = mpl.colormaps[color_map]

    # Color the cells:
    header = [[header_color for _ in range(len(rows[0]))]]
    if color_cols is None:
        return header + [[(1, 1, 1) for _ in range(len(rows[0]))] for __ in range(len(rows))]
    else:
        # Get the col maxes
        cols = [[] for _ in range(len(rows[0]))]
        for row in rows:
            for i, col in enumerate(row):
                if type(col) == float or type(col) == int and color_cols[i]:
                    cols[i].append(col)
                else:
                    cols[i].append(0)
        col_maxes = [max(_) for _ in cols]
        colors = header
        for row in rows:
            colors.append([])
            for i, col in enumerate(row):
                if color_cols[i]:
                    print(col)
                    colors[-1].append(my_cmap(col/col_maxes[i]))
                else:
                    colors[-1].append((1, 1, 1))
        return colors


def table(rows, column_names, color_cols=None, Show=False, header_color=(0.9, 0.9, 0.9), transpose=False):

    # Create figure and axis
    fig, ax = plt.subplots()

    # If the user wants to transpose the table
    if transpose:
        new_table = []
        for i in range(len(column_names)):
            row = [column_names[i]]
            for j in range(len(rows)):
                row.append(rows[j][i])
            new_table.append(row)
    else:
        new_table = [column_names] + rows
    graph_labels = {'EDTA_Mg': 'EDTA', 'cambrin': 'Cambrin', 'hairpin': 'Hairpin', 'p53tet': 'p53tet',
                     'streptavidin': 'STVDN', '3zp8_hammerhead': 'H-head', 'NCP': 'NCP'}
    new_table = [graph_labels[_] for _ in new_table[0]] + new_table[1:]

    # Get the colors for the cells
    cell_colors = color_cells(new_table[1:], header_color, color_cols)

    # Create table
    my_table = ax.table(cellText=new_table, loc='center', cellLoc='center', colLabels=None,
                        cellColours=cell_colors)

    # Hide axes
    ax.axis('off')

    # Set font size and style
    my_table.auto_set_font_size(False)
    my_table.set_fontsize(12)

    # Auto-size columns to fit the largest member
    for col in range(len(column_names)):
        my_table.auto_set_column_width(col)

    # Adjust layout
    my_table.scale(1, 1.5)  # Adjust scale for better readability

    # Show plot
    if Show:
        plt.show()
