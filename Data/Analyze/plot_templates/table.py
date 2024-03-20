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


def table(rows, column_names, color_cols=None, Show=False, header_color=(0.9, 0.9, 0.9)):

    # Create figure and axis
    fig, ax = plt.subplots()

    # Get the colors for the cells
    cell_colors = color_cells(rows, header_color, color_cols)

    # Create table
    my_table = ax.table(cellText=[column_names] + rows, loc='center', cellLoc='center', colLabels=None, cellColours=cell_colors)

    # Hide axes
    ax.axis('off')

    # Set font size and style
    my_table.auto_set_font_size(False)
    my_table.set_fontsize(12)

    # Auto-size columns to fit the largest member
    for col in range(len(column_names)):
        col_width = max([len(str(row[col])) for row in rows])
        my_table.auto_set_column_width(col)

    # Adjust layout
    my_table.scale(1, 1.5)  # Adjust scale for better readability

    # Show plot
    if Show:
        plt.show()
