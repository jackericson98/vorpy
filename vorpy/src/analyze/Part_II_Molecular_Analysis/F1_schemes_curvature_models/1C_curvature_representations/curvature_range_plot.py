import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.calculations.surf import calc_surf_func
from vorpy.src.calculations.curvature import mean_curvature
from vorpy.src.analyze.tools.plot_templates.line import line_plot
from vorpy.src.output import make_pdb_line


# Generate a range of ratios and distances
data = {}
count = 0
diff = 0
for j in range(10):
    # We want a range fromm -0.5 to 9.5
    diff += 0.05 * j
    distance = round(-1.0 + 0.1 * j + diff, 3)
    data[distance] = {}
    for i in range(20):
        # Calculate the large atoms radius
        large_atom_rad = (i + 2) * 0.5

        # This needs to account for the radius of the large atom
        large_atom_loc = 1 + large_atom_rad + distance
        # Find the ceter point of the surface
        center_point = 1 + 0.5 * distance
        # plot_atoms([[0, 0, 0], [large_atom_loc, 0, 0]], [1, large_atom_rad], Show=True)


        # Get the surf func
        func = calc_surf_func([0, 0, 0], 1, [large_atom_loc, 0, 0], large_atom_rad)

        # Calculate the curvature
        if large_atom_rad == 1:
            curvature = 0
        else:
            curvature = mean_curvature(func, [center_point, 0, 0])
        # Record the data
        data[distance][large_atom_rad] = curvature
        count += 1
        print('\rProgress: {:.2f} %'.format(100 * round(count / 200)), end='')

# Plot the data
distance = [_ for _ in data]
radii = [_ for _ in data[-1.0]]

cmap = plt.cm.rainbow
norm = Normalize(vmin=min(distance), vmax=max(distance))
sm = ScalarMappable(norm=norm, cmap=cmap)
line_colors = [cmap(norm(_)) for _ in distance]



ys = [[data[_][__] for __ in data[_]] for _ in data]
line_plot([[__ for __ in data[_]] for _ in data], ys, colors=line_colors, Show=True, x_label='Radii Ratio', y_label='Mean Curvature',
          labels=distance, legend_label_size=30, legend_title='Surface Distance', legend_title_size=30,
          title='Curvature by Distance and Radii', tick_val_size=30, x_label_size=30, y_label_size=30, linewidth=3, colorbar=sm, figsize=(8, 6))

