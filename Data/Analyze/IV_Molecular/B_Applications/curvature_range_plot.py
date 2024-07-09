from System.sys_funcs.calcs.surf import calc_surf_point_curv, calc_surf_func
from Data.Analyze.tools.plot_templates.line import line_plot
from Visualize.mpl_visualize import plot_atoms
from System.sys_funcs.output.atoms import make_pdb_line


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
        print('\n')
        print(make_pdb_line(tfact=1.0))
        print(make_pdb_line(x=large_atom_loc, tfact=large_atom_rad))
        # Calculate the curvature
        if large_atom_rad == 1:
            curvature = 0
        else:
            curvature = calc_surf_point_curv(func, [center_point, 0, 0])
        # Record the data
        data[distance][large_atom_rad] = curvature
        count += 1
        print('\rProgress: {:.2f} %'.format(100 * round(count / 200)), end='')

for _ in data:
    for __ in data[_]:
        print()
# Plot the data
distance = [_ for _ in data]
radii = [_ for _ in data[-1.0]]

ys = [[data[_][__] for __ in data[_]] for _ in data]
line_plot([[__ for __ in data[_]] for _ in data], ys, x_label='Ratio of Radii', y_label='Gaussian Curvature', labels=distance, legend_label_size=15, legend_title='Surface \nDistance', legend_title_size=15, title='Curvature by Distance and Radii', tick_val_size=15)

