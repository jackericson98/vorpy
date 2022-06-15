import matplotlib.pyplot as plt
from build_mesh import *


# Plot spheres function. Plots the spheres specified
def plot_atoms(atoms, colors=None, fig=None, ax=None, Show=False, dfo=None, grid=False):

    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
            Show = True  # If no outside figure is specified, then the figure needs to be shown from within
        ax = fig.add_subplot(projection="3d")

    # Set the colors of the spheres. Defaults to blue with a white base sphere
    if colors is None:
        colors = ['w'] + ['b' for _ in range(len(atoms) - 1)]
    # If not all colors are specified make the rest blue
    if len(colors) < len(atoms):
        colors = colors + ['b' for _ in range(len(atoms) - len(colors))]

    # If the number of spheres to plot is more than 80, then plot them as points rather than spheres. Marker size?
    if len(atoms) > 30:
        for sphere in atoms:
            ax.scatter(sphere.loc[0], sphere.loc[1], sphere.loc[2])
    else:
        res = 5 - len(atoms) // 20
        # Find u, v values that span phi and theta
        u, v = np.mgrid[0:2 * np.pi:res*8j, 0:np.pi:res*4j]
        rads, locs = [], []
        # Plot each sphere
        for i in range(len(atoms)):
            # Get radius and location values for the spheres
            rads.append(atoms[i].rad)
            locs.append(atoms[i].loc)
            # Get x, y, z data for the wireframe
            x = rads[i] * np.cos(u) * np.sin(v) + locs[i][0]
            y = rads[i] * np.sin(u) * np.sin(v) + locs[i][1]
            z = rads[i] * np.cos(v) + locs[i][2]

            # Plot the sphere
            ax.plot_wireframe(x, y, z, color=colors[i])

    # Set plot parameters
    if dfo is not None:
        ax.set_xlim(-dfo, dfo)
        ax.set_ylim(-dfo, dfo)
        ax.set_zlim(-dfo, dfo)
    # Set the grid if indicated
    if grid:
        ax.set_xlabel("X axis")
        ax.set_ylabel("Y axis")
        ax.set_zlabel("Z axis")
    else:
        ax.grid()
        ax.axis('off')
    # Show the figure if need be
    if Show:
        plt.show()


# Plot meshes function. Plots the meshes specified. If no meshes have been made user can specify spheres instead
def plot_surfs(atoms=None, surfs=None, fig=None, ax=None, Show=False, dfo=None, grid=False):
    plane = True
    if atoms is not None:
        surfs = []
        for pair in atoms:
            if pair[0].rad != pair[1].rad:
                plane = False
            surfs.append(make_meshes(pair[0], [pair[1]]))

    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
            Show = True  # If no outside figure is specified, then the figure needs to be shown from within
        ax = fig.add_subplot(projection="3d")

    for mesh in surfs:
        x, y, z = [], [], []
        for point in mesh:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])

        # Plot the mesh data values. If the mesh is a plane, plot it as scatter plot.
        if plane:
            ax.scatter(x, y, z, color='y', alpha=0.5)
        else:
            ax.plot_trisurf(x, y, z, color='y', alpha=0.5)
    if dfo is not None:
        ax.set_xlim(-dfo, dfo)
        ax.set_ylim(-dfo, dfo)
        ax.set_zlim(-dfo, dfo)

    # Set axes limits and labels
    # ax.set_xlabel('X Axis')
    # ax.set_ylabel('Y Axis')
    # ax.set_zlabel('Z Axis')
    ax.axis('off')
    ax.grid(grid)
    if Show:
        plt.show()


# Plot vertices function. Plots the vertices of a network.
def plot_verts(verts, fig=None, ax=None, Show=False, dfo=None, grid=False):
    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
            Show = True  # If no outside figure is specified, then the figure needs to be shown from within
        ax = fig.add_subplot(projection="3d")

    for vert in verts:
        ax.scatter(vert.loc[0], vert.loc[1], vert.loc[2])

    if Show:
        plt.show()
