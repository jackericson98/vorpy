import matplotlib.pyplot as plt
from build_mesh import *
from objects import Atom


# Set up plot function. Used to set the parameters for the plot
def setup_plot(num_col, colors=None, fig=None, ax=None, dfo=None, grid=False, alpha=None):
    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
            Show = True  # If no outside figure is specified, then the figure needs to be shown from within
        ax = fig.add_subplot(projection="3d")

    # Set the colors of the spheres. Defaults to blue with a white base sphere
    if colors is None:
        colors = ['b' for _ in range(num_col)]
    # If not all colors are specified make the rest blue
    if len(colors) < num_col:
        colors = colors + ['b' for _ in range(num_col - len(colors))]
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
    # Set alpha
    if alpha is None:
        alpha = 1
    return fig, ax, colors, alpha


# Plot spheres function. Plots the spheres specified
def plot_atoms(atoms, colors=None, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None):
    # Set up the plot
    fig, ax, colors, alpha = setup_plot(len(atoms), colors, fig, ax, dfo, grid, alpha)

    # If the number of atoms to plot is more than 80, then plot them as points rather than spheres.
    if len(atoms) > 80:
        for sphere in atoms:
            ax.scatter(sphere.loc[0], sphere.loc[1], sphere.loc[2])
    # Plot the spheres as wireframes
    else:
        # Set the resolution of the spheres
        res = 5 - len(atoms) // 20
        # Find u, v values that span phi and theta
        u, v = np.mgrid[0:2 * np.pi:res*8j, 0:np.pi:res*4j]
        # Plot each sphere
        for i in range(len(atoms)):
            # Get x, y, z data for the wireframe
            x = atoms[i].rad * np.cos(u) * np.sin(v) + atoms[i].loc[0]
            y = atoms[i].rad * np.sin(u) * np.sin(v) + atoms[i].loc[1]
            z = atoms[i].rad * np.cos(v) + atoms[i].loc[2]
            # Plot the sphere
            ax.plot_wireframe(x, y, z, color=colors[i], alpha=alpha)

    # Show the figure if need be
    if Show:
        plt.show()


# Plot surfaces function. Plots the surfaces given
def plot_surfs(surfs, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None):
    # Set up the plot
    fig, ax, colors, alpha = setup_plot(len(surfs), colors, fig, ax, dfo, grid, alpha)

    # Plot the surfaces
    for surf in surfs:
        x, y, z = [], [], []
        for point in surf.points:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        ## ax.plot_trisurf(x, y, z)
        ax.scatter(x, y, z, s=[0.1 for i in range(len(x))], alpha=0.1)

    # Show the figure
    if Show:
        plt.show()


# Plot edges function. Plots the edges given as lines
def plot_edges(edges, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None):
    # Set up the plot
    fig, ax, colors, alpha = setup_plot(len(edges), colors, fig, ax, dfo, grid, alpha)

    # Plot the edges
    for edge in edges:
        xs, ys, zs = [], [], []
        for point in edge.points:
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])
        # Plot the points
        ax.plot(xs, ys, zs)
    # Show the figure
    if Show:
        plt.show()


# Plot vertices function. Plots the vertices of a network.
def plot_verts(verts, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None):
    # Set up the plot
    fig, ax, colors, alpha = setup_plot(len(verts), colors, fig, ax, dfo, grid, alpha)

    # Plot each vertex
    for i in range(len(verts)):
        # Plot the point
        ax.scatter(verts[i].loc[0], verts[i].loc[1], verts[i].loc[2], c=colors[i])

    # Show if the plot needs to be shown
    if Show:
        plt.show()
