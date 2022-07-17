import matplotlib.pyplot as plt
from Cells.build_mesh import *
from System.system import Atom


# Set up plot function. Used to set the parameters for the plot
def setup_plot(fig=None, ax=None, dfo=None, grid=False, alpha=None):
    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
            Show = True  # If no outside figure is specified, then the figure needs to be shown from within
        ax = fig.add_subplot(projection="3d")
    ax.set_facecolor('k')
    ax.axis('auto')
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
    return fig, ax, alpha


# Plot spheres function. Plots the spheres specified
def plot_atoms(atoms, colors=None, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha)

    # Get the atoms colors
    if colors is None:
        atom_colors = {1.2: 'w', 1.52: 'r', 2.29: 'g', 1.55: 'b', 1.7: 'grey', 1.8: 'y'}
        colors = []
        for atom in atoms:
            try:
                colors.append(atom_colors[atom.rad])
            except KeyError:
                colors.append('pink')
    else:
        colors = colors + ['pink'for i in range(abs(len(atoms) - len(colors)))]
    # If the number of atoms to plot is more than 80, then plot them as points rather than spheres.
    if len(atoms) > 80:
        for i in range(len(atoms)):
            ax.scatter(atoms[i].loc[0], atoms[i].loc[1], atoms[i].loc[2], s=20, c=colors[i], alpha=alpha)
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
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha)
    # Set up the colors
    if colors is None:
        colors = ['y' for i in range(len(surfs))]
    elif len(colors) < len(surfs):
        colors = colors + ['y' for i in range(len(surfs) - len(colors))]
    # Plot the surfaces
    for i in range(len(surfs)):
        x, y, z = [], [], []
        # Add the surface points
        for point in surfs[i].points + surfs[i].edge_points:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        # ax.plot_trisurf(x, y, z, alpha=alpha)
        ax.scatter(x, y, z, s=[0.1 for j in range(len(x))], alpha=alpha, c=[colors[i] for k in range(len(x))])
    # Show the figure
    if Show:
        plt.show()


# Plot edges function. Plots the edges given as lines
def plot_edges(edges, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha)
    # Set the color if it is not indicated already
    if colors is None:
        colors = ['grey' for i in range(len(edges))]
    elif len(colors) < len(edges):
        colors = colors + ['grey' for i in range(len(edges) - len(colors))]
    # Plot the edges
    for i in range(len(edges)):
        xs, ys, zs = [edges[i].verts[0].loc[0]], [edges[i].verts[0].loc[1]], [edges[i].verts[0].loc[2]]
        for point in edges[i].points:
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])
        xs.append(edges[i].verts[1].loc[0])
        ys.append(edges[i].verts[1].loc[1])
        zs.append(edges[i].verts[1].loc[2])
        # Plot the points
        ax.plot(xs, ys, zs, c=colors[i])
    # Show the figure
    if Show:
        plt.show()


# Plot vertices function. Plots the vertices of a network.
def plot_verts(verts, plot_spheres=False, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha)
    # Default color is red
    if colors is None:
        colors = ['r' for i in range(len(verts))]
    # Plot each vertex
    for i in range(len(verts)):
        # Plot the point
        ax.scatter(verts[i].loc[0], verts[i].loc[1], verts[i].loc[2], c=colors[i])
    # Plot the inscribed spheres
    if plot_spheres:
        spheres = []
        for i in range(len(verts)):
            spheres.append(Atom(verts[i].loc, verts[i].rad))
        plot_atoms(spheres, fig=fig, ax=ax)
    # Show if the plot needs to be shown
    if Show:
        plt.show()
