import matplotlib.pyplot as plt
import numpy as np


# Set up plot function. Used to set the parameters for the plot
def setup_plot(fig=None, ax=None, dfo=None, grid=False, alpha=None, bg_color=None):
    # Create a new subplot if one isn't specified
    if ax is None:
        # Create new figure if one isn't specified
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
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
        if bg_color:
            ax.set_facecolor(bg_color)
        else:
            ax.set_facecolor('k')
    # Set alpha
    if alpha is None:
        alpha = 0.5
    return fig, ax, alpha


# Plot spheres function. Plots the spheres specified
def plot_atoms(atoms=None, atom_list=None, colors=None, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None,
               bg_color=None, res=4):

    # Give an option for not atom objects (lists) to be plotted
    locs, rads = [], []
    if atom_list is None:
        for i in range(len(atoms)):
            locs.append(atoms[i].loc)
            rads.append(atoms[i].rad)
    else:
        for i in range(len(atom_list)):
            locs.append(atom_list[i][0])
            rads.append(atom_list[i][1])
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha, bg_color)
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
        colors = colors + ['pink'for _ in range(abs(len(locs) - len(colors)))]
    # If the number of atoms to plot is more than 80, then plot them as points rather than spheres.
    if len(locs) > 80:
        for i in range(len(locs)):
            ax.scatter(locs[i][0], locs[i][1], locs[i][2], s=20, c=colors[i], alpha=alpha)
    # Plot the spheres as wireframes
    else:
        # Set the resolution of the spheres
        f = 5 - len(locs) // 20
        # Find u, v values that span phi and theta
        u, v = np.mgrid[0:2 * np.pi:f*res*2j, 0:np.pi:f*res*1j]
        # Plot each sphere
        for i in range(len(locs)):
            # Get x, y, z data for the wireframe
            x = rads[i] * np.cos(u) * np.sin(v) + locs[i][0]
            y = rads[i] * np.sin(u) * np.sin(v) + locs[i][1]
            z = rads[i] * np.cos(v) + locs[i][2]
            # Plot the sphere
            ax.plot_wireframe(x, y, z, color=colors[i], alpha=alpha)
    # Show the figure if need be
    if Show:
        plt.show()


# Plot vertices function. Plots the vertices of a network.
def plot_verts(verts, spheres=False, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None,
               bg_color=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha, bg_color)
    # Default color is red
    if colors is None:
        colors = ['r' for _ in range(len(verts))]
    # Plot each vertex
    for i in range(len(verts)):
        # Plot the point
        ax.scatter(verts[i].loc[0], verts[i].loc[1], verts[i].loc[2], c=colors[i])
    # Plot the inscribed spheres
    if spheres:
        spheres = []
        for i in range(len(verts)):
            spheres.append([verts[i].loc, verts[i].rad])
        plot_atoms(atom_list=spheres, fig=fig, ax=ax, colors=['grey'], alpha=1)
    # Show if the plot needs to be shown
    if Show:
        plt.show()


# Plot edges function. Plots the edges given as lines
def plot_edges(edges, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None, bg_color=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha, bg_color)
    # Set the color if it is not indicated already
    if colors is None:
        colors = ['grey' for _ in range(len(edges))]
    elif len(colors) < len(edges):
        colors = colors + ['grey' for _ in range(len(edges) - len(colors))]
    # Plot the edges
    for i in range(len(edges)):
        xs, ys, zs = [edges[i].pv0[0]], [edges[i].pv0[1]], [edges[i].pv0[2]]
        for point in edges[i].points:
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])
        xs.append(edges[i].pv1[0])
        ys.append(edges[i].pv1[1])
        zs.append(edges[i].pv1[2])

        # Plot the points
        ax.plot(xs, ys, zs, c=colors[i])
    # Show the figure
    if Show:
        plt.show()


# Plot surfaces function. Plots the surfaces given
def plot_surfs(surfs, simps=False, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None,
               bg_color=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha, bg_color)
    # Set up the colors
    if colors is None:
        colors = ['w' for _ in range(len(surfs))]
    elif len(colors) < len(surfs):
        colors = colors + ['w' for _ in range(len(surfs) - len(colors))]
    # Plot the surfaces
    for i in range(len(surfs)):
        x, y, z = [], [], []
        for point in surfs[i].points:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        # If simplices are requested get them or make them
        if simps:
            # Plot the simps using matplotlib tri_surf
            ax.plot_trisurf(x, y, z, triangles=surfs[i].tris, alpha=alpha, color=colors[i])
        # Otherwise, plot the points
        else:
            ax.scatter(x, y, z, s=[0.1 for _ in range(len(x))], alpha=alpha, c=[colors[i] for _ in range(len(x))])
    # Show the figure
    if Show:
        plt.show()


# Plot simplices function
def plot_simps(surf, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None, bg_color=None):
    # Set up the plot
    fig, ax, alpha = setup_plot(fig, ax, dfo, grid, alpha, bg_color)
    # Go through each triangle in the surfaces list of simplices
    for simp in surf.simps.triangles:
        p0, p1, p2 = surf.points[simp[0]], surf.points[simp[1]], surf.points[simp[2]]
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], c='w', linewidth=.1)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], c='w', linewidth=.1)
        ax.plot([p2[0], p0[0]], [p2[1], p0[1]], [p2[2], p0[2]], c='w', linewidth=.1)
    # Show the figure
    if Show:
        plt.show()
