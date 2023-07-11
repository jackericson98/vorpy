import matplotlib.pyplot as plt
from System.sys_funcs.calcs.calcs import calc_circ
import numpy as np


# Set up plot function. Used to set the parameters for the plot
def setup_plot(fig=None, ax=None, dfo=None, grid=False, bg_color=None):
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
    return fig, ax


# Plot spheres function. Plots the spheres specified
def plot_atoms(alocs, arads, colors=None, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None,
               bg_color=None, res=4):

    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Get the atoms colors
    if colors is None:
        colors = ['pink'for _ in range(abs(len(alocs)))]
    # If the number of atoms to plot is more than 80, then plot them as points rather than spheres.
    if len(alocs) > 80:
        for i in range(len(alocs)):
            ax.scatter(alocs[i][0], alocs[i][1], alocs[i][2], s=20, c=colors[i], alpha=alpha)
    # Plot the spheres as wireframes
    else:
        # Set the resolution of the spheres
        f = max(3 - len(alocs) // 40, 1)
        # Find u, v values that span phi and theta
        u, v = np.mgrid[0:2 * np.pi:f*res*2j, 0:np.pi:f*res*1j]
        # Plot each sphere
        for i in range(len(alocs)):
            # Get x, y, z data for the wireframe
            x = arads[i] * np.cos(u) * np.sin(v) + alocs[i][0]
            y = arads[i] * np.sin(u) * np.sin(v) + alocs[i][1]
            z = arads[i] * np.cos(v) + alocs[i][2]
            # Plot the sphere
            ax.plot_wireframe(x, y, z, color=colors[i], alpha=alpha)
    # Show the figure if need be
    if Show:
        plt.show()


# Plot vertices function. Plots the vertices of a network.
def plot_verts(vlocs, vrads, spheres=False, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None,
               bg_color=None):
    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Default color is red
    if colors is None:
        colors = ['b' for _ in range(len(vlocs))]
    # Plot each vertex
    for i in range(len(vlocs)):
        # Plot the point
        ax.scatter(vlocs[i][0], vlocs[i][1], vlocs[i][2], c=colors[i])
    # Plot the inscribed spheres
    if spheres:
        plot_atoms(alocs=vlocs, arads=vrads, fig=fig, ax=ax, colors=colors, alpha=alpha)
    # Show if the plot needs to be shown
    if Show:
        plt.show()


# Plot edges function. Plots the edges given as lines
def plot_edges(epnts, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None, bg_color=None, center=None):
    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Set the color if it is not indicated already
    if colors is None:
        colors = ['grey' for _ in range(len(epnts))]
    elif len(colors) < len(epnts):
        colors = colors + ['grey' for _ in range(len(epnts) - len(colors))]
    # Plot the edges
    for i in range(len(epnts)):
        xs, ys, zs = [], [], []
        for point in epnts[i]:
            xs.append(point[0])
            ys.append(point[1])
            zs.append(point[2])

        # Plot the points
        ax.plot(xs, ys, zs, c=colors[i], linewidth=5)
    if center is not None:
        ax.plot(center[0], center[1], center[2], c='r', marker='x', markersize=20)
    # Show the figure
    if Show:
        plt.show()


# Plot surfaces function. Plots the surfaces given
def plot_surfs(spnts, stris, simps=True, fig=None, ax=None, Show=False, dfo=None, grid=False, colors=None, alpha=None,
               bg_color=None):
    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Set up the colors
    if colors is None:
        colors = ['w' for _ in range(len(spnts))]
    elif len(colors) < len(spnts):
        colors = colors + ['w' for _ in range(len(spnts) - len(colors))]
    # Plot the surfaces
    for i in range(len(spnts)):
        x, y, z = [], [], []
        for point in spnts[i]:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        # If simplices are requested get them or make them
        if simps:
            # Plot the simps using matplotlib tri_surf
            ax.plot_trisurf(x, y, z, triangles=stris[i], alpha=alpha, color=colors[i])
        # Otherwise, plot the points
        else:
            ax.scatter(x, y, z, s=[0.1 for _ in range(len(x))], alpha=alpha, c=[colors[i] for _ in range(len(x))])
    # Show the figure
    if Show:
        plt.show()


# Plot simplices function.
def plot_simps(spnts, stris, fig=None, ax=None, Show=False, dfo=None, grid=False, alpha=None, bg_color=None):
    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Go through each triangle in the surfaces list of simplices
    for tri in stris:
        p0, p1, p2 = spnts[tri[0]], spnts[tri[1]], spnts[tri[2]]
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], c='w', linewidth=.1)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], c='w', linewidth=.1)
        ax.plot([p2[0], p0[0]], [p2[1], p0[1]], [p2[2], p0[2]], c='w', linewidth=.1)
    # Show the figure
    if Show:
        plt.show()


# Plot network function. Plots the network items
def plot_net(net, group=None, plot_all=False, atoms=False, verts=False, edges=False, surfs=False, fig=None, ax=None, grid=False,
             bg_color='white', dfo=None, Show=True, a_alpha=1, v_alpha=1, e_alpha=1, s_alpha=1):
    # Check for a figure or an ax
    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
    # Set up the plot
    fig, ax = setup_plot(fig, ax, dfo, grid, bg_color)
    # Check to see if the group is availible
    if group is not None:
        my_atoms = group.atoms
    else:
        my_atoms = net.atoms['num']
    # Atoms
    if atoms or plot_all:
        alocs = [net.atoms['loc'][i] for i in my_atoms]
        arads = [net.atoms['rad'][i] for i in my_atoms]
        plot_atoms(alocs=alocs, arads=arads, fig=fig, ax=ax, alpha=a_alpha)
    # Vertices
    if verts or plot_all:
        plot_verts(net.verts['vloc'], net.verts['vrad'], fig=fig, ax=ax, colors=['r' for _ in range(len(net.verts))], alpha=v_alpha, spheres=True)
    # Edges
    if edges or plot_all:
        plot_edges(net.edges['points'], fig=fig, ax=ax, colors=['w' for _ in range(len(net.edges))], alpha=e_alpha)
    # Surfaces
    if surfs or plot_all:
        plot_surfs(net.surfs['points'], net.surfs['tris'], fig=fig, ax=ax, colors=[np.random.rand(3) for _ in range(len(net.surfs))], alpha=s_alpha)
    # Show the plot
    if Show:
        plt.show()
