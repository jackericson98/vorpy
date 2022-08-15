from System.system import System
from Visualize.visualize import plot_atoms, plot_verts, plot_surfs, plot_edges
import matplotlib.pyplot as plt

# Create atom objects from sets of points
atoms = [[[0, 0, 0], 2]]

dist = 5
rad = 2

atoms += [[[dist, 0, 0], rad], [[-dist, 0, 0], rad], [[0, dist, 0], rad], [[0, -dist, 0], rad], [[0, 0, dist], rad],
          [[0, 0, -dist], rad]]

sys = System(atoms)

# Build the surfaces
sys.net.build(min_dist=.5)


# Set up the plot
fig = plt.figure()
ax = fig.add_subplot(projection="3d")


ax.set_title("Basic Cube Cell")
# Plot the elements of the network
plot_atoms(sys.atoms, fig=fig, ax=ax, dfo=2, alpha=0.1)
plot_verts(sys.net.verts, fig=fig, ax=ax, colors=['r' for i in range(8)])
plot_edges(sys.net.edges, fig=fig, ax=ax)
plot_surfs(sys.net.surfs, simps=True, fig=fig, ax=ax, Show=True, dfo=10)

# Analysis checks:
sys.analyze()
print(sys.atoms[0].cell_vol)
print(sys.net.surfs[0].sa)

# Export the file
# sys.export()