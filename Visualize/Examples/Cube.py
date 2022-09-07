import os

from System.system import System
from Visualize.visualize import  plot_verts, plot_surfs, plot_edges
import matplotlib.pyplot as plt
os.chdir("../..")
# Create atom objects from sets of points
atoms = [[[0, 0, 0], 1.7]]

dist = 2
rad = 1.2

atoms += [[[dist, 0, 0], rad], [[-dist - 1, 0, 0], rad], [[0, dist, 0], rad], [[0, -dist, 0], rad], [[0, 0, dist], rad],
          [[0, 0, -dist], rad]]

sys = System(atoms)

# Build the surfaces
sys.build_network(get_verts=True)


# Analysis checks:
sys.name = "Cube_03"
# Export the system
# sys.export()
##################################################### Set up the plot ##################################################

fig = plt.figure()


ax = fig.add_subplot(projection="3d")
# Plot the elements of the network
# plot_atoms(sys.atoms, fig=fig, ax=ax, dfo=2, alpha=0.5)
plot_verts(sys.net.verts, fig=fig, ax=ax, colors=['r' for i in range(8)])
plot_edges(sys.net.edges, fig=fig, ax=ax)
plot_surfs(sys.net.surfs, simps=True, fig=fig, ax=ax, dfo=10)

plt.show()

sys.export(export_all=True)
