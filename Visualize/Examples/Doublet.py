from System.system import System, Atom
from Visualize.visualize import *

# Create the outer atoms
d, r = 30, 1
outer_atoms_locs = [[-d, 0, 0], [d, 0, 0], [0, d, 0], [0, -d, 0], [0, 0, -d], [0, 0, d]]
outer_atoms = [Atom(_, r) for _ in outer_atoms_locs]

# Case 1: 1 surface, 2 edges

c1_atoms = [Atom([-1.001, 0.001, 0], 0.95), Atom([1, 0, 0], 0.95), Atom([0, 2, 0], 1.5), Atom([0, -2, 0], 1.5)]
c1_sys = System(c1_atoms + outer_atoms)
c1_sys.build_network(output=False, max_vert=100, box_size=1.5, min_dist=0.5)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

plot_atoms(c1_sys.atoms, fig=fig, ax=ax)
plot_surfs(c1_sys.net.surfs, fig=fig, ax=ax, simps=True)

plot_verts(c1_sys.net.verts, fig=fig, ax=ax)


plt.show()