from System.system import System, Atom
from Visualize.mpl_visualize import *

# Create the outer atoms
d, r = 30, 1
outer_atoms_locs = [[-d, 0, 0], [d, 0, 0], [0, d, 0], [0, -d, 0], [0, 0, -d], [0, 0, d]]
outer_atoms = [Atom(_, r) for _ in outer_atoms_locs]

# Case 1: 1 surface, 2 edges

# ATOM     16  O6  MOL     1      19.710  19.200  11.260  1.00  0.00           O
# ATOM    256  OW  SOL   192      18.970  17.660  13.210  1.00  0.00           O
# ATOM    257  HW1 SOL   192      19.210  17.810  12.290  1.00  0.00           H
# ATOM    333  HW2 SOL   254      19.210  19.680  13.760  1.00  0.00           H

c1_atoms = [Atom([-1.01, -0.001, 0], 0.95), Atom([1, 0, 0], 0.95), Atom([0, 2, 0], 1.5), Atom([0, -2, 0], 1.5)]
# c1_atoms = [Atom([19.710, 19.200, 11.260], 1.5), Atom([18.970, 17.660, 13.210], 1.5), Atom([19.210, 17.810, 12.290], 1.3), Atom([19.210, 19.680, 13.760], 1.3)]
c1_sys = System(c1_atoms + outer_atoms)
c1_sys.build_network(max_vert=100, box_size=100, surf_res=0.5)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

plot_balls(c1_sys.atoms, fig=fig, ax=ax)
plot_surfs(c1_sys.net.surfs, fig=fig, ax=ax, simps=True)

plot_verts(c1_sys.net.verts, fig=fig, ax=ax)


plt.show()
