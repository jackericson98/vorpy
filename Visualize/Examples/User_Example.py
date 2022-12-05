import os
os.chdir("../..")
from System.system import *
from Visualize.mpl_visualize import *

vert_atoms = [Atom([1, 5, 5], 2.5), Atom([1, -5, 5], 2.5), Atom([1, 5, -5], 2.5), Atom([1, -5, -5], 2.5)]


cases = [System(atoms=[Atom([-1, 0, 0], .5), Atom([3, 0, 0], .5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([1, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([0.5, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([0.1, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([0, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([-0.25, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([-0.49, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([-0.5, 0, 0], 0.5)] + vert_atoms),
         System(atoms=[Atom([-1.5, 0, 0], 1.5), Atom([-0.55, 0, 0], 0.5)] + vert_atoms)]

for sys in cases:
    sys.build_network(0.1)
    if sys.net.surfs is not None:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection="3d", xlim=10)
        plot_verts(sys.net.surfs[0].verts, fig=fig, ax=ax)
        plot_atoms(sys.net.surfs[0].atoms, fig=fig, ax=ax, colors=['r', 'b'], alpha=.5)
        plot_edges(sys.net.surfs[0].edges, fig=fig, ax=ax)
        plot_surfs(sys.net.surfs[:1], simps=True, fig=fig, ax=ax, dfo=1.5)