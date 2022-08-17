import os
from System.system import System
from Visualize.visualize import plot_edges, plot_verts, plot_atoms, plot_surfs
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = "./Data/test_data/Complex1.pdb"

# Get the System
sys = System(m_file)

sys.build_network(.1, surfs=True)

# Plot the System
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_atoms(sys.atoms[:1], fig=fig, ax=ax, alpha=.1, colors=['w' for i in range(len(sys.atoms))])
plot_verts(sys.atoms[0].verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
plot_edges(sys.atoms[0].edges, fig=fig, ax=ax)
plot_surfs(sys.atoms[0].surfs, fig=fig, ax=ax, alpha=1, simps=True)
print(sys.info)
sys.analyze()
sys.export()
